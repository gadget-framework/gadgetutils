#' Build output summary suitable for feeding into gadgetplots
#'
#' @param model A G3 model built using \code{\link[gadget3]{g3_to_r}} or \code{\link[gadget3]{g3_to_tmb}}. The model must include \code{\link[gadget3]{g3a_report_detail}}.
#' @param params The fitted parameters post-optimisation
#' @param stock_rec_step Named list where each element contains the step(s) for a specific stock (name of the element). E.g. list(bli_imm = 1, bli_dummy = 2). If NULL or a stock is not included/misspelled, all steps are returned.  
#' @param steps Which steps to include in the annual output? Vector of int.
#' @param printatstart Should the stock standard be printed at the start of the timestep (1) or the end (0)
#' @return List of tibbles
#' @export
g3_fit <- function(model, 
                   params, 
                   stock_rec_step = NULL, 
                   steps = 1,
                   printatstart = 1){
  
  ## Checks
  stopifnot(is.list(params))
  if (!(printatstart %in% c(0,1))){
    stop("The printatstart argument must be '0' or '1' (class numeric)")
  }
  
  if (inherits(model, "g3_r")) {
    if (is.data.frame(params)) params <- params$value
      if ("report_detail" %in% names(params) && printatstart == 1) {
          params$report_detail <- 1L
          tmp <- attributes(model(params))
          data_env <- environment(model)
      } else {
          tmp <- NULL
      }
  } else if (inherits(model, "g3_cpp")) {
      if (is.data.frame(params) && "report_detail" %in% params$switch && printatstart == 1) {
          params['report_detail', 'value'] <- 1L
          obj_fun <- gadget3::g3_tmb_adfun(model, params, type = 'Fun')
          tmp <- obj_fun$report(gadget3::g3_tmb_par(params))
          data_env <- as.environment(attr(model, 'model_data'))
      } else {
          tmp <- NULL
      }
  } else {
      stop("Unknown model type: ", class(model))
  }

  # If model is missing actions, then add them and try again
  if (is.null(tmp)) {
    
    ## Incorporate reporting action
    re_actions <- c(attr(model, 'actions'),
                    list(gadget3::g3a_report_detail(attr(model, 'actions'))))
    
    ## Mortality and harvest rates are calculated from the stock info. (biomass/numbers)
    ## at the beginning of the timestep and the numbers/biomass consumed in the timestep
    ## Therefore if we want to print the stock info. at the end of the timestep
    ## this needs to be in addition to, rather than instead of, the reporting
    ## at the begging of the timestep
    
    ## Solution is to calculate mortality within gadget3 rather than from the 
    ## reported output from gadget3
    
    if (printatstart == 0){
      re_actions <- c(re_actions,
                      list(gadget3::g3a_report_history(attr(model, 'actions'), 
                                                       var_re = c('__num$', '__wgt$'),
                                                       out_prefix = 'endprint_',
                                                       run_at = 11)))
    }
    
    ## Compile model
    model <- gadget3::g3_to_r(re_actions)
    
    ## Run model
    if (is.data.frame(params)) params <- params$value
    params$report_detail <- 1L
    tmp <- attributes(model(params))
    data_env <- environment(model)
  }

  ## Add run information back to reporting output, if not already there
  if (is.null(tmp$model_params)) tmp$model_params <- params
  if (is.null(tmp$model_data)) tmp$model_data <- data_env

  return(g3_fit_inner(tmp, stock_rec_step = stock_rec_step, steps = steps))
} 

# Generate fit report from reporting output of g3_fit()
g3_fit_inner <- function(tmp,
                         stock_rec_step = NULL,
                         steps = 1) {
  data_env <- tmp$model_data

  ## Returning parameters as they are put in for now
  ## so model(fit$params) will work without any fiddling
  params <- tmp$model_params
  out_params <- params
  
  ## Calculate the step size as a proportion
  step_lengths <- tmp$step_lengths
  step_size <- 1/length(step_lengths)

  # If dend_ or endprint_ are available, then printatstart was 0
  printatstart <- if (any(grepl("^(dend|endprint)_(.+)__wgt$", names(tmp)))) 0 else 1

  ## Stock, predator, and fleet names
  sp_names <- stock_pred_names(tmp)
  fleet_names <- sp_names$preds[!(sp_names$preds %in% sp_names$stocks)]
  
  ##############################################################################
  
  ## Catch distributions
  catchdist <- g3f_catchdistribution(tmp)
  
  
  ## Sparse distributions
  if (any(grepl('^nll_(a|c)sparse_', names(tmp)))){
    
    sp_df <- tibble::tibble(year = NA_real_, step = NA_real_, area = NA_character_,
                            age = NA_integer_, length = NA_real_, obs_mean = NA_real_,
                            obs_stddev = NA_real_, obs_n = NA_real_,
                            model_sum = NA_real_, model_sqsum = NA_real_, 
                            model_mean = NA_real_, model_n = NA_real_)
    
    ## Observed and model search strings
    sp_re_obs <- 'nll_(asparse|csparse)_([A-Za-z]+)_(.+)__(area$|age$|length$|year$|step$|obs_mean$|obs_stddev|obs_n)'
    sp_re <-     'nll_(asparse|csparse)_([A-Za-z]+)_(.+)__(model_sum$|model_sqsum$|model_n$|nll$)'
    sp_re_all <- 'nll_(asparse|csparse)_([A-Za-z]+)_(.+)__(age$|length$|year$|step$|obs_mean$|obs_stddev|obs_n$|model_sum$|model_sqsum$|model_n$|nll$)'
    
    sparse_nll_names <- names(mget(ls(data_env)[grep(sp_re_obs, ls(data_env))], envir = data_env))
    
    sparsedist <- 
      mget(ls(data_env)[grep(sp_re_obs, ls(data_env))], envir = data_env) %>%
      purrr::map(~tibble::tibble(value = as.numeric(.))) %>% 
      dplyr::bind_rows(.id = 'comp') %>%
      dplyr::bind_rows(
        tmp[grep(sp_re, names(tmp))] %>% 
          purrr::map(~tibble::tibble(value = as.numeric(.))) %>% 
          dplyr::bind_rows(.id = 'comp')
      ) %>% 
      dplyr::mutate(data_type = gsub(sp_re_all, '\\1', .data$comp),
                    function_f = gsub(sp_re_all, '\\2', .data$comp),
                    component =  gsub(sp_re_all, '\\3', .data$comp),
                    column = gsub(sp_re_all, '\\4', .data$comp)) %>% 
      #component = gsub('__', '.', .data$component)) %>% 
      dplyr::select(-.data$comp) %>% 
      ## Take weights from parameters
      dplyr::full_join(
        params$value[unique(gsub('__(.+)$', '_weight', gsub('^nll_', '', sparse_nll_names)))] %>% 
          purrr::map(~tibble::tibble(weight = as.numeric(.))) %>% 
          dplyr::bind_rows(.id = 'comp') %>% 
          dplyr::mutate(data_type = gsub('(asparse|csparse)_([A-Za-z]+)_(.+)_weight$', '\\1', .data$comp),
                        function_f = gsub('(asparse|csparse)_([A-Za-z]+)_(.+)_weight$', '\\2', .data$comp),
                        component =  gsub('(asparse|csparse)_([A-Za-z]+)_(.+)_weight$', '\\3', .data$comp)) %>% 
          dplyr::select(-.data$comp)
      , by = c('data_type', 'function_f', 'component'))
    
    sparsedist <- 
      c(
        split(sparsedist, sparsedist$component) %>% 
          purrr::map(function(x){
            x %>% 
              dplyr::group_by(.data$column) %>% 
              dplyr::mutate(row = dplyr::row_number()) %>% 
              tidyr::pivot_wider(names_from = .data$column, values_from = .data$value, values_fill = NA) %>% 
              dplyr::ungroup() %>% 
              dplyr::select(-row)
          }),
        list(sp_df)) %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(model_mean = .data$model_sum / .data$model_n) %>% 
      tidyr::drop_na(.data$component) %>% 
      dplyr::select(.data$year, .data$step, .data$area, .data$data_type, 
                    .data$function_f, .data$component, .data$age, .data$length, 
                    .data$obs_mean, .data$obs_stddev, .data$obs_n, 
                    .data$model_sum, .data$model_mean, .data$model_sqsum, 
                    .data$model_n, .data$weight, .data$nll)
    
    ## Fix lin reg nlls
    if (any(sparsedist$function_f == 'linreg')){
      sparsedist <- 
        do.call('rbind', 
                lapply(
                  split(sparsedist, sparsedist$component),
                  function(x){
                    if (unique(x$function_f) == 'linreg'){
                      x$nll <- 0
                      x$nll[nrow(x)] <- 
                        tmp[grep(paste0('nll_(asparse|csparse)_linreg_', 
                                        unique(x$component), '__(nll$)'), 
                                 names(tmp))][[1]][['nll']] 
                    }
                    return(x)
                  })
                )
    }
    
  }else{
    sparsedist <- NULL
  }
  
  
  ## Survey or other indices 
  sidat <- g3f_sidat(tmp)
  
  ## Suitability
  if (any(grepl('^suit_(.+)_(.+)__report', names(tmp)))){
    suit_re <- paste0(
        "^suit",
        "_(\\Q", paste(sp_names$stocks, collapse = "\\E|\\Q"), "\\E)",
        "_(.+)",
        "__report" )
    
    # New-style g3a_suitability_report() output
    suitability <-
      tmp[grep(suit_re, names(tmp))] %>%
      purrr::map(as.data.frame.table, stringsAsFactors = F, responseName = 'suit') %>%
      dplyr::bind_rows(.id = 'comp') %>%
      dplyr::mutate(stock = gsub(suit_re, '\\1', .data$comp),
                    fleet = gsub(suit_re, '\\2', .data$comp)) %>%
      extract_year_step() %>% 
      dplyr::select(.data$stock, .data$fleet, 
                    dplyr::matches("year|step|area|age|length|predator_length"),
                    .data$suit) %>%
      identity()
    
    # Fix length and age if present
    if ('length' %in% names(suitability)){
      suitability <-
        suitability %>% 
        split_length() %>%
        dplyr::group_by(.data$stock, .data$fleet) %>% 
        dplyr::group_modify(~replace_inf(.x)) %>% 
        dplyr::ungroup() %>%
        dplyr::mutate(length = (.data$upper + .data$lower)/2) %>% 
        dplyr::select(-c(.data$lower, .data$upper))
    }
    
    if ('age' %in% names(suitability)){
      suitability$age <- as.numeric(gsub('age', '', suitability$age))
    }
  } else if (any(grepl('detail_(.+)__suit_', names(tmp)))){

    suitability <-
      tmp[grep('detail_(.+)__suit_', names(tmp))] %>%
      purrr::map(as.data.frame.table, stringsAsFactors = F) %>%
      dplyr::bind_rows(.id = 'comp') %>%
      dplyr::mutate(stock = gsub('detail_(.+)__suit_(.+)$', '\\1', .data$comp),
                    fleet = gsub('detail_(.+)__suit_(.+)$', '\\2', .data$comp)) %>% 
                    #area = as.numeric(.data$area)) %>%
                    #length = gsub('len','', .data$length) %>% as.numeric(),
                    #age = gsub('age','', .data$age) %>% as.numeric()) %>%
      split_length() %>%
      dplyr::group_by(.data$stock, .data$fleet) %>% 
      dplyr::group_modify(~replace_inf(.x)) %>% 
      dplyr::ungroup() %>%
      extract_year_step() %>%
      dplyr::rename(suit = .data$Freq) %>%
      dplyr::select(.data$stock, .data$fleet, 
                    dplyr::matches("year|step|area|age|length|predator_length"),
                    .data$suit) %>%
      tibble::as_tibble()

  }else{
    warning("Cannot find suitability tables, no suitability data will be included in results")
    suitability <- data.frame(
      stock = "x",
      fleet = "x",
      year = 1,
      step = 1,
      area = "a",
      age = 1,
      length = 1,
      predator_length = 1,
      suit = 1 )[c(),,drop = FALSE]
    }
  
  ## --------------------------------------------------------------------------------
  
  ## Likelihood
  if (any(grepl('^nll_', names(tmp)))){
    
    likelihood <-
      ## Catch and abundance distributions values
      tmp[grep('^nll_(a|c)dist_(.+)__(wgt$|num$)', names(tmp))] %>%
      purrr::map(~tibble::tibble(time=names(.), value = as.numeric(.))) %>% 
      dplyr::bind_rows(.id='lik_comp') %>% 
      dplyr::mutate(measurement = gsub('.+(wgt|num)','\\1', .data$lik_comp),
                    component = gsub('nll_(cdist|adist)_([A-Za-z]+)_(.+)__(wgt$|num$|weight$)', '\\3', .data$lik_comp),
                    data_type = gsub('nll_(cdist|adist)_([A-Za-z]+)_(.+)__(wgt$|num$|weight$)', '\\1_\\2', .data$lik_comp)) %>% 
      dplyr::select(-.data$lik_comp) %>% 
      dplyr::left_join(
        tmp[grep('^nll_(a|c)dist_(.+)__weight$', names(tmp))] %>%
          purrr::map(~tibble::tibble(time=names(.), weight = as.numeric(.))) %>% 
          dplyr::bind_rows(.id='lik_comp') %>% 
          dplyr::mutate(component = gsub('nll_(cdist|adist)_([A-Za-z]+)_(.+)__(wgt$|num$|weight$)', '\\3', .data$lik_comp),
                        data_type = gsub('nll_(cdist|adist)_([A-Za-z]+)_(.+)__(wgt$|num$|weight$)', '\\1_\\2', .data$lik_comp)) %>% 
          dplyr::select(-.data$lik_comp),
      by = c('component', 'data_type', 'time')) %>% 
      ## Understocking
      dplyr::bind_rows(
        tmp[grep('^nll_understocking', names(tmp))] %>%
          purrr::map(~tibble::tibble(time=names(.), value = as.numeric(.))) %>% 
          dplyr::bind_rows(.id='lik_comp') %>% 
          dplyr::mutate(component = 'understocking',
                        data_type = 'model_preystocks',
                        measurement = 'wgt') %>% 
          dplyr::select(-.data$lik_comp) 
      ) %>% 
      extract_year_step() %>% 
      dplyr::select(.data$year, .data$step, .data$component, .data$data_type, 
                    .data$measurement, .data$weight, .data$value)
    
    ## Sparse data
    if (!is.null(sparsedist)){
      likelihood <-
        likelihood %>% 
        dplyr::bind_rows(
          sparsedist %>% 
            dplyr::group_by(.data$year, .data$step, .data$component, 
                            .data$function_f, .data$data_type, .data$weight) %>% 
            dplyr::summarise(value = sum(.data$nll), .groups = 'drop') %>% 
            dplyr::mutate(data_type = paste(.data$data_type, 
                                            .data$function_f, sep = '_'),
                          measurement = '') %>% 
            dplyr::select(.data$year, .data$step, .data$measurement, 
                          .data$component, .data$data_type, .data$weight, 
                          .data$value)  
        )
      sparsedist <- sparsedist %>% dplyr::select(-c(.data$weight, .data$nll))
    }
  }else{
    likelihood <- NULL
  }
  
  ## ---------------------------------------------------------------------------
  ## Stock-recruitment
  ## ---------------------------------------------------------------------------
  
  stock.recruitment <- g3f_stock.recruitment(tmp, stock_rec_step)
  
  ## ---------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  
  ## Stock reports
  weight_reports <- 
    tmp[grepl("^(dstart|detail)_(.+)__wgt$", names(tmp))] %>% 
    purrr::map(as.data.frame.table, stringsAsFactors = FALSE, responseName = 'weight') %>% 
    dplyr::bind_rows(.id='comp') %>% 
    dplyr::mutate(stock = gsub("^(dstart|detail)_(.+)__wgt$", '\\2', .data$comp)) %>% 
    dplyr::select(-.data$comp) 
  
  if (any(grepl("__cons$", names(tmp)))) {
    # detail_(prey)_(pred)__cons arrays
    consumption_re <- paste0(
        "^detail",
        "_(\\Q", paste(sp_names$stocks, collapse = "\\E|\\Q"), "\\E)",
        "_(.+)",
        "__cons" )
  } else {
    # Pre-predation stype __predby_ arrays
    consumption_re <- paste0(
        "^detail",
        "_(\\Q", paste(sp_names$stocks, collapse = "\\E|\\Q"), "\\E)",
        "__predby",
        "_(.+)$" )
  }

  fleet_reports <- 
    tmp[grepl(consumption_re, names(tmp))] %>%
    purrr::map(as.data.frame.table, stringsAsFactors = FALSE, responseName = 'biomass_consumed') %>% 
    dplyr::bind_rows(.id='comp') %>% 
    dplyr::mutate(stock = gsub(consumption_re, '\\1', .data$comp),
                  fleet = gsub(consumption_re, '\\2', .data$comp)) %>%
    dplyr::select(-.data$comp) %>% 
    dplyr::left_join(weight_reports, by = c("time", "area", "stock", "age", "length")) %>% 
    split_length() %>% 
    dplyr::group_by(.data$stock, .data$fleet) %>% 
    dplyr::group_modify(~replace_inf(.x)) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(avg.length = (.data$lower + .data$upper)/2) %>% 
    dplyr::mutate(number_consumed = 
                    ifelse(.data$biomass_consumed == 0, 0, .data$biomass_consumed / .data$weight)) %>%
    extract_year_step() %>% 
    tibble::as_tibble()
  
  ## Abundance
  num_reports <- 
    tmp[grepl("^(dstart|detail)_(.+)__num$", names(tmp))] %>% 
    purrr::map(as.data.frame.table, stringsAsFactors = FALSE, responseName = 'abundance') %>% 
    dplyr::bind_rows(.id='comp') %>% 
    dplyr::mutate(stock = gsub("^(dstart|detail)_(.+)__num$", '\\2', .data$comp)) %>% 
    dplyr::select(-.data$comp) %>% 
    dplyr::left_join(weight_reports, by = c("time", "area", "stock", "age", "length")) %>% 
    split_length() %>%
    dplyr::group_by(.data$stock) %>% 
    dplyr::group_modify(~replace_inf(.x)) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(avg.length = (.data$lower + .data$upper)/2) %>% 
    extract_year_step() %>% 
    tibble::as_tibble()
  
  ## Stock full
  stock.full <- 
    num_reports %>%
    dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$avg.length) %>% 
    dplyr::summarise(number = sum(.data$abundance), 
                     mean_weight = sum(.data$abundance*.data$weight)/sum(.data$abundance), .groups = 'drop') %>% 
    tidyr::replace_na(list(mean_weight = 0)) %>% 
    dplyr::rename(length = .data$avg.length)
  
  ## Stock std
  stock.std <- 
    num_reports %>% 
    dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$age) %>% 
    dplyr::summarise(number = sum(.data$abundance),
                     mean_length = sum(.data$avg.length*.data$abundance)/sum(.data$abundance),
                     stddev_length = sqrt(sum((.data$avg.length-.data$mean_length)^2*.data$abundance)/sum(.data$abundance)),
                     mean_weight = sum(.data$abundance*.data$weight)/sum(.data$abundance), .groups = 'drop') %>% 
    dplyr::mutate(age = gsub('age', '', .data$age) %>% as.numeric())
  
  ## Stock prey - fleet consumption first
  stock.prey <- 
    fleet_reports %>%
    dplyr::filter(.data$fleet %in% fleet_names) %>% 
    dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$age) %>%
    dplyr::summarise(number_consumed = sum(.data$number_consumed),
                     biomass_consumed = sum(.data$biomass_consumed), .groups = 'drop') %>% 
    dplyr::mutate(age = gsub('age', '', .data$age) %>% as.numeric()) %>% 
    dplyr::left_join(stock.std %>% 
                       dplyr::select(.data$year, .data$step, .data$area, .data$stock, .data$age, .data$number),
                     by = c("year", "step", "area", "stock", "age")) %>% 
    dplyr::mutate(mortality = -log(1 - .data$number_consumed / .data$number)/step_size) %>% 
    tidyr::replace_na(list(mortality = 0)) 
  
  ## Add in biological consumption if included in model
  if (length(fleet_names) < length(sp_names$preds)){
    
    stock.prey <- 
      stock.prey %>% 
      dplyr::left_join(
        fleet_reports %>% 
          dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$age) %>%
          dplyr::summarise(number_consumed_total = sum(.data$number_consumed),
                           biomass_consumed_total = sum(.data$biomass_consumed), .groups = 'drop') %>% 
          dplyr::mutate(age = gsub('age', '', .data$age) %>% as.numeric()),
      by = c('year','step','area','stock','age')) %>% 
      dplyr::mutate(mortality_total = -log(1 - .data$number_consumed_total / .data$number)/step_size) %>% 
      tidyr::replace_na(list(mortality_total = 0))
    
  }
  
  ## -------------------------------
  ## Predation - fleets & biological
  ## -------------------------------
  
  ## Predator - prey
  if (length(sp_names$preds) > 0){
    
    # predator.prey <- 
    #   fleet_reports %>%
    #   dplyr::left_join(suitability %>% dplyr::rename(avg.length = length)) %>%
    #   dplyr::left_join(num_reports %>% dplyr::select(-c(.data$upper,.data$lower,.data$avg.length)) %>% 
    #                      dplyr::rename(bioweight = .data$weight), by = c('year', 'step', 'stock', 'area', 'age', 'length')) %>% 
    #   dplyr::mutate(suit = dplyr::case_when(number_consumed == 0 ~ 0, .default = .data$suit)) %>% 
    #   dplyr::mutate(length = .data$avg.length,
    #                 harv.bio = .data$abundance * .data$bioweight * .data$suit) %>% 
    #   dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$fleet, .data$length) %>%
    #   dplyr::summarise(number_consumed = sum(.data$number_consumed),
    #                    biomass_consumed = sum(.data$biomass_consumed),
    #                    harv.bio = sum(.data$harv.bio)) %>% 
    #   dplyr::left_join(stock.full, by = c("year", "step", "area", "stock", "length")) %>% 
    #   dplyr::mutate(mortality = -log(1 - .data$number_consumed / .data$number)/step_size) %>% 
    #   tidyr::replace_na(list(mortality = 0)) %>% 
    #   dplyr::ungroup() %>% 
    #   dplyr::mutate(tb = .data$number * .data$mean_weight) %>% 
    #   dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$fleet) %>%
    #   dplyr::mutate(suit = .data$harv.bio / .data$tb,
    #                 suit = ifelse(is.finite(.data$suit), .data$suit, 0)) %>%
    #   dplyr::rename(predator = .data$fleet, prey = .data$stock) %>% 
    #   dplyr::select(-c(.data$number, .data$mean_weight)) %>% 
    #   dplyr::ungroup()  
    
    not_all_na <- function(x) any(!is.na(x))
   
    predator.prey <- 
      lapply(split(fleet_reports, list(fleet_reports$stock, fleet_reports$fleet)), function(x, suits){
      if (nrow(x) == 0) return(NULL)
        suppressMessages(
          x %>% 
            dplyr::left_join(
              suits %>% 
                dplyr::rename(avg.length = length) %>% 
                dplyr::filter(.data$stock %in% unique(x$stock),
                              .data$fleet %in% unique(x$fleet)) %>% 
                dplyr::select(dplyr::where(not_all_na))
              )  
          )
      }, suits = suitability) %>% 
      dplyr::bind_rows() %>% 
      dplyr::left_join(num_reports %>% dplyr::select(-c(.data$upper,.data$lower,.data$avg.length)) %>% 
                         dplyr::rename(bioweight = .data$weight), by = c('year', 'step', 'stock', 'area', 'age', 'length')) %>% 
      dplyr::mutate(suit = dplyr::case_when(number_consumed == 0 ~ 0, .default = .data$suit)) %>% 
      dplyr::mutate(length = .data$avg.length,
                    harv.bio = .data$abundance * .data$bioweight * .data$suit) %>% 
      dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$fleet, .data$length) %>%
      dplyr::summarise(number_consumed = sum(.data$number_consumed),
                       biomass_consumed = sum(.data$biomass_consumed),
                       harv.bio = sum(.data$harv.bio)) %>% 
      dplyr::left_join(stock.full, by = c("year", "step", "area", "stock", "length")) %>% 
      dplyr::mutate(mortality = -log(1 - .data$number_consumed / .data$number)/step_size) %>% 
      tidyr::replace_na(list(mortality = 0)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(tb = .data$number * .data$mean_weight) %>% 
      dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$fleet) %>%
      dplyr::mutate(suit = .data$harv.bio / .data$tb,
                    suit = ifelse(is.finite(.data$suit), .data$suit, 0)) %>%
      dplyr::rename(predator = .data$fleet, prey = .data$stock) %>% 
      dplyr::select(-c(.data$number, .data$mean_weight)) %>% 
      dplyr::ungroup() 
    
   
    ## Fleet info
    fleet.info <- 
      predator.prey %>%
      dplyr::select(.data$year,
                    .data$step, 
                    .data$area,
                    fleet = .data$predator, 
                    stock = .data$prey, 
                    .data$length, 
                    .data$harv.bio,
                    .data$biomass_consumed) %>% 
      dplyr::group_by(.data$year, .data$step, .data$area, .data$fleet) %>%
      dplyr::summarise(harv.bio = sum(.data$harv.bio),
                       amount = sum(.data$biomass_consumed), .groups = 'drop')  %>%
      dplyr::group_by(.data$year, .data$step, .data$area, .data$fleet) %>%
      dplyr::mutate(amount = ifelse(is.na(.data$amount), 0, .data$amount),
                    harv.rate = ifelse(.data$amount == 0, 0, .data$amount / .data$harv.bio)) %>% 
      dplyr::ungroup()
    
    ## -------------------------------
    ## Predation - fleets
    ## -------------------------------
    
    if (length(fleet_names) > 0){
      
      ## Suitability over fleets                 
      harv.suit <- 
        predator.prey %>%
        dplyr::filter(.data$predator %in% fleet_names) %>% 
        dplyr::filter(.data$biomass_consumed > 0) %>%
        dplyr::group_by(.data$year, .data$step, .data$prey, .data$length) %>%
        dplyr::summarise(suit = sum(.data$biomass_consumed * .data$suit) / sum(.data$biomass_consumed), .groups = 'drop') %>%
        dplyr::rename(stock = .data$prey)
      
      ## TO-DO: Add weighted mortality 
      f.by.year <-
        stock.prey %>% 
        dplyr::left_join(
          stock.prey %>% 
            dplyr::group_by(.data$stock) %>% 
            dplyr::summarise(age.min = max(.data$age), age.max=max(.data$age), .groups = 'drop'), 
          by = "stock") %>% 
        dplyr::group_by(.data$stock, .data$year, .data$area) %>%
        dplyr::summarise(catch = sum(.data$biomass_consumed),
                         num.catch = sum(.data$number_consumed),
                         F = mean(.data$mortality[.data$age >= .data$age.min & .data$age <= .data$age.max]), .groups = 'drop')
      
    }else{
      harv.suit <- stock.full[,c('year','step','stock','length')] %>% unique.data.frame() %>% transform(suit = NA)
      f.by.year <- stock.full[,c('year','area','stock')] %>% unique.data.frame() %>% transform(F = NA)
    }
  }else{
    fleet.info <- NULL
  }
  
  ## Re-calculate stock.std and stock.full if we are printing at the end of the timestep
  if (printatstart == 0){
    
    ## Stock reports
    weight_reports <- 
      tmp[grepl('endprint_(.+)__wgt', names(tmp))] %>% 
      purrr::map(as.data.frame.table, stringsAsFactors = FALSE, responseName = 'weight') %>% 
      dplyr::bind_rows(.id='comp') %>% 
      dplyr::mutate(stock = gsub('endprint_(.+)__wgt', '\\1', .data$comp)) %>% 
      dplyr::select(-.data$comp) 
    
    ## Abundance
    num_reports <- 
      tmp[grepl('endprint_(.+)__num', names(tmp))] %>% 
      purrr::map(as.data.frame.table, stringsAsFactors = FALSE, responseName = 'abundance') %>% 
      dplyr::bind_rows(.id='comp') %>% 
      dplyr::mutate(stock = gsub('endprint_(.+)__num', '\\1', .data$comp)) %>% 
      dplyr::select(-.data$comp) %>% 
      dplyr::left_join(weight_reports, by = c("time", "area", "stock", "age", "length")) %>% 
      split_length() %>%
      dplyr::group_by(.data$stock) %>% 
      dplyr::group_modify(~replace_inf(.x)) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(avg.length = (.data$lower + .data$upper)/2) %>% 
      extract_year_step() %>% 
      tibble::as_tibble()
    
    ## Stock full
    stock.full <- 
      num_reports %>%
      dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$avg.length) %>% 
      dplyr::summarise(number = sum(.data$abundance), 
                       mean_weight = mean(.data$weight)) %>% 
      dplyr::ungroup() %>% 
      dplyr::rename(length = .data$avg.length)
    
    ## Stock std
    stock.std <- 
      num_reports %>% 
      dplyr::group_by(.data$year, .data$step, .data$area, .data$stock, .data$age) %>% 
      dplyr::summarise(number = sum(.data$abundance),
                       mean_length = sum(.data$avg.length*.data$abundance)/sum(.data$abundance),
                       stddev_length = sqrt(sum((.data$avg.length-.data$mean_length)^2*.data$abundance)/sum(.data$abundance)),
                       mean_weight = sum(.data$abundance*.data$weight)/sum(.data$abundance)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(age = gsub('age', '', .data$age) %>% as.numeric())
    
  }
  
  ## Annual output
  res.by.year <- 
    stock.full %>%
    dplyr::filter(.data$step %in% steps) %>%
    dplyr::left_join(harv.suit,
                     by = c('stock', 'year', 'step', 'length')) %>%
    dplyr::group_by(.data$stock, .data$year, .data$area, .data$step) %>%
    dplyr::summarise(total.number = sum(.data$number),
                     total.biomass = sum(.data$number * .data$mean_weight),
                     harv.biomass = sum(.data$number * .data$suit * .data$mean_weight, na.rm = TRUE), .groups = 'drop') %>%
    dplyr::left_join(f.by.year,
                     by = c('stock', 'year', 'area')) %>% 
    dplyr::left_join(stock.recruitment %>% 
                       dplyr::group_by(.data$stock, .data$year, .data$area) %>% 
                       dplyr::summarise(recruitment = sum(.data$recruitment), .groups = 'drop'),
                     by = c('stock', 'year', 'area')) %>% 
    dplyr::ungroup()
  
  out <- list(
    sidat = sidat,
    stockdist = catchdist$stockdist,
    catchdist.fleets = catchdist$catchdist.fleets,
    sparsedist = sparsedist,
    #suitability = predator.prey %>% 
    #  dplyr::select(.data$year,.data$step,area=.data$area, stock=.data$prey,
    #                fleet=.data$predator,.data$length,.data$suit) %>% 
    #  dplyr::ungroup(),
    suitability = suitability %>% 
      dplyr::bind_rows(
        tibble::tibble(stock = NA_character_, fleet = NA_character_,
                       year = NA_real_, step = NA_real_, area = NA_character_,
                       age = NA_integer_, length = NA_real_, predator_length = NA_character_)
        ) %>% utils::head(-1),
    likelihood = likelihood,
    stock.prey = stock.prey,
    stock.std = stock.std,
    stock.full = stock.full,
    fleet.info = fleet.info,
    predator.prey = predator.prey,
    stock.recruitment = stock.recruitment,
    res.by.year = res.by.year,
    params = out_params,
    score = data.frame(nll = tmp$nll)
  )
  
  ## Add 'summary' attribute if it exists
  if ('summary' %in% names(attributes(out_params))){
    attributes(out)$summary <- attr(out_params, 'summary')
  }
  
  class(out) <- c('gadget.fit',class(out))
  return(out)
  
}

#' @title Extracts year and step columns
#' @description Gadget3 reports have a time column in the format 'YEAR-STEP'. This function creates a 'year' and 'step' column by splitting the existing 'time' column. This is useful for processing gagdet3 reports.
#' @param data dataframe of gadget3 reports. Must include a 'time' column.
#' @export
extract_year_step <- function(data){
  if (!("time" %in% names(data))) return(data)
  
  if (any(!grepl('-', data$time))){
    data[!grepl('-', data$time), 'time'][[1]] <- 
      paste0(data[!grepl('-', data$time), 'time'][[1]], '-0')
  }
  
  data <- 
    within(data, {
    step = gsub('^(\\d+)-(\\d+)$', "\\2", time) |> as.numeric();
    year = gsub('^(\\d+)-(\\d+)$', "\\1", time) |> as.numeric()
  })
  
  return(data[,names(data) != "time"])
  
}

## -----------------------------------------------------------------------------
## Internal functions used by g3_fit
## -----------------------------------------------------------------------------

split_age <- function(data){
  if (!('age' %in% names(data))) return(data)
  
  tmp <-
    data %>% 
    dplyr::mutate(lower_age = gsub('(.+):(.+)', '\\1', .data$age),
                  upper_age = gsub('(.+):(.+)', '\\2', .data$age)) 
  return(tmp)
}

split_length <- function(data){
  
  if (!('length' %in% names(data))) return(data)
  
  data <- 
    data |> 
    within( {
      upper = gsub('(.+):(.+)', '\\2', length) |> as.numeric();
      lower = gsub('(.+):(.+)', '\\1', length) |> as.numeric()
    } )
  
  return(data)
  
}

replace_inf <- function(data, group_col = NULL){
  
  if (!('upper' %in% names(data))) return(data)
  
  inner_fun <- function(x){
    ind <- is.infinite(x$upper)
    if (any(ind)){
      replacement <- abs(diff(sort(unique(x$upper[!ind]), decreasing = TRUE)[1:2]))
      if (!is.na(replacement)) x$upper[ind] <- x$lower[ind] + replacement
    } 
    return(x)
  }
  
  if (is.null(group_col)) return(inner_fun(data))
  else {
    return(
      do.call("rbind", lapply(split(data, data[[group_col]]), FUN = inner_fun))
    )
  }
  
  return(data)
}

## Need to check what happens if no predation
stock_pred_names <- function(reports, stock_names = NULL){
  
  inner_fun <- function(reports, re, pos){
    
    tmp_names <- regmatches(names(reports), regexec(re, names(reports)))
    tmp_names <- vapply(tmp_names, function(x){
      if (length(x) == pos) x[[pos]] else as.character(NA)
    }, as.character(1))
    return(
      unique(tmp_names[!is.na(tmp_names)])
    )
    
  }
  
  if (is.null(stock_names)) 
    stock_names <- inner_fun(reports, "^(?:dstart|detail)_(.+)__num$", 2)
  
  pred_names <- 
    inner_fun(reports,
              paste0("^detail",
                     "_(\\Q", paste(stock_names, collapse = "\\E|\\Q"), "\\E)",
                     "_(.+)",
                     "__cons"), 3)
  if (length(pred_names) == 0) pred_names <- NA_character_
  
  return(list(stocks = stock_names, preds = pred_names))
  
}

## Helper to add columns
add_missing_columns <- function(data, columns){
  
  for (i in names(columns)){
    if (!i %in% names(data))  data[[i]] <- columns[[i]]
  }
  return(data)
}
