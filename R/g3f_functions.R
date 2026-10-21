#' Summarises the fit of a gadget3 model to catch distributions
#'
#' @param reports Reported output of a G3 model
#' @param by_stock Should the stock distribution be compiled?
#' @param by_fleet Should the catch distribution per fleet be compiled?
#' @param reports Reported output of a G3 model
#' @return list containing two data frames: catchdist.fleets and stockdist
#' @export
g3f_catchdistribution <- function(reports, by_stock = TRUE, by_fleet = TRUE){
  
  cd_str <- '^cdist_([A-Za-z]+)_(.+)_(model|obs)__(num|wgt)$'
  
  if (any(grepl(cd_str, names(reports)))){
    
    cd_data <- 
      do.call("rbind", 
              lapply(stats::setNames(names(reports)[grep(cd_str, names(reports))],
                                     names(reports)[grep(cd_str, names(reports))]),
                     function(x){
                       out <- as.data.frame.table(reports[[x]], stringsAsFactors = FALSE)
                       out$report_name <- x
                       if (!"stock" %in% names(out)) out$stock <- "-999"
                       if (!"stock_re" %in% names(out)) out$stock_re <- "-999"
                       if (!"predator_length" %in% names(out)) out$predator_length <- "-999"
                       return(out)
                      }
                     )
              )
    
    cd_data <- 
      within(cd_data, {
        data_function = gsub(cd_str, '\\1', report_name);
        origin = gsub(cd_str, '\\3', report_name);
        name = gsub(cd_str, '\\2', report_name)})
    
    ## Pivot predictions wider and merge in parameters
    cd_data <- 
      stats::reshape(
        cd_data, 
        timevar = "origin",
        idvar = c("time", "area", "age", "length",
                  "stock", "stock_re", "predator_length",
                  "name", "data_function"),
        direction = "wide") 
    
    ## Inner function to process arrays
    catchdist_inner <- function(data){
      
      if (all(data$stock == "-999" & data$stock_re == "-999"))
        cols <- c("time", "area", "name")
      else
        cols <- c("time", "area", "name", "length", "age", "predator_length")
      
      data <- 
        merge(x = data,
              y = 
                data |> 
                stats::aggregate(
                  as.formula(paste("cbind(Freq.model, Freq.obs) ~ ", 
                                   paste(cols, collapse = "+"))),
                  FUN = sum,
                  na.rm = TRUE),
              by = cols)
      
      data$predicted <- data$Freq.model.x / data$Freq.model.y
      data$observed <- data$Freq.obs.x / data$Freq.obs.y
      
      data <-
        data |> 
        split_length() |> 
        replace_inf(group_col = c("stock")) |> 
        within(data, {
          stock = gsub("-999", NA_character_, stock);
          stock_re = gsub("-999", NA_character_, stock_re);
          predator_length = gsub("-999", NA_character_, predator_length)
          #residuals = ifelse(observed == 0, NA_real_, observed - predicted)
        }) |> 
        extract_year_step() |> 
        identity()
      
      ## Noting working inside "within" for some reason...
      data$residuals <- ifelse(data$observed == 0, 
                               NA_real_, 
                               data$observed - data$predicted)
      
      if (all(c("upper", "lower") %in% names(data))){
        data$avg.length <- (data$upper + data$lower)/2
      }
      
      names(data)[match(c("Freq.model.x", "Freq.obs.x",
                          "Freq.model.y", "Freq.obs.y"), names(data))] <- 
        c("pred.num", "obs.num",
          "total.pred", "total.catch")
      
      return(data[order(data$name, data$year, 
                        data$step, data$stock, 
                        data$avg.length),])
      
    }
    
    ## Catch distributions by stock
    nastock_index <- cd_data$stock == "-999" & cd_data$stock_re == "-999"
    
    if (all(nastock_index) | !by_stock){
      stockdist <- NULL
    } else {
      stockdist <- catchdist_inner(cd_data[!nastock_index,])
      
      names(stockdist)[match(c("predicted", "observed", "avg.length"), 
                             names(stockdist))] <- 
        c("pred.ratio", "obs.ratio", "length")
      
      stockdist <- stockdist[,c("name", "year", "step", "area", 
                                "predator_length", "stock", "stock_re", 
                                "lower", "upper", "length", "age", 
                                "obs.num", "obs.ratio", "pred.num", "pred.ratio")]
      
    }
    
    ## Catch distributions by fleet
    if (all(!nastock_index) | !by_fleet){
      catchdist.fleets <- NULL
    } else {
      catchdist.fleets <- catchdist_inner(cd_data[nastock_index,])
      
      catchdist.fleets <- catchdist.fleets[,c("name", "year", "step", "area", 
                                              "stock", "stock_re", "lower", 
                                              "upper", "avg.length", "age", 
                                              "obs.num", "total.catch", 
                                              "observed", "pred.num",
                                              "total.pred", "predicted", 
                                              "residuals")]
    }
    
    } else {
      stockdist <- NULL
      catchdist.fleets <- NULL
    }
  
  return(list(catchdist.fleets = catchdist.fleets, stockdist = stockdist))
         
  
}

#' Summarises the fit of a gadget3 model to distribution data from landings/samples
#'
#' @param reports Reported output of a G3 model
#' @return Data frame or NULL
#' @export
g3f_catchdist.fleets <- function(reports) g3f_catchdistribution(reports, by_stock = FALSE)$catchdist.fleets

#' Summarises the fit of a gadget3 model to stock distribution taken from landings/samples
#'
#' @param reports Reported output of a G3 model
#' @return Data frame or NULL
#' @export
g3f_stockdist <- function(reports) g3f_catchdistribution(reports, by_fleet = FALSE)$stockdist

#' Summarises the fit of a gadget3 model to sparse data
#'
#' @param reports Reported output of a G3 model
#' @return Data frame or NULL
#' @export
g3f_sparsedata <- function(reports, data_env){
  
  ## Sparse distributions
  if (any(grepl('^nll_(a|c)sparse_', names(reports)))){
    
    ## Observed and model search strings
    sp_re_obs <- 'nll_(asparse|csparse)_([A-Za-z]+)_(.+)__(area|age|length|year|step|obs_mean|obs_stddev|obs_n)$'
    sp_re <-     'nll_(asparse|csparse)_([A-Za-z]+)_(.+)__(model_sum|model_sqsum|model_n)$'
    sp_re_all <- 'nll_(asparse|csparse)_([A-Za-z]+)_(.+)__(age|length|year|step|obs_mean|obs_stddev|obs_n|model_sum|model_sqsum|model_n)$'
    
    ## Observations taken from model environment, model output taken from reports
    sparse_reports <- 
      c(mget(ls(data_env)[grep(sp_re_obs, ls(data_env))], envir = data_env),
        reports[grep(sp_re, names(reports))])
    
    sparsedist <- 
      do.call("rbind",
              lapply(stats::setNames(names(sparse_reports),
                                     names(sparse_reports)), function(x){
                                       out <- data.frame(value = sparse_reports[[x]])
                                       out$column <- gsub(sp_re_all, "\\4", x)
                                       out$id <- gsub("__(.+)$", "", x)
                                       return(out)
                                     }))
    
    
    sparsedist <- 
      do.call("rbind",
              lapply(split(sparsedist, sparsedist$id), function(x){
                x_re <- "nll_(asparse|csparse)_([A-Za-z]+)_(.+)$"
                out <- 
                  utils::unstack(x, "value ~ column") |> 
                  within( {
                    component = gsub(x_re, '\\3', unique(x$id));
                    function_f = gsub(x_re, '\\2', unique(x$id));
                    data_type = gsub(x_re, '\\1', unique(x$id));
                    })
                
                out <- 
                  add_missing_columns(
                    out,
                    list(year = NA_real_, step = NA_real_, area = NA_character_,
                         age = NA_integer_, length = NA_real_,
                         obs_mean = NA_real_, obs_stddev = NA_real_, 
                         obs_n = NA_real_, model_sum = NA_real_,
                         model_sqsum = NA_real_, model_n = NA_real_)
                    )
                
                return(out)
                })
             )
      
    sparsedist$model_mean <- sparsedist$model_sum / sparsedist$model_n
    row.names(sparsedist) <- 1:nrow(sparsedist)
    sparsedist <- sparsedist[,c("year","step","area","data_type","function_f",
                                "component","age","length","obs_mean",
                                "obs_stddev","obs_n","model_sum",
                                "model_mean","model_sqsum","model_n")]
    
  }else{
    sparsedist <- NULL
  }
  return(sparsedist)
}


#' Summarises the fit of a gadget3 model to survey indices
#'
#' @param reports Reported output of a G3 model
#' @return Data frame or NULL
#' @export
g3f_sidat <- function(reports){
  
  if (any(grepl('^.+_surveyindices_.+__num$|^.+_surveyindices_.+__wgt$', names(reports)))){
    
    si_str <- '(c|a)dist_([A-Za-z]+)_([A-Za-z]+)_(.+)_(model|obs)__(params|num|wgt)(.[0-9]+)?$'
    
    sidat_params <- 
      as.data.frame(
        do.call("rbind", 
                reports[grepl('^.+_surveyindices_.+__params$', names(reports))])
        ) 
    names(sidat_params) <- c("intercept", "slope")
    
    sidat_params <- 
      within(sidat_params, {
        index = gsub(si_str, '\\2', row.names(sidat_params));
        type = gsub(si_str, '\\3', row.names(sidat_params));
        fleet = gsub(si_str, '\\4', row.names(sidat_params))
      })
    
    sidat <- 
      do.call("rbind",
              lapply(
                reports[grep('^(a|c)dist_surveyindices_.+__(num|wgt)$',names(reports))],
                FUN = as.data.frame.table, stringsAsFactors = FALSE
                )
      ) 
    
    sidat <- 
      within(sidat, {
        index = gsub(si_str, '\\2', row.names(sidat));
        type = gsub(si_str, '\\3', row.names(sidat));
        fleet = gsub(si_str, '\\4', row.names(sidat));
        origin = gsub(si_str, '\\5', row.names(sidat));
        name = gsub(si_str, '\\2.\\4', row.names(sidat))})
    
    ## Pivot wider and merge in parameters
    sidat <- 
      stats::reshape(
        sidat, 
        timevar = "origin",
        idvar = c("length", "time", "area", "name", "fleet", "type", "index"),
        direction = "wide"
        ) |> 
      merge(sidat_params, by = c("index", "type", "fleet")) |> 
      within( {
        predicted = ifelse(type == "log", 
                           exp(intercept)*Freq.model^slope,
                           intercept + Freq.model*slope)
      }) |> 
      extract_year_step() |> 
      split_length() 
    
    names(sidat)[match(c("Freq.model", "Freq.obs"), names(sidat))] <- c("number", "observed")
    row.names(sidat) <- 1:nrow(sidat)
    
    return(sidat[,c("name", "year", "step", "area", "length", "lower", "upper", 
                    "fleet", "index", "type", "intercept", "slope",
                    "observed", "number", "predicted")])
    
    
  }else{
    sidat <- NULL
  }
  
}

#' Summarises the renewal and spawning reports for a gadget 3 model
#'
#' @param reports Reported output of a G3 model
#' @param rec_step The step(s) in which recruitment occurs. The renewal and spawning reports are not reset to zero in the steps following a recruitment event, so only the recruitment steps should be included. Either:
#'   \itemize{
#'     \item a vector of steps applied to all stocks, e.g. \code{rec_step = 2} or \code{rec_step = 1:4};
#'     \item a named list of steps for specific stocks (the names of the list), e.g. \code{rec_step = list(imm1 = 1, imm2 = 1:2)}. All steps are returned for stocks not included in the list;
#'     \item NULL, all steps are returned for all stocks.
#'   }
#' @return Data frame or NULL
#' @export
g3f_stock.recruitment <- function(reports,
                                  rec_step = 1){

  sr_str <- 'detail_(.+)__(spawnednum|renewalnum)$'

  if (!(is.null(rec_step) || is.numeric(rec_step) ||
        (is.list(rec_step) && !is.null(names(rec_step)) && all(names(rec_step) != "")))){
    stop("rec_step must be NULL, a vector of steps or a named list of steps for each stock")
  }

  if (any(grepl(sr_str, names(reports)))){

    report_rec_names <- names(reports)[grepl(sr_str, names(reports))]

    if (is.list(rec_step)){
      missing_stocks <- setdiff(names(rec_step), gsub(sr_str, "\\1", report_rec_names))
      if (length(missing_stocks) > 0){
        warning("No renewal or spawning reports found for the following stocks in rec_step: ",
                paste(missing_stocks, collapse = ", "))
      }
    }

    stock_rec <-
      lapply(
        stats::setNames(report_rec_names, report_rec_names),
        function(x){

          ## Note, spawn always enters the youngest age class, but renewal can
          ## enter any (i think) therefore filter out all age classes with
          ## 0 numbers summed across years... seems more robust than just
          ## filtering out zeros
          age_ind <- which(gadget3::g3_array_agg(reports[[x]],
                                                 agg = "sum",
                                                 margins = "age") > 0)

          if (length(age_ind) == 0) return(NULL)

          ## Recruitment steps for this stock (NULL for all steps)
          steps <- if (is.list(rec_step)) rec_step[[gsub(sr_str, "\\1", x)]] else rec_step

          out <- gadget3::g3_array_agg(reports[[x]],
                                       agg = "sum",
                                       margins = c("year", "step", "area", "age"),
                                       age = names(age_ind))

          out |>
            as.data.frame.table(stringsAsFactors = FALSE, responseName = "recruitment") |> 
            within({
              stock = gsub(sr_str, '\\1', x);
              age = as.numeric(gsub("age", "", age));
              year = as.numeric(year);
              step = as.numeric(step)
            }) -> out

          if (!is.null(steps)) out <- out[out$step %in% steps, ]

          return(out[,c("stock", "year", "step", "area", "age", "recruitment")])
          
        })
    
    rec_out <- do.call("rbind", stock_rec)
    if (!is.null(rec_out)) row.names(rec_out) <- NULL
    
  }else rec_out <- NULL

  return(rec_out)
}

#' Extracts the suitabilities from the reports from a gadget3 model
#'
#' @param reports Reported output of a G3 model
#' @param stock_names Vector of stock names for which the suits will be retrieved. If NULL, all stocks reported will be included.
#' @return Data frame 
#' @export
g3f_suitability <- function(reports, stock_names = NULL){
  
  if (is.null(stock_names))
    stock_names <- stock_pred_names(reports)$stocks 
  
  suit_re <- NA
  # New-style g3a_suitability_report() output
  if (any(grepl('^suit_(.+)_(.+)__report$', names(reports)))){
    
    suit_re <- paste0(
      "^suit",
      "_(\\Q", paste(stock_names, collapse = "\\E|\\Q"), "\\E)",
      "_(.+)",
      "__report" )
    
  # Old-style output
  } else if (any(grepl('detail_(.+)__suit_(.+)$', names(reports)))){
    
    suit_re <- paste0(
      "^detail",
      "_(\\Q", paste(stock_names, collapse = "\\E|\\Q"), "\\E)",
      "__suit",
      "_(.+)" )
    
  }
    
  if (!is.na(suit_re)){
    
    report_suit_names <- names(reports)[grep(suit_re, names(reports))]
    
    suitability <-
      do.call("rbind",
              lapply(
                stats::setNames(report_suit_names, report_suit_names),
                function(x){
                  out <- as.data.frame.table(reports[[x]], 
                                             stringsAsFactors = FALSE,
                                             responseName = "suit")
                  out$report_name <- x
                  return(out)
                }
              ))
    
    suitability <- 
      within(suitability, {
        stock = gsub(suit_re, "\\1", report_name);
        fleet = gsub(suit_re, "\\2", report_name)
      }) |> 
      extract_year_step() |> 
      split_length() |> 
      replace_inf(group_col = "stock") |> 
      identity()
    
    if (all(c("upper", "lower") %in% names(suitability))){
      suitability$length <- (suitability$upper + suitability$lower)/2
    }
    
    if ('age' %in% names(suitability)){
      suitability$age <- as.numeric(gsub('age', '', suitability$age))
    }
    
    row.names(suitability) <- 1:nrow(suitability)
    
    suit_cols <- c("year","step","area","stock","fleet",
                   "predator_length","length","age","suit")
    
    suitability <- suitability[, intersect(suit_cols, names(suitability))]
    
  } else {
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
  return(suitability)
}

#' Summarises the likelihood components of a gadget3 model by year and step
#'
#' @param reports Reported output of a G3 model
#' @param data_env The model data environment, e.g. \code{environment(model)} for an R model. Required for sparse data likelihood components.
#' @param params The model parameters (data frame with a 'value' column). Required for the weights of sparse data likelihood components.
#' @return Data frame or NULL. One row per likelihood component, year and step, containing the unweighted negative log-likelihood ('value') and the component weight.
#' @export
g3f_likelihood <- function(reports, data_env = NULL, params = NULL){
  
  lik_re <- "^nll_(adist|cdist)_([A-Za-z]+)_(.+)__"
  
  ## Likelihood
  if (any(grepl(lik_re, names(reports)))){
    
    lik_keys <- grep(paste0(lik_re, "(wgt|num)$"), names(reports), value = TRUE) 
    
    lik_vals <- 
      do.call("rbind", 
              lapply(stats::setNames(lik_keys, lik_keys), function(x){
                  as.data.frame.table(reports[[x]], 
                                      stringsAsFactors = FALSE, 
                                      responseName = "value") |> 
                  within( {
                    measurement = gsub(lik_re, "\\4", x);
                    id = gsub("(wgt|num)$", "", x)
                  } ) |> 
                  identity()
              }))
      
    lik_keys <- grep(paste0(lik_re, "weight$"), names(reports), value = TRUE)
    
    lik_vals <- 
      merge(x = lik_vals,
            y = do.call("rbind", 
                        lapply(stats::setNames(lik_keys, lik_keys), function(x){
                          as.data.frame.table(reports[[x]], 
                                              stringsAsFactors = FALSE, 
                                              responseName = "weight") |> 
                            within( {
                              id = gsub("weight$", "", x)
                            }) |> 
                            identity()
                        })),
            by = c("time", "id"))
    
    lik_vals <- 
      lik_vals |> 
      within( {
        component = gsub(lik_re, "\\3", id);
        data_type = gsub(lik_re, "\\1_\\2", id)
      } )
    
    likelihood <-
      lik_vals[,names(lik_vals) != "id"] |> 
      rbind(
        reports$nll_understocking__wgt |> 
          as.data.frame.table(stringsAsFactors = FALSE,
                              responseName = "value") |> 
          within( {
            component = "understocking";
            data_type = "model_preystocks";
            measurement = "wgt"
            weight = NA_real_  ## NEED TO ADD THIS TO REPORTING
          } )
      ) |> 
      extract_year_step()
      
    
    ## Sparse data
    if (any(grepl('^nll_(a|c)sparse_', names(reports)))){
      
      if (is.null(data_env) || is.null(params)){
        
        warning("The sparse distribution likelihood cannot be compiled without the data_env and params arguments")
        sp_data <- NULL
        
      } else {
      
        sp_re_obs <- 'nll_(asparse|csparse)_([A-Za-z]+)_(.+)__(year|step)$'
        
        ## Observations taken from model environment, model output taken from reports
        sparse_reports <- mget(ls(data_env)[grep(sp_re_obs, ls(data_env))], 
                               envir = data_env)
        
        ## Year and steps for each component
        sp_data <- 
          do.call("rbind",
                  lapply(stats::setNames(names(sparse_reports),
                                         names(sparse_reports)), function(x){
                                           out <- data.frame(value = sparse_reports[[x]])
                                           out$column <- gsub(sp_re_obs, "\\4", x)
                                           out$id <- gsub("__(.+)$", "", x)
                                           return(out)
                                         }))
        
        ## NLL for each component
        sp_data <- 
          do.call("rbind",
                  lapply(split(sp_data, sp_data$id), function(x){
                    x_re <- "nll_(asparse|csparse)_([A-Za-z]+)_(.+)$"
                    out <- 
                      utils::unstack(x, "value ~ column") |> 
                      within( {
                        component = gsub(x_re, '\\3', unique(x$id));
                        data_type = gsub(x_re, '\\1_\\2', unique(x$id));
                        value = 0;
                        weight = 0;
                        measurement = NA_character_
                      })
                    
                    ## Now add in the NLL
                    nll_re <- paste0(unique(x$id), "__nll")
                    par_re <- paste0(gsub("^nll_", "", unique(x$id)), "_weight")
                    
                    ## If linear regression, append to final timestep
                    if (grepl("_linreg_", nll_re)){
                      out$value[nrow(out)] <- reports[[nll_re]][["nll"]]
                      out$weight[nrow(out)] <- params$value[[par_re]]
                      out <- out[nrow(out),]
                    } else {
                      out$value <- reports[[nll_re]]
                      out$weight <- params$value[[par_re]]
                    }
                    
                    return(out)
                  })
          )
          
      }
      likelihood <- likelihood |> rbind(sp_data)
    }
    
    likelihood <- likelihood[order(likelihood$component),
                             c("year", "step", "component", 
                               "data_type", "measurement", "weight", 
                               "value")]
    
    return(likelihood)
    
  }else 
    return(NULL) 
}

#' Collates the abundance and mean weight reports for each stock in a gadget3 model
#'
#' @param reports Reported output of a G3 model
#' @return Data frame or NULL. One row per stock, year, step, area, age and length group, containing the abundance and mean weight.
#' @export
g3f_stock_reports <- function(reports){
  
  sr_re <- "^(dstart|detail)_(.+)__(num|wgt)$"
  
  if (!any(grepl(sr_re, names(reports)))) return(NULL)
  
  num_reports <- reports_to_df(reports, 
                               re = "^(dstart|detail)_(.+)__num$",
                               responseName = "abundance",
                               labels = c(stock = "\\2"))
  
  wgt_reports <- reports_to_df(reports, 
                               re = "^(dstart|detail)_(.+)__wgt$",
                               responseName = "weight",
                               labels = c(stock = "\\2"))
  
  stock_reports <- 
    merge(x = num_reports,
          y = wgt_reports,
          by = setdiff(intersect(names(num_reports), names(wgt_reports)), 
                       c("abundance", "weight")),
          all.x = TRUE) |> 
    tidy_stock_reports()
  
  cols <- c("year", "step", "area", "stock", "age", 
            "length", "lower", "upper", "avg.length", 
            "abundance", "weight")
  
  return(
    order_by(stock_reports[intersect(cols, names(stock_reports))],
             c("stock", "year", "step", "area", "age", "lower"))
  )
}

#' Collates the consumption of each stock by each predator (fleets and predator stocks) in a gadget3 model
#'
#' @param reports Reported output of a G3 model
#' @param stock_names Vector of (prey) stock names. If NULL, all stocks reported will be included.
#' @return Data frame or NULL. One row per prey stock, predator (column 'fleet'), year, step, area, age and length group (and predator length/age for predator stocks), containing the biomass and numbers consumed and the mean weight of the prey.
#' @export
g3f_fleet_reports <- function(reports, stock_names = NULL){
  
  if (is.null(stock_names))
    stock_names <- stock_pred_names(reports)$stocks
  
  stock_re <- paste0("_(\\Q", paste(stock_names, collapse = "\\E|\\Q"), "\\E)")
  
  if (any(grepl("__cons$", names(reports)))){
    # detail_(prey)_(pred)__cons arrays
    consumption_re <- paste0("^detail", stock_re, "_(.+)__cons$")
  } else {
    # Pre-predation stype __predby_ arrays
    consumption_re <- paste0("^detail", stock_re, "__predby_(.+)$")
  }
  
  if (!any(grepl(consumption_re, names(reports)))) return(NULL)
  
  cons_reports <- reports_to_df(reports, 
                                re = consumption_re,
                                responseName = "biomass_consumed",
                                labels = c(stock = "\\1", fleet = "\\2"))
  
  wgt_reports <- reports_to_df(reports, 
                               re = paste0("^(dstart|detail)", stock_re, "__wgt$"),
                               responseName = "weight",
                               labels = c(stock = "\\2"))
  
  fleet_reports <- 
    merge(x = cons_reports,
          y = wgt_reports,
          by = setdiff(intersect(names(cons_reports), names(wgt_reports)), 
                       c("biomass_consumed", "weight")),
          all.x = TRUE) |> 
    tidy_stock_reports()
  
  fleet_reports$number_consumed <- 
    ifelse(fleet_reports$biomass_consumed == 0, 
           0, 
           fleet_reports$biomass_consumed / fleet_reports$weight)
  
  if ("predator_age" %in% names(fleet_reports)){
    fleet_reports$predator_age <- as.numeric(gsub("age", "", fleet_reports$predator_age))
  }
  
  cols <- c("year", "step", "area", "stock", "fleet", "age", 
            "length", "lower", "upper", "avg.length", 
            "predator_length", "predator_age", 
            "biomass_consumed", "number_consumed", "weight")
  
  return(
    order_by(fleet_reports[intersect(cols, names(fleet_reports))],
             c("stock", "fleet", "year", "step", "area", "age", "lower"))
  )
}

## Common processing of stock and fleet reports: length groups, year/step and age
tidy_stock_reports <- function(data){
  
  data <- 
    data |> 
    split_length() |> 
    replace_inf(group_col = "stock") |> 
    extract_year_step()
  
  if (all(c("upper", "lower") %in% names(data))){
    data$avg.length <- (data$lower + data$upper)/2
  }
  
  if ("age" %in% names(data)){
    data$age <- as.numeric(gsub("age", "", data$age))
  }
  
  return(data)
}

#' Summarises the abundance and mean weight of each stock by length group
#'
#' @param reports Reported output of a G3 model
#' @param stock_reports Output of \code{\link{g3f_stock_reports}}
#' @return Data frame or NULL
#' @export
g3f_stock.full <- function(reports, 
                           stock_reports = g3f_stock_reports(reports)){
  
  if (is.null(stock_reports)) return(NULL)
  
  stock_reports$biomass <- stock_reports$abundance * stock_reports$weight
  
  stock_full <- 
    aggregate_by(stock_reports, 
                 by = c("year", "step", "area", "stock", "avg.length"), 
                 cols = c("abundance", "biomass"))
  
  stock_full$mean_weight <- stock_full$biomass / stock_full$abundance
  stock_full$mean_weight[is.na(stock_full$mean_weight)] <- 0
  
  names(stock_full)[match(c("avg.length", "abundance"), names(stock_full))] <- 
    c("length", "number")
  
  return(
    order_by(stock_full[intersect(c("year", "step", "area", "stock", 
                                     "length", "number", "mean_weight"),
                                   names(stock_full))],
             c("year", "step", "area", "stock", "length"))
  )
}

#' Summarises the abundance, length and weight of each stock by age
#'
#' @param reports Reported output of a G3 model
#' @param stock_reports Output of \code{\link{g3f_stock_reports}}
#' @return Data frame or NULL
#' @export
g3f_stock.std <- function(reports, 
                          stock_reports = g3f_stock_reports(reports)){
  
  if (is.null(stock_reports)) return(NULL)
  
  by <- intersect(c("year", "step", "area", "stock", "age"), names(stock_reports))
  grp <- do.call(paste, c(stock_reports[by], sep = "\r"))
  
  stock_reports <- 
    within(stock_reports, {
      mean_length <- ave(abundance * avg.length, grp, FUN = sum) / ave(abundance, grp, FUN = sum)
      sqdev <- abundance * (avg.length - mean_length)^2
      lengths <- abundance * avg.length
      biomass <- abundance * weight
    })
  
  stock_std <- 
    aggregate_by(stock_reports, 
                 by = by, 
                 cols = c("abundance", "lengths", "sqdev", "biomass")) |> 
    within( {
      mean_weight <- biomass / abundance
      stddev_length <- sqrt(sqdev / abundance)
      mean_length <- lengths / abundance
      number <- abundance
    })
  
  return(
    order_by(stock_std[c(by, "number", "mean_length", "stddev_length", "mean_weight")], by)
  )
}

#' Summarises the consumption of, and the resulting mortality on, each stock by age
#'
#' @param reports Reported output of a G3 model
#' @param fleet_reports Output of \code{\link{g3f_fleet_reports}}
#' @param stock_std Output of \code{\link{g3f_stock.std}}
#' @param fleet_names Names of the fleets (predators that are not stocks)
#' @return Data frame or NULL. Consumption and mortality are by fleets only. If the model includes predator stocks, the consumption and mortality by all predators are included in the '_total' columns.
#' @export
g3f_stock.prey <- function(reports, 
                           fleet_reports = g3f_fleet_reports(reports),
                           stock_std = g3f_stock.std(reports),
                           fleet_names = get_fleet_names(reports)){
  
  if (is.null(fleet_reports) || is.null(stock_std)) return(NULL)
  
  by <- intersect(c("year", "step", "area", "stock", "age"), names(fleet_reports))
  cols <- c("number_consumed", "biomass_consumed")
  
  ## Consumption by fleets
  stock_prey <- 
    aggregate_by(fleet_reports[fleet_reports$fleet %in% fleet_names, ], 
                 by = by, 
                 cols = cols) |> 
    merge(stock_std[c(by, "number")], by = by, all.x = TRUE)
  
  stock_prey$mortality <- 
    calc_mortality(stock_prey$number_consumed, stock_prey$number, 
                   get_step_size(reports, stock_prey$step))
  
  ## Add in biological consumption if included in model
  if (!all(fleet_reports$fleet %in% fleet_names)){
    
    total <- aggregate_by(fleet_reports, by = by, cols = cols)
    names(total)[match(cols, names(total))] <- paste0(cols, "_total")
    
    stock_prey <- merge(stock_prey, total, by = by, all.x = TRUE)
    stock_prey$mortality_total <- 
      calc_mortality(stock_prey$number_consumed_total, stock_prey$number, 
                     get_step_size(reports, stock_prey$step))
    
  }
  
  return(order_by(stock_prey, by))
}

#' Summarises the consumption of each stock (prey) by each predator (fleets and predator stocks) by length group
#'
#' @param reports Reported output of a G3 model
#' @param fleet_reports Output of \code{\link{g3f_fleet_reports}}
#' @param stock_reports Output of \code{\link{g3f_stock_reports}}
#' @param stock_full Output of \code{\link{g3f_stock.full}}
#' @param suitability Output of \code{\link{g3f_suitability}}
#' @return Data frame or NULL
#' @export
g3f_predator.prey <- function(reports,
                              fleet_reports = g3f_fleet_reports(reports),
                              stock_reports = g3f_stock_reports(reports),
                              stock_full = g3f_stock.full(reports, stock_reports),
                              suitability = g3f_suitability(reports)){
  
  if (is.null(fleet_reports) || is.null(stock_reports)) return(NULL)
  
  names(suitability)[names(suitability) == "length"] <- "avg.length"
  
  ## Merge in the suitabilities of each prey-predator combination, joining on 
  ## the columns used by that suitability (e.g. no predator_length for fleets)
  predator_prey <- 
    split(fleet_reports, 
          list(fleet_reports$stock, fleet_reports$fleet), 
          drop = TRUE) |> 
    lapply(function(x){
      suits <- suitability[suitability$stock %in% x$stock & 
                             suitability$fleet %in% x$fleet, , drop = FALSE]
      if (nrow(suits) == 0) return(transform(x, suit = NA_real_))
      suits <- suits[, vapply(suits, function(y) any(!is.na(y)), logical(1)), drop = FALSE]
      return(merge(x, suits, all.x = TRUE))
    }) |> 
    rbind_fill()
  
  ## Abundance of the prey
  by <- intersect(c("year", "step", "stock", "area", "age", "length"), names(stock_reports))
  predator_prey <- merge(predator_prey, stock_reports[c(by, "abundance")], by = by, all.x = TRUE)
  
  predator_prey$suit[predator_prey$number_consumed %in% 0] <- 0
  
  ## TO-CHECK: for predator stocks, fleet_reports has a row per predator length
  ## and age group, so each prey cell (abundance * weight) is counted once per 
  ## predator group when harv.bio is summed below. This may inflate harv.bio,
  ## suit (and fleet.info) for predator stocks. Fleets are unaffected.
  predator_prey$harv.bio <- predator_prey$abundance * predator_prey$weight * predator_prey$suit
  
  predator_prey <- 
    aggregate_by(predator_prey,
                 by = c("year", "step", "area", "stock", "fleet", "avg.length"),
                 cols = c("number_consumed", "biomass_consumed", "harv.bio"))
  names(predator_prey)[names(predator_prey) == "avg.length"] <- "length"
  
  predator_prey <- 
    merge(predator_prey, stock_full, 
          by = intersect(c("year", "step", "area", "stock", "length"), names(stock_full)), 
          all.x = TRUE) |> 
    within( {
      mortality <- calc_mortality(number_consumed, number, get_step_size(reports, step))
      tb <- number * mean_weight
      suit <- harv.bio / tb
      suit <- ifelse(is.finite(suit), suit, 0)
    })
  
  names(predator_prey)[match(c("fleet", "stock"), names(predator_prey))] <- c("predator", "prey")
  
  cols <- c("year", "step", "area", "prey", "predator", "length", 
            "number_consumed", "biomass_consumed", "harv.bio", 
            "mortality", "tb", "suit")
  
  return(
    order_by(predator_prey[intersect(cols, names(predator_prey))], 
             c("year", "step", "area", "prey", "predator", "length"))
  )
}

#' Summarises the catch, harvestable biomass and harvest rate of each predator (fleets and predator stocks)
#'
#' @param reports Reported output of a G3 model
#' @param predator_prey Output of \code{\link{g3f_predator.prey}}
#' @return Data frame or NULL
#' @export
g3f_fleet.info <- function(reports, 
                           predator_prey = g3f_predator.prey(reports)){
  
  if (is.null(predator_prey)) return(NULL)
  
  fleet_info <- 
    aggregate_by(predator_prey, 
                 by = c("year", "step", "area", "predator"),
                 cols = c("harv.bio", "biomass_consumed"))
  
  names(fleet_info)[match(c("predator", "biomass_consumed"), names(fleet_info))] <- 
    c("fleet", "amount")
  
  fleet_info$amount[is.na(fleet_info$amount)] <- 0
  fleet_info$harv.rate <- ifelse(fleet_info$amount == 0, 0, fleet_info$amount / fleet_info$harv.bio)
  
  return(order_by(fleet_info, c("year", "step", "area", "fleet")))
}

#' Summarises the suitability of each stock to all fleets combined, by length group
#'
#' @param reports Reported output of a G3 model
#' @param predator_prey Output of \code{\link{g3f_predator.prey}}
#' @param fleet_names Names of the fleets (predators that are not stocks)
#' @param stock_full Output of \code{\link{g3f_stock.full}}. Only used to return a table of NA suitabilities if there are no fleets.
#' @return Data frame. The suitability of each fleet is weighted by its catch.
#' @noRd
g3f_harv.suit <- function(reports,
                          predator_prey = g3f_predator.prey(reports),
                          fleet_names = get_fleet_names(reports),
                          stock_full = g3f_stock.full(reports)){
  
  cols <- c("year", "step", "area", "stock", "length")
  
  if (is.null(predator_prey) || !any(predator_prey$predator %in% fleet_names)){
    if (is.null(stock_full)) return(NULL)
    return(transform(unique(stock_full[cols]), suit = NA))
  }
  
  harv_suit <- predator_prey[which(predator_prey$predator %in% fleet_names & 
                                     predator_prey$biomass_consumed > 0), ]
  harv_suit$weighted_suit <- harv_suit$biomass_consumed * harv_suit$suit
  
  harv_suit <- 
    aggregate_by(harv_suit, 
                 by = c("year", "step", "area", "prey", "length"),
                 cols = c("weighted_suit", "biomass_consumed"))
  
  harv_suit$suit <- harv_suit$weighted_suit / harv_suit$biomass_consumed
  names(harv_suit)[names(harv_suit) == "prey"] <- "stock"
  
  return(order_by(harv_suit[c(cols, "suit")], cols))
}

#' Summarises the annual catch and fishing mortality of each stock
#'
#' @description The fishing mortality (F) of each stock is the annual F averaged over a range of ages (Fbar),
#' calculated separately for each stock and area. F in a step is \eqn{-log(1 - C/N)}, where \eqn{C} is the number
#' caught by fleets and \eqn{N} is the number at the beginning of the step. F at age is summed over the steps of the
#' year, and averaged (unweighted) over the age range. To calculate Fbar for groups of stocks, e.g. the immature
#' and mature stocks of a species, see \code{\link{g3_fit.update_fbar}}.
#'
#' @param reports Reported output of a G3 model
#' @param stock_prey Output of \code{\link{g3f_stock.prey}}
#' @param stock_full Output of \code{\link{g3f_stock.full}}. Only used to return a table of NA mortalities if there are no fleets.
#' @param f_by_age The age range over which F is averaged. Either:
#'   \itemize{
#'     \item NULL (default), the oldest age of each stock, i.e. assuming the top of a logistic selection curve;
#'     \item a vector of the minimum and maximum age (or all ages, e.g. \code{4:7}) applied to all stocks;
#'     \item a named list with an age range for specific stocks, e.g. \code{list(imm = c(1, 2), mat = c(5, 6))}. 
#'       The oldest age is used for stocks not included in the list.
#'   }
#'   F is NA for stocks with no ages within the range.
#' @return Data frame
#' @export
g3f_f.by.year <- function(reports, 
                          stock_prey = g3f_stock.prey(reports),
                          stock_full = g3f_stock.full(reports),
                          f_by_age = NULL){
  
  by <- c("stock", "year", "area")
  
  if (!(is.null(f_by_age) || is.numeric(f_by_age) || 
        (is.list(f_by_age) && !is.null(names(f_by_age)) && all(names(f_by_age) != "")))){
    stop("f_by_age must be NULL, a vector of ages or a named list of ages for each stock")
  }
  
  if (is.null(stock_prey) || nrow(stock_prey) == 0){
    if (is.null(stock_full)) return(NULL)
    return(transform(unique(stock_full[by]), F = NA))
  }
  
  ## Age range of each stock, the oldest age by default
  stocks <- unique(stock_prey$stock)
  if (is.list(f_by_age) && length(setdiff(names(f_by_age), stocks)) > 0){
    warning("The following stocks in f_by_age have no fleet catches: ",
            paste(setdiff(names(f_by_age), stocks), collapse = ", "))
  }
  age_range <- 
    lapply(stats::setNames(stocks, stocks), function(s){
      ages <- if (is.list(f_by_age)) f_by_age[[s]] else f_by_age
      if (is.null(ages)) ages <- max(stock_prey$age[stock_prey$stock == s])
      return(range(ages))
    })
  
  in_range <- 
    stock_prey$age >= vapply(age_range, min, numeric(1))[stock_prey$stock] & 
    stock_prey$age <= vapply(age_range, max, numeric(1))[stock_prey$stock]
  
  no_ages <- setdiff(stocks, stock_prey$stock[in_range])
  if (length(no_ages) > 0){
    warning("F is NA for the following stocks, which have no ages within f_by_age: ",
            paste(no_ages, collapse = ", "))
  }
  
  f_by_year <-
    aggregate_by(stock_prey, by = by, cols = c("biomass_consumed", "number_consumed")) |>
    merge(calc_fbar(stock_prey[in_range, ], by = c("stock", "area")), by = by, all.x = TRUE)
  
  names(f_by_year)[match(c("biomass_consumed", "number_consumed"), names(f_by_year))] <- 
    c("catch", "num.catch")
  
  return(order_by(f_by_year, by))
}

#' Summarises the annual abundance, biomass, catch, fishing mortality and recruitment of each stock
#'
#' @param reports Reported output of a G3 model
#' @param steps Which steps to include in the annual output? Vector of int.
#' @param stock_full Output of \code{\link{g3f_stock.full}}
#' @param harv_suit Suitability of each stock to all fleets combined (weighted by catch), by year, step, area and length group. Computed internally by default.
#' @param f_by_year Output of \code{\link{g3f_f.by.year}}
#' @param stock_recruitment Output of \code{\link{g3f_stock.recruitment}}
#' @return Data frame or NULL
#' @export
g3f_res.by.year <- function(reports,
                            steps = 1,
                            stock_full = g3f_stock.full(reports),
                            harv_suit = g3f_harv.suit(reports, stock_full = stock_full),
                            f_by_year = g3f_f.by.year(reports, stock_full = stock_full),
                            stock_recruitment = g3f_stock.recruitment(reports)){
  
  if (is.null(stock_full)) return(NULL)
  
  by <- c("stock", "year", "area")
  
  res_by_year <- 
    stock_full[stock_full$step %in% steps, ] |> 
    merge(harv_suit, by = c("stock", "year", "step", "area", "length"), all.x = TRUE) |> 
    within( {
      harv.biomass <- number * suit * mean_weight
      harv.biomass[is.na(harv.biomass)] <- 0
      total.biomass <- number * mean_weight
      total.number <- number
    }) |> 
    aggregate_by(by = c(by, "step"), 
                 cols = c("total.number", "total.biomass", "harv.biomass")) |> 
    merge(f_by_year, by = by, all.x = TRUE)
  
  if (is.null(stock_recruitment)){
    res_by_year$recruitment <- rep(NA_real_, nrow(res_by_year))
  } else {
    res_by_year <- 
      merge(res_by_year, 
            aggregate_by(stock_recruitment, by = by, cols = "recruitment"),
            by = by, all.x = TRUE)
  }
  
  return(order_by(res_by_year, c(by, "step")))
}
