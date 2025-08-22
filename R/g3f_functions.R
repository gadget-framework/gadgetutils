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
        within(data, {
          stock = gsub("-999", NA_character_, stock);
          stock_re = gsub("-999", NA_character_, stock_re);
          predator_length = gsub("-999", NA_character_, predator_length);
          residuals = ifelse(observed == 0, NA, observed - predicted)
        }) |> 
        extract_year_step() |> 
        split_length(mean_length_col = "avg.length", replace_inf = TRUE) 
      
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
#' @param stock_rec_step Named list where each element contains the step(s) for a specific stock (name of the element). E.g. list(bli_imm = 1, bli_dummy = 2). If NULL or a stock is not included/misspelled, all steps are returned.  
#' @return Data frame or NULL
#' @export
g3f_stock.recruitment <- function(reports,
                                  stock_rec_step = NULL){
  
  sr_str <- 'detail_(.+)__(spawnednum|renewalnum)$'
  
  if (any(grepl(sr_str, names(reports)))){
    
    stock_rec <- 
      lapply(
        stats::setNames(names(reports)[grepl(sr_str, names(reports))],
                        names(reports)[grepl(sr_str, names(reports))]),
        function(x){
          
          ## Note, spawn always enters the youngest age class, but renewal can 
          ## enter any (i think) therefore filter out all age classes with 
          ## 0 numbers summed across years... seems more robust than just
          ## filtering out zeros
          age_ind <- which(g3_array_agg(reports[[x]],
                                        agg = "sum",
                                        margins = "age") > 0)
          
          if (length(age_ind) == 0) return(NULL)
          
          ## Were renewal/spawning steps specified?
          ## If yes
          if (gsub(sr_str, "\\1", x) %in% names(stock_rec_step)){
            out <- gadget3::g3_array_agg(reports[[x]], 
                                         agg = "sum", 
                                         margins = c("year", "step", "area", "age"), 
                                         step = stock_rec_step[[gsub(sr_str, "\\1", x)]],
                                         age = names(age_ind))  
          }else{
            ## If no
            out <- gadget3::g3_array_agg(reports[[x]], 
                                         agg = "sum", 
                                         margins = c("year", "step", "area", "age"),
                                         age = names(age_ind))
          }
          
          out |> 
            as.data.frame.table(stringsAsFactors = FALSE, responseName = "recruitment") |> 
            within({
              stock = gsub(sr_str, '\\1', x);
              age = as.numeric(gsub("age", "", age));
              year = as.numeric(year);
              step = as.numeric(step)
            }) -> out
          
          return(out[,c("stock", "year", "step", "area", "age", "recruitment")])
          
        })
    
    rec_out <- do.call("rbind", stock_rec)
    row.names(rec_out) <- 1:nrow(rec_out)
    
  }else rec_out <- NULL

  return(rec_out)
}

g3f_suitability <- function(reports){
  
  suit_str <- 
  
  ## Suitability
  if (any(grepl('^suit_(.+)_(.+)__report', names(tmp)))){
    suit_re <- paste0(
      "^suit",
      "_(\\Q", paste(stock_names, collapse = "\\E|\\Q"), "\\E)",
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
  
}