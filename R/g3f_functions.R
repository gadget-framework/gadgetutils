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
#' @param stock_rec_step Named list where each element contains the step(s) for a specific stock (name of the element). E.g. list(bli_imm = 1, bli_dummy = 2). If NULL or a stock is not included/misspelled, all steps are returned.  
#' @return Data frame or NULL
#' @export
g3f_stock.recruitment <- function(reports,
                                  stock_rec_step = NULL){
  
  sr_str <- 'detail_(.+)__(spawnednum|renewalnum)$'
  
  if (any(grepl(sr_str, names(reports)))){
    
    report_rec_names <- names(reports)[grepl(sr_str, names(reports))]
    
    stock_rec <- 
      lapply(
        stats::setNames(report_rec_names, report_rec_names),
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

