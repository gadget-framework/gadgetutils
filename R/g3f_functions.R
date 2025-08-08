#' Summarises the fit of a gadget 3 model to survey indices
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
      split_length() |> 
      replace_inf()
    
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