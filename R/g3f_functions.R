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
