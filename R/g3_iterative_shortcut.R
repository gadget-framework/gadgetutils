#' @title Shortcut weights 
#' @description Estimates the weights for gadget3 likelihood components using shortcut (ie not iterative re-weighting) methods
#' @param model A gadget3 model of class 'g3_cpp' or 'g3_r'
#' @param params A gadget3 parameter dataframe
#' @details This function attempts to estimate the weights of gadget3 likelihood components. The shortcut weights are calculated by taking the inverse of the residual variance. How the residual variance is calculated depends on the type of likelihood distribution, i.e., whether it is a catch or abundance distribution. The residual variance calculations are performed by the helper functions 'estimate_weights_cdist' (catch distributions) and 'estimate_weights_adist' (abundance distributions) 
#' @export
estimate_weights <- function(model, params){
  
  if (inherits(model, 'g3_r')){
    model <- gadget3::g3_to_tmb(attr(model, 'actions'))
  }
  if (!inherits(params, 'data.frame')){
    stop('Expecting params to be a data.frame')
  }
  
  ## Collate reports from the model
  adfun <- gadget3::g3_tmb_adfun(model, params, type = 'Fun')
  res <- adfun$report(adfun$par)
  
  ## Catch distribution weights
  out <- 
    lapply(res[grepl('^cdist_(sumofsquares|multinomial)_(.+)_obs__(wgt$|num$)', names(res))],
           function(x){
             estimate_weights_cdist(as.data.frame.table(x, stringsAsFactors = FALSE))
           }) %>% 
    dplyr::bind_rows(.id = 'comp') %>% 
    dplyr::bind_rows(
      lapply(res[grepl('(^adist|^cdist)_surveyindices_(.+)_obs__(wgt$|num$)', names(res))],
             function(x){
               estimate_weights_adist(as.data.frame.table(x, stringsAsFactors = FALSE))
             }) %>% 
        dplyr::bind_rows(.id = 'comp')
    ) %>% 
    dplyr::mutate(weight = 1/.data$ss,
                  comp = gsub('obs__(num$|wgt$)', 'weight', .data$comp))
  
  return(out)
}

#' @title Residual variance for catch distributions
#' @description Calculates the residual variance for catch distributions using a GLM
#' @param dat Catch distribution model report (as a dataframe)
#' @export
estimate_weights_cdist <- function(dat){
  tmp <- 
    dat %>% 
    dplyr::filter(.data$Freq > 0) %>% 
    dplyr::group_by(.data$time) %>% 
    dplyr::mutate(p = .data$Freq/sum(.data$Freq)) 
  
  # Adding this hack to avoid crashes when age column is not present in data. 
  # Note that there is no grouping by area, nor over stocks or time (but stocks and time are intentional, I think). 
  # Note also that the order of character vector under will define the column order. Remove this comment when you have read it. -Mikko
  group_cols <- intersect(c("age", "length"), colnames(tmp))
  
  tmp %>% 
    dplyr::left_join(
      tmp %>% 
        dplyr::group_by(dplyr::across(tidyselect::all_of(group_cols))) %>% 
        dplyr::summarise(phat=mean(.data$p), .groups = 'drop'),
      by = group_cols
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::summarise(ss = sum((.data$p-.data$phat)^2/dplyr::n()), .groups = 'drop')
}

#' @title Residual variance for abundance distributions
#' @description This function fits a loess model to abundance data and returns the residual variance of the model fit
#' @param dat Abundance distribution model report (as a dataframe)
#' @export
estimate_weights_adist <- function(dat){
  tmp <-   
    dat %>% 
    extract_year_step() %>% 
    dplyr::mutate(time = .data$year)  
  tmp %>% ## TODO Mutate time to include 'step'
    modelr::add_predictions(
      model = stats::loess(log(Freq) ~ time, span = 0.75, data = tmp)) %>% 
    dplyr::summarise(ss=sum((log(.data$Freq)-.data$pred)^2)/dplyr::n())
}


