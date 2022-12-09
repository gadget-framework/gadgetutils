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
    dplyr::mutate(weight = 1/ss,
                  comp = gsub('obs__(num$|wgt$)', 'weight', .$comp))
  
  return(out)
}

#' @export
estimate_weights_cdist <- function(dat){
  tmp <- 
    dat %>% 
    dplyr::filter(Freq > 0) %>% 
    dplyr::group_by(time) %>% 
    dplyr::mutate(p = Freq/sum(Freq)) 
  
  # Adding this hack to avoid crashes when age column is not present in data. Note that there is no grouping by area, nor over stocks or time (but stocks and time are intentional, I think). Note also that the order of character vector under will define the column order. Remove this comment when you have read it. -Mikko
  group_cols <- intersect(c("age", "length"), colnames(tmp))
  
  tmp %>% 
    dplyr::left_join(
      tmp %>% 
        dplyr::group_by(dplyr::across(tidyselect::all_of(group_cols))) %>% 
        dplyr::summarise(phat=mean(p), .groups = 'drop'),
      by = group_cols
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::summarise(ss = sum((p-phat)^2/dplyr::n()), .groups = 'drop')
}

#' @export
estimate_weights_adist <- function(dat){
  
  dat %>% 
    extract_year_step() %>% 
    dplyr::mutate(time = year) %>% ## TODO Mutate time to include 'step'
    modelr::add_predictions(
      model = loess(log(Freq) ~ time, data=., span = 0.25)) %>% 
    dplyr::summarise(ss=sum((log(Freq)-pred)^2)/dplyr::n())
}


