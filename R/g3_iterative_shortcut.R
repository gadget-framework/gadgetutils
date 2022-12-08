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
    bind_rows(
      lapply(res[grepl('(^adist|^cdist)_surveyindices_(.+)_obs__(wgt$|num$)', names(res))],
             function(x){
               estimate_weights_adist(as.data.frame.table(x, stringsAsFactors = FALSE))
             }) %>% 
        dplyr::bind_rows(.id = 'comp')
      ) %>% 
    mutate(weight = 1/ss,
           comp = gsub('obs__(num$|wgt$)', 'weight', .$comp))
  
  return(out)
}

#' @export
estimate_weights_cdist <- function(dat){
  tmp <- 
    dat %>% 
    filter(Freq > 0) %>% 
    group_by(time) %>% 
    mutate(p = Freq/sum(Freq)) 
  
  if ()
  tmp %>% 
    left_join(tmp %>% 
                group_by(age,length) %>% 
                summarise(phat=mean(p), .groups = 'drop'), 
              by = c('age', 'length')) %>% 
    ungroup() %>% 
    summarise(ss = sum((p-phat)^2/n()))
}

#' @export
estimate_weights_adist <- function(dat){
  
  dat %>% 
    extract_year_step() %>% 
    dplyr::mutate(time = year) %>% ## TODO Mutate time to include 'step'
    modelr::add_predictions(model = loess(log(Freq) ~ time, 
                                          data=.,                                   
                                          span = 0.25)) %>% 
    summarise(ss=sum((log(Freq)-pred)^2)/n())
}


