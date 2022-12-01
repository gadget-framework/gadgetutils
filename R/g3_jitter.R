#' Jitters a parameters values
#'
#' @param param A row from a g3 tmb parameter template
#' @param jitter_fraction The fraction of jittering for a value
#' @return Jittered set of parameters
#' @export
jitter_param <- function(param, jitter_fraction = 0.05){
  
  ## A function to jitter initial parameters
  if (is.infinite(param$lower) || is.infinite(param$upper)){
    warning(paste("Parameter", names(param), "not jittered as bounds are infinite."))
    return(param)
  }
  if (jitter_fraction <= 0){
    stop("The jitter_fraction parameter should be greater than zero.")
  }
  
  ## Is the parameter exponentiated
  exp_param <- grepl('_exp$', param$switch)
  
  old_val <- param$value[[1]]
  if (exp_param) old_val <- exp(old_val)
  
  ## Jitter shift term using a uniform distribution centered
  jst <- runif(n = 1, min = 0 - jitter_fraction, max = 0 + jitter_fraction)
  new_val <- old_val + jst * old_val
  
  ## Transform back if exponentiated
  if (exp_param) new_val <- log(new_val)
  
  ## Within bounds?
  if (new_val > param$upper){
    warning(paste0('The jittered param ', param$switch, 'was above the upper bound'))
  } 
  else{
    if (new_val < param$lower){
      warning(paste0('The jittered param ', param$switch, 'was below the lower bound'))
    } 
  }

  ## Convert to value
  param$value[[1]] <- new_val
  return(param)
  
}

#' Jitters a parameters values
#'
#' @param param A g3 tmb parameter template
#' @param jitter_fraction The fraction of jittering for a value
#' @param patterns_to_ignore Regular expression of parameters to avoid jittering
#' @return Jittered set of parameters
#' @export
g3_jitter <- function(params, jitter_fraction = 0.1, patterns_to_ignore = '_weight'){
  
  if (!all(c('value', 'lower', 'upper') %in% names(params))){
    print("The columns 'value', 'lower', 'upper' were not all found in the 'params' parameter.")
    stop(invisible())
  }
  
  tmp <- lapply(split(params, params$switch), function(x, jitter_fraction){
  #  print(x)
    if (any(grepl(patterns_to_ignore, names(x))) || !x$optimise){ return(x) }
    else{ return(jitter_param(x, jitter_fraction = jitter_fraction)) }

  }, jitter_fraction = jitter_fraction)
  
  out <- do.call('rbind', tmp)
  return(out[match(params$switch, out$switch),])
}

