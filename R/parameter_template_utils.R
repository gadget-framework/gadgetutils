#' Set initial guesses for parameter values in g3_cpp parameter template
#'
#' @param params Parameter template from g3_to_tmb()
#' @param pattern Regular expression of parameter names in template to change
#' @param value Initial value for all parameters matching pattern
#' @param lower Lower bound to set for all parameters matching pattern
#' @param upper Upper bound to set for all parameters matching pattern
#' @param optimise 1 if parameters matching pattern should be optimised, zero otherwise.
#' @return Updated parameter template
#' @export
g3_init_guess <- function(params, pattern, 
                          value = 0, lower = -999, upper = 999, 
                          optimise = 0){
  
  
  if (!inherits(params, 'data.frame')){
    stop("The 'params' argument should be a data.frame")
  }
 
  ## Check whether parameter is exponentiated
  is_param_exp <- any(grepl(paste0(pattern, '_exp'), params$switch))
  if (is_param_exp) pattern <- paste0(pattern, '_exp')
 
  ## Are the parameters time- or age-varying
  is_param_varying <- any(grepl(paste0(pattern, '\\.[0-9]'), params$switch))
  
  ## Create new pattern for parameters varying with eg. time
  v_pattern <- ifelse(is_param_varying, 
                      paste0(pattern, '\\.[0-9]'),
                      pattern)
  
  ## Make sure parameter values are within bounds
  if (value <= lower && optimise == 1){ 
    warning(paste("The 'value' provided for", pattern,
                  "is <= the 'lower' bound and is therefore adjusted to fall within the bounds"))
    value <- max(lower, value) + 0.01*(upper - lower)
  }
  else{
    if (value >= upper && optimise == 1){ 
      warning(paste("The 'value' provided for", pattern,
                    "is >= the 'upper' bound and is therefore adjusted to fall within the bounds"))
      value <- min(upper, value) - 0.01*(upper - lower)
    }  
  }
  
  ## Need to log value if parameter is exponentiated
  if (is_param_exp){ 
    value <- log(value)
    lower <- log(lower)
    upper <- log(upper)
  }
  
  ## Fill in the horizontal template
  params[grepl(v_pattern, params$switch), 'value'] <- value
  params[grepl(v_pattern, params$switch), 'lower'] <- lower
  params[grepl(v_pattern, params$switch), 'upper'] <- upper
  params[grepl(v_pattern, params$switch), 'parscale'] <- diff(c(lower, upper))
  params[grepl(v_pattern, params$switch), 'optimise'] <- as.logical(optimise)
  
  return(params)
}

#' @export
g3_add_parscale <- function(parameters){
  par <- g3_tmb_parscale(parameters)
  out <- utils::relist(par, unclass(parameters$value[parameters$optimise]))
  parameters$parscale[match(names(out), parameters$switch)] <- out
  return(parameters)
}

#' @export
transform_bounded <- function(params, ...){
  .Defunct()
}

#' Convert value into it's normalised form for use with bounded()
#'
#' @param x Raw value
#' @param lower Lower bound for x
#' @param upper Upper bound for x
#' @return Normalised value for use with bounded()
#' @export
value_from_bounds <- function(x, lower, upper){
  if (x == lower && x == upper) return(x)
  else  return(log((upper - lower)/(x - lower) - 1))
}

#' Convert normalised bounded form back to value
#'
#' @param x Input normalised value
#' @param lower Lower bound for x
#' @param upper Upper bound for x
#' @return Raw value for x without bounding
#' @export
eval_bounded <- function(x, lower, upper){
  return(lower + (upper - lower)/(1 + exp(x)))
}





