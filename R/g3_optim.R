#' Wrapper for stats::optim that returns the g3_cpp parameter template with optimised values
#'
#' @param model A g3 model of class 'g3_cpp', i.e., model <- g3_to_tmb(actions)
#' @param params The parameter template for a g3_cpp classed model
#' @param use_parscale Should parscale be used g3_tmb_parscale(params), see \code{\link[stats]{optim}}
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param print_status Print step comments, does not suppress comments from optim call (see control$trace)
#' @param print_id A character string appended to the print statements (useful for when running g3_optim in parallel)
#' @param ... Further arguments to be passed gadget3::g3_tmb_adfun, see \code{\link[gadget3]{g3_tmb_adfun}}
#' @return A g3_cpp parameter template with optimised values
#'
#' @export
g3_optim <- function(model, 
                     params, 
                     use_parscale = TRUE,
                     method = 'BFGS',
                     control = list(),
                     print_status = FALSE,
                     print_id = '',
                     ...){
  
  ## Checks
  if (!inherits(params, 'data.frame')){
    cat("\nPlease provide a parameter data.frame, i.e. attr(tmb_model, 'parameter_template'\n")
    stop()
  }
  if (!inherits(model, 'g3_cpp')){
    cat("\nPlease provide a 'g3_cpp' classed model, i.e. resulting from g3_to_tmb(actions)\n")
    stop()
  }
  if (!is.logical(use_parscale) || !is.logical(print_status)){
    cat("\nThe 'use_parscale' and 'print_status' arguments should be TRUE or FALSE\n")
    stop()
  }
  
  ## Prefix for printing
  if (print_status && nchar(print_id > 0)) print_id <- paste(' for', print_id)
  
  ## Create the objective function
  if (print_status) echo_message('##  Creating objective function', print_id)
  obj_fun <- gadget3::g3_tmb_adfun(model, params, ...)
  
  ## Configure bounds for optimisation
  if (method == 'L-BFGS-B'){
    parhigh <- gadget3::g3_tmb_upper(params)
    parlow <- gadget3::g3_tmb_lower(params)
  }
  else{
    parhigh <- Inf
    parlow <- -Inf
  }
 
  ## Add some defaults to the control list if they were not provided
  if (!('maxit' %in% names(control))) control$maxit <- 1000
  if (!('trace' %in% names(control))) control$trace <- 2
  if (!('reltol' %in% names(control))) control$reltol <- .Machine$double.eps^2
  
  ## Using parscale?
  if (use_parscale){
    control <- c(control, list(parscale = gadget3::g3_tmb_parscale(params)))
  }
  
  ## Run optimiser
  if (print_status) echo_message('##  Running optimisation', print_id)
  fit_opt <- optim(par = obj_fun$par, fn = obj_fun$fn, gr = obj_fun$gr,
                   method = method,
                   lower = parlow, upper = parhigh,
                   control = control)
  
  if (print_status) echo_message('##  Optimisation finished', 
                                 print_id, 
                                 '. Convergence = ',
                                 ifelse(fit_opt$convergence == 0, TRUE, FALSE))
  
  ## Optimised parameters  
  p <- gadget3::g3_tmb_relist(params, fit_opt$par)
  params$value[names(p)] <- p
  
  ## Add summary of input/output to an attribute of params
  attributes(params)$summary <- data.frame(method = method,
                                           max_iterations = control$maxit,
                                           reltol = control$reltol,
                                           function_calls = fit_opt$counts[1],
                                           gradient_calls = fit_opt$counts[2],
                                           convergence = ifelse(fit_opt$convergence == 0, TRUE, FALSE),
                                           score = fit_opt$value)
  
  return(params)
  
}
  