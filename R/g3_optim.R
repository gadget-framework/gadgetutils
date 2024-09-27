#' Wrapper for stats::optim that returns the g3_cpp parameter template with optimised values
#'
#' @param model Either a g3 model of class 'g3_cpp', i.e., model <- gadget3::g3_to_tmb(actions), or an objective function from gadget3::g3_tmb_adfun(model, param)
#' @param params The parameter template for a g3_cpp classed model. 
#' @param use_parscale Should parscale be used g3_tmb_parscale(params), see \code{\link[stats]{optim}}
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param print_status Print step comments, does not suppress comments from optim call (see control$trace)
#' @param print_id A character string appended to the print statements (useful for when running g3_optim in parallel)
#' @param check_finite boolean Should we check that initial/final conditions produce finite results before running?
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
                     print_id = 'g3_optim',
                     check_finite = TRUE,
                     ...){
  
  ## Checks
  if (!inherits(params, 'data.frame')){
    cat("\nPlease provide a parameter data.frame, i.e. attr(tmb_model, 'parameter_template'\n")
    stop()
  }
  if (!inherits(model, 'g3_cpp') & !all(c('par','fn','gr') %in% names(model))){
    cat("\nPlease provide a 'g3_cpp' classed model, i.e. resulting from g3_to_tmb(actions), or an objective function resulting from g3_tmb_adfun()\n")
    stop()
  }
  if (!is.logical(use_parscale) || !is.logical(print_status)){
    cat("\nThe 'use_parscale' and 'print_status' arguments should be TRUE or FALSE\n")
    stop()
  }
  
  ## Prefix for printing
  if (nchar(print_id > 0)) print_id <- paste(' for', print_id)
  
  ## Create the objective function
  if (inherits(model, 'g3_cpp')){
    if (print_status) echo_message('##  Creating objective function', print_id)
    obj_fun <- gadget3::g3_tmb_adfun(model, params, ...)  
  }else{
    if (print_status) echo_message('##  Using the provided objective function', print_id)
    obj_fun <- model
  }
  
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
    x <- gadget3::g3_tmb_parscale(params)
    if (any(is.na(x))){
      # Put FALSE/TRUE back into parameter list (to get the pre-cpp names), extract names of TRUE values
      par_names <- names(Filter(
          function(y) identical(y, TRUE),
          g3_tmb_relist(params, is.na(x)) ))
      stop('Error NAs detected in parscale', print_id, ": ", paste(par_names, collapse = ", "))
    }
  }

  if (check_finite) {
    if (print_status) echo_message('##  Checking initial conditions', print_id)

    x <- obj_fun$fn(obj_fun$par)
    if (!is.finite(x)) {
      warning("Initial conditions", print_id, " produce a non-finite likelihood. Use g3a_trace_nan() to help locate when this happened: ", x)
    } else {
      x <- obj_fun$gr(obj_fun$par)
      if (any(!is.finite(x))) {
        # Flatten (params x 1) arrays, if that's what's returned
        if (length(dim(x)) > 1) x <- colSums(x, na.rm = FALSE)
        names(x) <- names(obj_fun$par)
        # Put FALSE/TRUE back into parameter list (to get the pre-cpp names), extract names of TRUE values
        par_names <- names(Filter(
            function(y) identical(y, TRUE),
            g3_tmb_relist(params, !is.finite(x)) ))
        warning("Initial conditions", print_id, " produce a non-finite gradient in: ", paste(par_names, collapse = ", "))
      }
    }
  }
  
  ## Run optimiser
  if (print_status) echo_message('##  Running optimisation', print_id)
  fit_opt <- try({
    
    optim(par = obj_fun$par, fn = obj_fun$fn, gr = obj_fun$gr,
          method = method,
          lower = parlow, upper = parhigh,
          control = control)
    
  }, silent = FALSE)
  
  ## Check whether the optimisation crashed
  if (inherits(fit_opt, 'try-error')){
    model.opt.fail <- TRUE
    warning("Optimisation", print_id, " failed: ", paste(as.character(fit_opt), collapse = " "))
    
    ## Construct failed fit_opt object
    fit_opt <- list(par = gadget3::g3_tmb_par(params),
                    counts = c(NA, NA),
                    convergence = NA,
                    value = NA)
    fit_opt$par[] <- NA
  }
  else{
    model.opt.fail <- FALSE
    if (print_status) echo_message('##  Optimisation finished', 
                                   print_id, 
                                   '. Convergence = ',
                                   ifelse(fit_opt$convergence == 0, TRUE, FALSE))
    
    if (check_finite) {
      if (print_status) echo_message('##  Checking optimised parameters', print_id)

      x <- obj_fun$fn(fit_opt$par)
      if (!is.finite(x)) {
        warning("Optimised parameters", print_id, " produce a non-finite likelihood. Use g3a_trace_nan() to help locate when this happened: ", x)
      } else {
        x <- obj_fun$gr(fit_opt$par)
        if (any(!is.finite(x))) {
          # Flatten (params x 1) arrays, if that's what's returned
          if (length(dim(x)) > 1) x <- colSums(x, na.rm = FALSE)
          names(x) <- names(fit_opt$par)
          # Put FALSE/TRUE back into parameter list (to get the pre-cpp names), extract names of TRUE values
          par_names <- names(Filter(
              function(y) identical(y, TRUE),
              g3_tmb_relist(params, !is.finite(x)) ))
          warning("Optimised parameters", print_id, " produce a non-finite gradient in: ", paste(par_names, collapse = ", "))
        }
      }
    }
  }
    
  ## Optimised parameters  
  p <- gadget3::g3_tmb_relist(params, fit_opt$par)
  params$value[names(p)] <- p
    
  ## Add summary of input/output to an attribute of params
  attributes(params)$summary <- 
    data.frame(type = 'g3.optim',
               id = print_id,
               method = method,
               maxiter = control$maxit,
               reltol = format(control$reltol, scientific = TRUE),
               optim_complete = ifelse(model.opt.fail, 0, 1),
               fn_calls = fit_opt$counts[1],
               gd_calls = fit_opt$counts[2],
               convergence = ifelse(model.opt.fail, NA, 
                                    ifelse(fit_opt$convergence == 0, TRUE, FALSE)),
               score = fit_opt$value,
               return_complete = '',
               comment = ''
               )
  
  return(params)
    
}
 

#' Wrapper for g3_optim. This is an internal function used by the g3_* functions
#'
#' @param model Either a g3 model of class 'g3_cpp', i.e., model <- gadget3::g3_to_tmb(actions), or an objective function from gadget3::g3_tmb_adfun(model, param)
#' @param params A list of parameter templates for a g3_cpp classed model.
#' @param use_parscale Should parscale be used g3_tmb_parscale(params), see \code{\link[stats]{optim}}
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param serial_compile g3_tmb_adfun will be run in serial mode (i.e., not in parallel), potentially helping with memory issues
#' @param mc.cores number of cores used, defaults to the number of available cores
#' @param ... Further arguments to be passed gadget3::g3_tmb_adfun, see \code{\link[gadget3]{g3_tmb_adfun}}
#' @return A g3_cpp parameter template with optimised values
#'
run_g3_optim <- function(model, params,
                         use_parscale, method, control,
                         serial_compile, mc.cores, ...){
  
  ## Compiling the model prior to g3_optim?
  if (serial_compile){
    echo_message('##  COMPILING MODEL AND CREATING ADFUN IN SERIAL\n')
    
    objfns <- lapply(stats::setNames(names(params), names(params)), function(x){
      echo_message(' - Compiling component ', x)
      return(gadget3::g3_tmb_adfun(model, params[[x]], ...))
    }) 
    
  }else{
    objfns <- list(model)
  }
  
  ## Now run g3_optim
  ## In serial?
  if (mc.cores == 1){
    
    out <- lapply(stats::setNames(names(params), names(params)), function(x){
      
      if (length(objfns) > 1) md <- x
      else md <- 1
      
      return(
        g3_optim(model = objfns[[md]],
                 params = params[[x]],
                 use_parscale = use_parscale,
                 method = method,
                 control = control,
                 print_status = TRUE,
                 print_id = x)  
      )
    })
  }else{
    gc()  # NB: Trigger garbage collect before forking, so we don't inadventently copy unused memory
    out <- parallel::mclapply(stats::setNames(names(params), names(params)), function(x){
      
      if (length(objfns) > 1) md <- x
      else md <- 1
      
      # parallel::mclapply eats warnings, so we have to capture them and replay them outside
      ww <- list()
      withCallingHandlers({
        out <- g3_optim(model = objfns[[md]],
                 params = params[[x]],
                 use_parscale = use_parscale,
                 method = method,
                 control = control,
                 print_status = TRUE,
                 print_id = x)
      }, warning = function (w) { ww <<- c(ww, list(w)) })

      return(list(out = out, ww = ww))
    }, mc.cores = mc.cores)
    out <- lapply(out, function(o) {
      if (inherits(o, "try-error")) {
        warning("Optimisation failed: ", as.character(o))
        return(o)
      }
      if (identical(names(o), c("out", "ww"))) {
        # Re-trigger warnings outside of mclapply
        for (w in o$ww) warning(w)
        return(o$out)
      }
      return(o)
    })
  }
  return(out)
}

#' Internal function to check whether g3_optim has succeeded or crashed. Used by the gadgetutils g3_* functions.
#' 
#' \code{check_params_out} checks for NULLs and try-errors in a list of g3 parameters
#' @name check_params_out
#' @param param_list List of g3 parameters
#' @param id_name replacement name for 'id' in the output
#' @return A data.frame that has collated the 'summary' attribute of each element of param_list.
#'
check_params_out <- function(param_list, id_name = NULL){
  
  stopifnot(is.list(param_list))
  
  classtype <- class(param_list)[grepl('^g3.', class(param_list))]
  if (length(classtype) == 0) classtype <- NA
  
  parsum <- lapply(names(param_list), function(x){
    
    ## NULL
    if (is.null(param_list[[x]])){
      
      warning(paste0('Optimisation failed: a NULL was produced for component: ', x))
      
      ## Construct the summary attribute if it does not exist
      if (!is.null(attr(param_list[[x]], 'summary'))){
        tmp <- attr(param_list[[x]], 'summary')
      }else{
        tmp <- data.frame(type = classtype,
                          id = x,
                          method = NA, maxiter = NA, reltol = NA,
                          optim_complete = NA, fn_calls = NA, gd_calls = NA,
                          convergence = NA, score = NA, 
                          return_complete = 0,
                          comment = 'Optimisation crashed: NULL issue')
      }
      
      tmp$type <- classtype
      tmp$id <- x
      tmp$return_complete = 0
      
    }else{
      ## try-error
      if (inherits(param_list[[x]], "try-error")){
        
        warning(paste0('Optimisation failed: a try-error class was produced for component: ', x))
        
        tmp <- data.frame(type = classtype,
                          id = x,
                          method = NA, maxiter = NA, reltol = NA,
                          optim_complete = NA, fn_calls = NA, gd_calls = NA,
                          convergence = NA, score = NA,
                          return_complete = 0,
                          comment = paste0('Optimisation crashed: try-error message: ', 
                                           attr(param_list[[x]], 'condition')))
      }else{
        
        tmp <- attr(param_list[[x]], 'summary')
        tmp$type <- classtype
        tmp$id <- x
        tmp$return_complete <- 1
        
      }
    }
    return(tmp)
  })
  
  out <- dplyr::bind_rows(parsum) %>% `rownames<-`(NULL)
  if (!is.null(id_name)) names(out)[names(out) == 'id'] <- id_name
  return(out)
  
}

