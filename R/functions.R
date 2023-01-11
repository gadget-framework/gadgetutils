## -----------------------------------------------------------------------------
##
## Miscellaneous functions used in gadgetutils
##
## -----------------------------------------------------------------------------

#' @export
echo_message <- function(...){
  system(sprintf("echo '%s'", paste0(..., collapse = '')))
}

save_obj <- function(..., file) {
  x <- list(...)
  save(list = names(x), file = file, envir = list2env(x))
}

## A number of g3_* functions in gadgetutils return a list of optimised parameter
## data.frames. Sometimes NULL will be returned from mclapply if a process 
## does not complete. In these cases it would be useful to issue a warning and substitute
## in the input parameters to continue the optimisation
#' @export
check_null_params <- function(params_out, param_in){
  
  ## Get index of NULL elements
  params_null <- sapply(params_out, FUN = is.null)
  if (all(params_null)) stop("None of the optimised parameters were returned")
  
  ## Loop over list elements to see if NULL was returned
  for (i in seq_along(params_out)){
    
    if (!is.null(params_out[[i]])){
      attr(params_out[[i]], 'summary') <- cbind(attr(params_out[[i]], 'summary'),
                                                data.frame(return_complete = 1))
    }
    else{
      
      ## NULL was returned
      warning(paste("A NULL was returned for component", names(params_out[i])))
      
      ## First sort out summary attribute
      sumat <- attributes(params_out[[i]])
      if (is.null(sumat)){
        sumat <- attr(params_out[[which(!params_null)[1]]], 'summary')
        sumat[,c('optim_complete', 'fn_calls', 'gd_calls', 'convergence', 'score')] <- NA
      }
      sumat$return_complete <- 0
      
      ## Second change parameters
      params_out[[i]] <- params_in[[i]]
      params_out[[i]]$value[params_out[[i]]$optimise == 1] <- NA
      
      ## Update attribute
      attr(params_out[[i]], 'summary') <- sumat
      
    }
    
  }
  return(params_out)
}