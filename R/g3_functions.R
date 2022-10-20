#' Perform multiple optimisation runs of a model, removing a year for each
#'
#' @param dir Directory to store output from RETRO runs
#' @param model A G3 model, produced by g3_to_tmb()
#' @param params Initial parameters to use with the model
#' @param num.years How many years back to go
#' @return List of results for each optimisation run
#' @export
g3_retro <- function(dir, model, params, num.years = 5){
  
  ## Create directories for each retro
  if (!dir.exists(file.path(dir, 'RETRO'))) dir.create(file.path(dir, 'RETRO'))
  
  retro_models <- vector('list', length = num.years)
  for (i in 1:num.years){
    ## Model and parameters
    new_params <- attr(model, 'parameter_template')
    new_params$value <- params$value[names(new_params$value)]
    new_params$optimise <- params[match(new_params$switch, params$switch), 'optimise']

    if (!('retro_years' %in% rownames(new_params))) {
        stop("No 'retro_years' parameter, cannot apply g3_retro. Make sure g3a_time is included in your model, with retro_years as to the default value")
    }
    new_params['retro_years', 'value'] <- i
    
    retro_models[[i]] <- list()
    retro_models[[i]]$model <- model
    retro_models[[i]]$param <- new_params
    
    save_obj(tmb_model = model,
             file = file.path(dir, 'RETRO', paste0('model_', i, '.Rdata')))
    save_obj(param_init = new_params,
             file = file.path(dir, 'RETRO', paste0('param_init_', i, '.Rdata')))
    
  }
  
  ## Inner function to run the model - using weights from year 0 model
  run_func <- function(param, tmb_model){
    
    ## Compile objective function
    # NB: We need to recreate the obj_fun to use new unoptimised parameters, i.e. retro_years.
    obj_fun <- g3_tmb_adfun(tmb_model, param)
    
    ## Optimise model
    fit.opt <- optim(obj_fun$par,
                     obj_fun$fn,
                     obj_fun$gr,
                     method = 'BFGS',
                     control = list(trace = 2,
                                    maxit = 1000,
                                    reltol = .Machine$double.eps^2))
    
    ## Output parameters
    p <- g3_tmb_relist(param, fit.opt$par)
    param$value[names(p)] <- p
    return(param)

  }
  
  out <- parallel::mclapply(retro_models, function(x) run_func(x$param, x$model), mc.cores = parallel::detectCores())
  
  ## Save final params
  for (i in 1:length(out)){
    save_obj(param_final = out[[i]], file = file.path(dir, 'RETRO', paste0('param_final_', i, '.Rdata')))
  }
  return(out)
}


save_obj <- function(..., file) {
  x <- list(...)
  save(list = names(x), file = file, envir = list2env(x))
}
