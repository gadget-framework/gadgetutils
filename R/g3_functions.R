#' @export
g3_retro <- function(dir, model, params, num.years = 5){
  
  ## Create directories for each retro
  if (!dir.exists(file.path(dir, 'RETRO'))) dir.create(file.path(dir, 'RETRO'))
  
  ## Get final year from model environment
  model_actions <- attr(model, 'actions')
  time_action_pos <- locate_g3_actions(model_actions, 'time')
  time_action <- model_actions[[time_action_pos$pos]]
  
  ## Set parameters for g3a time
  start_year <- get('start_year', environment(time_action[['000']]))
  end_year <- get('end_year', environment(time_action[['000']]))
  steps <- get('step_lengths', environment(time_action[['000']]))
  
  retro_models <- vector('list', length = num.years)
  for (i in 1:num.years){
    
    ## Replace existing time action using peel year
    new_actions <- model_actions
    new_actions[[time_action_pos$pos]] <- g3a_time(start_year = start_year,
                                                   end_year = end_year - i,               
                                                   steps = steps)
    
    
    ## Model and parameters
    new_model <- g3_to_tmb(new_actions)
    new_params <- attr(new_model, 'parameter_template')
    new_params$value <- params$value[names(new_params$value)]
    new_params$optimise <- params[match(new_params$switch, params$switch), 'optimise']
    
    retro_models[[i]] <- list()
    retro_models[[i]]$model <- new_model
    retro_models[[i]]$param <- new_params
    
    save_obj(tmb_model = new_model,
             file = file.path(dir, 'RETRO', paste0('model_', i, '.Rdata')))
    save_obj(param_init = new_params,
             file = file.path(dir, 'RETRO', paste0('param_init_', i, '.Rdata')))
    
  }
  
  ## Inner function to run the model - using weights from year 0 model
  run_func <- function(param, tmb_model){
    
    ## Compile objective function
    obj_fun <- g3_tmb_adfun(tmb_model, param)
    
    ## Optimise model
    fit.opt <- optim(g3_tmb_par(param),
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
