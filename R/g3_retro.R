#' Analytical retrospective
#'
#' \code{g3_retro} runs an analytical retrospective model fitting run. 
#' @name g3_retro
#' @param gd Directory to store output
#' @param outdir Directory name within gd to store run outputs
#' @param model A G3 model, produced by g3_to_r() or g3_tmb_adfun()
#' @param params Initial parameters to use with the model, this should be a TMB parameter template i.e. attr(tmb_model, 'parameter_template')
#' @param nyears How many years back to go
#' @param use_parscale Logical indicating whether optim(control$parscale) should be used
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param serial_compile g3_tmb_adfun will be run in serial mode (i.e., not in parallel), potentially helping with memory issues
#' @param mc.cores The number of cores to use, defaults to the number available
#' @return List optimised TMB parameter data.frames, the name of each list element corresponds to the number of years peeled in the retrospective analysis.
#' @export
g3_retro <- function(gd, outdir = 'RETRO',
                     model, params,
                     nyears = 5,
                     use_parscale = TRUE,
                     method = 'BFGS',
                     control = list(),
                     serial_compile = FALSE,
                     mc.cores = parallel::detectCores()){
  
  ## Some checks:
  ## We want the TMB parameter template
  if (!inherits(params, 'data.frame')){
    stop("Error: expected a TMB parameter template for 'params'")
  }
  
  ## Need the 'retro_years' parameter
  if (!('retro_years' %in% rownames(params))) {
    stop("No 'retro_years' parameter, cannot apply g3_retro. Make sure g3a_time is included in your model, with the 'retro_years' parameter set to its default value")
  }
  
  ## Compile the TMB model if the R version is passed in...
  if (inherits(model, 'g3_r')){
    model <- gadget3::g3_to_tmb(actions = attr(model, 'actions'))
  }
  
  ## Create output directory if it does not exist
  out_path <- file.path(gd, outdir)
  if (!exists(out_path)) dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  
  ## ---------------------------------------------------------------------------
  ## Create a parameter data.frame for each year peeled
  ## ---------------------------------------------------------------------------
  
  retro_params_in <- vector('list', length = nyears)
  for (i in 1:nyears){
    retro_params_in[[i]] <- params
    retro_params_in[[i]]$value$retro_years <- i
  }

  retro_params_in <- stats::setNames(retro_params_in, 1:nyears)
  
  ## Save input parameters
  save(retro_params_in, file = file.path(out_path, 'retro_params_in.Rdata'))
  
  ## ---------------------------------------------------------------------------
  ## Optimise the models
  ## ---------------------------------------------------------------------------
  
  echo_message('##  RUNNING ANALYTICAL RETROSPECTIVE\n')
  
  retro_params_out <- run_g3_optim(model, retro_params_in,
                                   use_parscale, method, control,
                                   serial_compile, mc.cores)
  
  ## Add class
  class(retro_params_out) <- c('g3.retro', class(retro_params_out))
  
  ## Save and write parameters
  save(retro_params_out, file = file.path(out_path, 'retro_params_out.Rdata'))
  
  ## Write summary of optimisation settings and run details
  summary <- lapply(names(retro_params_out), function(x){
    return(
      cbind(data.frame(nyears_peeled = x),
            attr(retro_params_out[[x]], 'summary'), 
            stringsAsFactors = FALSE)
    )
  })
  
  do.call('rbind', summary) %>% write.g3.file(out_path, 'optim.summary.retro')
  
  for (i in names(retro_params_out)){
    write.g3.param(retro_params_out[[i]],
                   out_path,
                   paste0('retro.params.out.', i, 'year')
                   )
  }
  
  return(retro_params_out)
}
