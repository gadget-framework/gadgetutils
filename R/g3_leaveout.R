#' Leave-one-out analysis. 
#' 
#' Perform a leave-one-out analysis. This involves running a set of optimisations, 
#' each one of which has a component (or set of components) turned off (i.e. weight = 0).
#'
#' \code{g3_leaveout} performs a leave-one-out analysis.
#' @name g3_leaveout
#' @param gd Directory to store output
#' @param outdir Directory name within gd to store run outputs
#' @param model A G3 model, produced by g3_to_tmb() or g3_to_r()
#' @param params Initial parameters to use with the model, this should be a TMB parameter template i.e. attr(tmb_model, 'parameter_template')
#' @param grouping List of component names to group together. This allows you to leave out multiple components at a time, e.g. length and age-length distributions from a particular fleet
#' @param use_parscale Logical indicating whether optim(control$parscale) should be used
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param serial_compile g3_tmb_adfun will be run in serial mode (i.e., not in parallel), potentially helping with memory issues
#' @param mc.cores The number of cores to use, defaults to the number available
#' @return A list of optimised parameter data frames (one for each group)
#' @export
g3_leaveout <- function(gd, outdir = 'LOCV', 
                        model, params,
                        grouping = list(),
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
  
  ## Compile the TMB model if the R version is passed in...
  if (inherits(model, 'g3_r')){
    model <- gadget3::g3_to_tmb(actions = attr(model, 'actions'))
  }
  
  ## Create output directory if it does not exist
  out_path <- file.path(gd, outdir)
  if (!exists(out_path)) dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  
  ## ---------------------------------------------------------------------------
  ## Sort out groupings
  ## ---------------------------------------------------------------------------
  
  ## Get all individual components from params object
  comps <- params[grepl('_weight$', params$switch), 'switch']
  comps <- gsub('^.dist_[a-z]+_|_weight', '', comps)
  
  ## Group components
  group_comps <- do.call('c', grouping)
  
  ## Check whether the supplied groupings all exist
  if (!all(group_comps %in% comps)){
    stop(paste0("The following components were provided in 'grouping' but not found in 'params': ",
                paste(group_comps[!group_comps %in% comps], collapse = ', ')))
  }
  
  ## Add ungrouped components to the list of groupings
  all_comps <- c(grouping,
                 split(comps[!(comps %in% group_comps)],
                       comps[!(comps %in% group_comps)]))
  
  ## Set up input parameters for each grouping that will be dropped
  leaveout_params_in <- lapply(all_comps, function(x, params){
    out <- params
    out$value[grepl(paste(paste0('_', x, '_weight'), collapse = '|'), params$switch)] <- 0
    return(out)
  }, params = params)
  
  ## Save output
  save(leaveout_params_in, file = file.path(out_path, 'leaveout_params_in.Rdata'))
  
  ## ---------------------------------------------------------------------------
  ## Run the optimisations
  ## ---------------------------------------------------------------------------
  
  echo_message('##  RUNNING LEAVE-ONE-OUT ANALYSIS\n')
  
  leaveout_params_out <- run_g3_optim(model, leaveout_params_in,
                                      use_parscale, method, control,
                                      serial_compile, mc.cores)
  
  ## Add class
  class(leaveout_params_out) <- c('g3.leaveout', class(leaveout_params_out))
  
  ## Save output
  save(leaveout_params_out, file = file.path(out_path, 'leaveout_params_out.Rdata'))
  
  ## Summary of optimisation settings and run details
  summary <- lapply(names(all_comps), function(x){
    return(
      cbind(data.frame(components_left_out = paste(all_comps[[x]], collapse = ', '),
                       group = x),
            attr(leaveout_params_out[[x]], 'summary'), stringsAsFactors = FALSE)
    )
  })
  
  do.call('rbind', summary) %>% write.g3.file(out_path, 'optim.summary.leaveout')
  
  return(leaveout_params_out)
  
}
