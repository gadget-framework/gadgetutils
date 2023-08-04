#' Jitter analysis
#' 
#' \code{g3_jitter} performs a jitter analysis.
#' @name g3_jitter
#' @param gd Directory to store output
#' @param outdir Directory name within gd to store run outputs
#' @param model A G3 model, produced by g3_to_tmb() or g3_to_r()
#' @param params Initial parameters to use with the model, this should be a TMB parameter template i.e. attr(tmb_model, 'parameter_template')
#' @param njits Number of jitters to run, defaults to 10
#' @param jitter_fraction The fraction of jittering for a value
#' @param pattern_to_ignore Regular expression of parameters to avoid jittering
#' @param within_bounds Logical, if TRUE, jittered values that fall outside parameter bounds (params$lower, params$upper) will be adjusted to fall within the bounds
#' @param use_parscale Logical indicating whether optim(control$parscale) should be used
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param serial_compile g3_tmb_adfun will be run in serial mode (i.e., not in parallel), potentially helping with memory issues
#' @param mc.cores The number of cores to use, defaults to the number available
#' @return A list of optimised parameter data frames (one for each jitter)
#' @export
g3_jitter <- function(gd, outdir = 'JITTER',
                      model, params, 
                      njits = 10, jitter_fraction = 0.1, 
                      pattern_to_ignore = NULL, within_bounds = TRUE,
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
  ## Jitter the params
  ## ---------------------------------------------------------------------------
  
  if (mc.cores > 1){
    jitpar_in <- parallel::mclapply(stats::setNames(1:njits, 1:njits), 
                                    function(x){
                                      return(
                                        jitter_params(params, 
                                                      jitter_fraction, 
                                                      pattern_to_ignore, 
                                                      within_bounds)
                                      )
                                    }, mc.cores = mc.cores) 
  }
  else{
    jitpar_in <- lapply(stats::setNames(1:njits, 1:njits), 
                        function(x){
                          return(
                            jitter_params(params, 
                                          jitter_fraction, 
                                          pattern_to_ignore, 
                                          within_bounds)
                          )
                        })
  }
  
  
  ## Save parameters
  save(jitpar_in, file = file.path(out_path, 'jitpar_in.Rdata'))
  
  ## ---------------------------------------------------------------------------
  ## Run optimisation
  ## ---------------------------------------------------------------------------
  
  echo_message('##  RUNNING JITTER ANALYSIS\n')
  
  jitpar_out <- run_g3_optim(model, jitpar_in,
                             use_parscale, method, control,
                             serial_compile, mc.cores)
  
  ## Add class
  class(jitpar_out) <- c('g3.jitter', class(jitpar_out))
  
  ## Save output
  save(jitpar_out, file = file.path(out_path, 'jitpar_out.Rdata'))
  
  ## Summary of optimisation settings and run details
  check_params_out(jitpar_out, 'jitter_id') %>% 
    write.g3.file(out_path, 'optim.summary.jitter')
  
  return(jitpar_out)
  
}

#' Creates a jittered g3_tmb parameter template
#'
#' @param params A g3 tmb parameter template
#' @param jitter_fraction The fraction of jittering for a value
#' @param pattern_to_ignore Regular expression of parameters to avoid jittering
#' @param within_bounds Logical, if TRUE, jittered values that fall outside parameter bounds (params$lower, params$upper) will be adjusted to fall within the bounds
#' @return Jittered set of parameters
#' @export
jitter_params <- function(params, 
                          jitter_fraction = 0.1, 
                          pattern_to_ignore = NULL,
                          within_bounds = TRUE){
  
  ## Store unjittered parameter values
  params$old_value <- params$value
  
  ## Extract parameter values and bounds
  pars <- gadget3::g3_tmb_par(params)
  pars_down <- gadget3::g3_tmb_lower(params)
  pars_up <- gadget3::g3_tmb_upper(params)
  
  ## Are any parameters exponentiated?
  pars_exp <- which(grepl('_exp$', names(pars)))
  if (length(pars_exp) > 0) pars[pars_exp] <- exp(pars[pars_exp])
  
  ## Calculate a shift term for all parameters to be jittered
  jst <- stats::runif(n = length(pars), 
                      min = 0 - jitter_fraction, 
                      max = 0 + jitter_fraction)
  ## Jitter
  out <- pars + jst * pars
  
  ## Log any exponentiated values
  if (length(pars_exp) > 0) out[pars_exp] <- log(out[pars_exp])
  
  ## Check for boundary issues
  if (within_bounds){
    upi <- out >= pars_up
    dni <- out <= pars_down
    
    ## Adjust to 10% within boundary
    if (any(upi)) out[upi] <- pars_up[upi] - 0.1 *  (pars_up[upi] - pars[upi])
    if (any(dni)) out[dni] <- pars_down[dni] + 0.1 * (pars[dni] - pars_down[dni])
  }
  
  ## Insert back to parameter template
  p <- gadget3::g3_tmb_relist(params, out)
  params$value[names(p)] <- p
  
  ## Revert any values coming from 'pattern_to_ignore'
  if (!is.null(pattern_to_ignore)){
    ignore_pos <- grepl(pattern_to_ignore, params$switch)
    if (any(ignore_pos)){ params$value[ignore_pos] <- params$old_value[ignore_pos] }
    else warning("A 'pattern_to_ignore' was supplied however the pattern was not found in params$switch")
  }
  
  return(params[, names(params) != 'old_value'])
  
}
