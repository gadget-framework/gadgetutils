#' Sensitivity analysis
#' 
#' \code{g3_sensitivity} performs a sensitivity analysis for a single parameter.
#' @name g3_sensitivity
#' @param gd Directory to store output
#' @param outdir Directory name within gd to store run outputs
#' @param model A G3 model, produced by g3_to_tmb() or g3_to_r()
#' @param params Initial parameters to use with the model, this should be a TMB parameter template i.e. attr(tmb_model, 'parameter_template')
#' @param par_name A glob-like string to match parameter names, see the 'name_spec' argument for \code{\link[gadget3]{g3_init_val}}
#' @param par_values A vector of values to substitute into params
#' @param use_parscale Logical indicating whether optim(control$parscale) should be used
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param serial_compile g3_tmb_adfun will be run in serial mode (i.e., not in parallel), potentially helping with memory issues
#' @param return_input_parameters If TRUE, returns a list of modified parameter dataframes (i.e. without optimisations)
#' @param run_optimisations If TRUE, will return a list of fits using the modified (non-optimised) parameters
#' @param mc.cores The number of cores to use, defaults to the number available
#' @return A list of parameter dataframes (either optimised or not) or a list of fits based on the unoptimised parameters
#' @export
g3_sensitivity <- function(gd, outdir = 'SENS',
                           model, params,
                           par_name, par_values,
                           use_parscale = TRUE,
                           method = 'BFGS',
                           control = list(),
                           serial_compile = FALSE,
                           return_input_parameters = FALSE,
                           run_optimisations = TRUE,
                           mc.cores = parallel::detectCores()){
  
  ## Some checks:
  ## We want the TMB parameter template
  if (!inherits(params, 'data.frame')){
    stop("Error: expected a TMB parameter template for 'params'")
  }
  ## Compile the TMB model if the R version is passed in...
  if (inherits(model, 'g3_r') && run_optimisations){
    if (is.null(attr(model, 'actions'))) stop("Error: please provide the c++ model or the R model with 'actions' attribute")
    model <- gadget3::g3_to_tmb(actions = attr(model, 'actions'))
  }
  
  ## Create output directory if it does not exist
  out_path <- file.path(gd, outdir)
  if (!exists(out_path)) dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  
  ## Copying g3_init_val for contructing the regex
  # Parse par_name --> regex
  name_re <- paste0(vapply(strsplit(par_name, ".", fixed = TRUE)[[1]], function (part) {
    # [1979-1984] - range match
    m <- regmatches(part, regexec('^\\[(\\d+)[:-](\\d+)\\]$', part))
    if (all(vapply(m, length, numeric(1)) == 3)) {
      m <- m[[1]]
      return(paste0(
        '(?:',
        paste(seq(as.numeric(m[[2]]), as.numeric(m[[3]])), collapse = "|"),
        ')'))
    }
    
    # # - numeric match
    part <- gsub("#", "\\E\\d+\\Q", part, fixed = TRUE)
    
    # *  - string match
    part <- gsub("*", "\\E.*\\Q", part, fixed = TRUE)
    
    # | - or part
    part <- gsub("|", "\\E|\\Q", part, fixed = TRUE)
    
    # Make sure by default text in part is quoted, scope | above
    return(paste0('(?:\\Q', part, '\\E)'))
  }, character(1)), collapse = "\\.")
  
  name_re <- paste0('^', name_re,'(_exp)?', '$')
  m <- regmatches(params$switch, regexec(name_re, params$switch))
  
  matches <- sapply(m, length) > 0
  if (!any(matches)) {
    warning("g3_init_val('", par_name, "') didn't match any parameters")
    return(params)
  }
  
  # Make boolean vector for all places to auto_exp 
  auto_exp <- vapply(m, function(x) length(x) >= 2 && x[[2]] == '_exp', logical(1))
  
  if (any(params[matches, 'optimise']) && !return_input_parameters && run_optimisations){
    warning("The parameter(s) '", par_name, "' will be optimised")
    if (any(par_values < params[matches, 'lower']) || any(par_values > params[matches, 'upper']))
      stop("The selected parameter will be optimised, but the input 'par_values' lie outside the parameters lower and upper bounds")
  }
  
  ## Create list of input parameter dataframes
  senpar_in <- list()
  for (i in 1:length(par_values)){
    
    senpar_in[[paste(gsub('\\*\\.|\\.#', '', par_name), par_values[i], sep = '_')]] <- params
    senpar_in[[i]][matches, 'value'] <- par_values[i]
    if (any(auto_exp)) senpar_in[[i]][auto_exp, 'value'] <- sapply(senpar_in[[i]][auto_exp, 'value'], log)
    
  }
  
  ## Return manipulated parameters if not running optimisations
  if (return_input_parameters){
    echo_message('\n----------------------------- RETURNING INPUT PARAMETERS\n')
    return(senpar_in)
  }
  
  if (!run_optimisations){
    echo_message('\n-------------------------------------- COLLATING THE FIT\n')
    out <- parallel::mclapply(senpar_in, function(x, model){
      return( g3_fit(model, x) )
    }, model = model, mc.cores = mc.cores)
    return(out)
  }
  
  ## Save parameters
  save(senpar_in, file = file.path(out_path, 'senpar_in.Rdata'))
  
  ## ---------------------------------------------------------------------------
  ## Run optimisation
  ## ---------------------------------------------------------------------------
  
  echo_message('\n----------------------------- RUNNING SENSITIVITY ANALYSIS\n')
  
  senpar_out <- run_g3_optim(model, senpar_in,
                             use_parscale, method, control,
                             serial_compile, mc.cores)
  
  ## Add class
  class(senpar_out) <- c('g3.sens', class(senpar_out))
  
  ## Save output
  save(senpar_out, file = file.path(out_path, 'senpar_out.Rdata'))
  
  ## Summary of optimisation settings and run details
  check_params_out(senpar_out, 'sensitivity_id') %>% 
    write.g3.file(out_path, 'optim.summary.sensitivity')
  
  return(senpar_out)
  
}