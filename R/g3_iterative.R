#' Perform multiple optimisation runs of a model, reweighting with each run
#'
#' @param gd Directory to store output
#' @param wgts Directory name within gd to store run outputs
#' @param r_model A G3 model, produced by g3_to_r()
#' @param params.in Initial parameters to use with the model
#' @param grouping List of component names to optmise together
#' @return Final set of parameters
#' @importFrom rlang .data
#' @importFrom stats optim
#' @export
g3_iterative <- function(gd,
                         wgts = 'WGTS',
                         r_model, 
                         tmb_model, 
                         params.in, 
                         grouping = list(),
                         opt_method = 'BFGS',
                         use_parscale = TRUE){
  
  out_path <- file.path(gd, wgts)
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
  
  ## -------------- Iterative re-weighting setup -------------------------------
  
  ## Run the R model to get the initial results
  init_weights <- g3_lik_out(r_model, params.in)
  
  if (is.na(attr(init_weights, 'nll'))){
    stop('The gadget model did not run')
  }
  
  write.g3.file(init_weights, out_path, 'lik.init')
  
  ## Setup the initial parameter files
  init_params <- g3_iterative_setup(init_weights, grouping = grouping)
  
  ## Write groupings 
  write.g3.file(init_params$grouping, out_path, 'lik.groupings')
  
  ## Initial parameters (weights included in parameter template)
  for (i in names(init_params$params)){
    write.g3.param(init_params$params[[i]], 
                   out_path, 
                   paste0('params.init.stage1.', i))
  }
  
  ## -------------- Run first stage of iterative re-weighting  -----------------
  
  stage1_params <- parallel::mclapply(init_params$params, 
                                      function(x){g3_iterative_run(x, 
                                                                   tmb_model, 
                                                                   opt_method,
                                                                   use_parscale)}, 
                                      mc.cores = parallel::detectCores())
  
  ## Summary of optimisation settings and run details
  summary1 <- 
    lapply(stage1_params, function(x) attr(x, 'summary')) %>% 
    dplyr::bind_rows(.id = 'group')
  
  write.g3.file(summary1, out_path, 'optim.summary.stage1')
  
  for (i in names(stage1_params)){
    attr(stage1_params[[i]], 'summary') <- NULL
    write.g3.param(stage1_params[[i]],
                   out_path,
                   paste0('params.out.stage1.', i))
  }
  
  
  
  ## ------------ Update weights for second round of re-weighting --------------
  
  ## Update weights for second round of re-weighting
  int_params <- parallel::mclapply(stage1_params, 
                                   function(x){ g3_lik_out(r_model, x) }, 
                                   mc.cores = parallel::detectCores()) %>% 
    g3_iterative_final()
  
  ## Write initial params for stage 2
  ## Note: is it necessary to write all these params out?
  ## Split between weights and params?
  
  for (i in names(int_params)){
    write.g3.param(int_params[[i]],
                   out_path,
                   paste0('params.init.stage2.', i))
  }
  
  ## ----------- Second round of re-weighting ----------------------------------
  
  stage2_params <- parallel::mclapply(int_params, 
                                      function(x){g3_iterative_run(x, 
                                                                   tmb_model,
                                                                   opt_method,
                                                                   use_parscale)}, 
                                      mc.cores = parallel::detectCores())
  
  #save(stage2_params, file = file.path(out.dir, 'stage2_params.Rdata'))
  
  ## Summary of optimisation settings and run details
  summary2 <- 
    lapply(stage2_params, function(x) attr(x, 'summary')) %>% 
    dplyr::bind_rows(.id = 'group')
  
  write.g3.file(summary2, out_path, 'optim.summary.stage2')
  
  for (i in names(stage2_params)){
    attr(stage2_params[[i]], 'summary') <- NULL
    write.g3.param(stage2_params[[i]],
                   out_path,
                   paste0('params.out.stage2.', i))
  }
  
  ## ------------ Final parameter set ------------------------------------------
  
  final_lik <- 
    parallel::mclapply(stage2_params, 
                       function(x){ g3_lik_out(r_model, x) }, 
                       mc.cores = parallel::detectCores()) 
  final_score <- 
    final_lik %>% 
    dplyr::bind_rows(.id = 'group') %>% 
    dplyr::group_by(.data$group) %>% 
    dplyr::summarise(s = sum(.data$value*.data$weight))
  
  final_params <- 
    stage2_params[[final_score[which.min(final_score$s), 'group'][[1]]]]
  
  write.g3.param(final_params,
                 out_path,
                 'final.params')
  
  #save(final_params, file = file.path(out.dir, 'final_params.Rdata'))
  
  return(final_params)
}

#' @export
g3_iterative_setup <- function(lik_out,
                               grouping = list()){
  
  if (!inherits(grouping, 'list')){
    stop('The groupings argument should be a list')
  }
  
  if (!all(unlist(grouping) %in% gsub('^.dist_[a-z]+_', '', lik_out$comp))){
    stop('The specified grouping do not all match component names')
  }
  
  ## Setup groups for re-weighting (i.e. components to be optimised together)
  ## Components not grouped will be NA
  gps <- 
    grouping %>% 
    purrr::map(~tibble::tibble(comp = .)) %>% 
    dplyr::bind_rows(.id = 'group') %>% 
    dplyr::bind_rows(tibble::tibble(comp = NA_character_,
                                    group = NA_character_))
  
  ldf <- 
    lik_out %>% 
    dplyr::mutate(param_name = paste0(.data$comp, '_weight'),
                  comp = gsub('^.dist_[a-z]+_', '', .data$comp)) %>% 
    dplyr::left_join(gps, by = 'comp') %>% 
    dplyr::mutate(init_weight = 1/.data$value,
                  group = purrr::map2(.data$group, .data$comp, 
                                      ~tidyr::replace_na(.x,.y)) %>% unlist())
  
  ## Update parameter template
  param <- attr(lik_out, 'param')
  param$value[ldf$param_name] <- ldf$init_weight
  
  out.params <- 
    ldf %>% 
    split(.$group) %>% 
    purrr::map(
      function(x){
        param$value[x$param_name] <- 
          x$init_weight*1e4
        param
      }
    )
  return(list(params = out.params, grouping = ldf[,c('comp', 'group')]))
}

#' @export
g3_lik_out <- function(model, param){
  
  if(!('data.frame' %in% class(param) & 
       sum(c("switch", "value","optimise" ) %in% names(param))==3)){
    stop('Error in param, expected data.frame')
  }
  
  res <- model(param$value)
  out <- attributes(res)
  nll <- res[[1]]
  
  lik.out <- 
    out[grep('dist_.+_obs__num',names(out), value = TRUE)] %>% 
    purrr::map(~sum(.>0)) %>% 
    purrr::map(~tibble::tibble(df = .)) %>% 
    dplyr::bind_rows(.id = 'comp') %>% 
    dplyr::mutate(comp = gsub('_obs__num', '', .data$comp)) %>% 
    dplyr::left_join(
      out[grep('nll_.dist_.+__num',names(out), value = TRUE)] %>% 
        purrr::map(sum) %>% 
        purrr::map(~tibble::tibble(value = .)) %>% 
        dplyr::bind_rows(.id = 'comp') %>% 
        dplyr::mutate(comp = gsub('nll_(.+)__num$', '\\1', .data$comp)),
      by = 'comp') %>% 
    dplyr::left_join(param %>% 
                       dplyr::select(comp = .data$switch, weight = .data$value) %>% 
                       dplyr::mutate(comp = gsub('_weight', '', .data$comp),
                              weight = unlist(.data$weight)),
                     by = 'comp')
  attr(lik.out, 'param') <- param
  attr(lik.out, 'actions') <- attr(model,'actions')
  attr(lik.out, 'model_out') <- out
  attr(lik.out, 'nll') <- nll
  return(lik.out)
}

# g3_iterative_make_fun <- function(tmb_model,param){
#   obj_fun <- g3_tmb_adfun(tmb_model,param)
#   attr(obj_fun,'param') <- param
# }

#' @export
g3_iterative_run <- function(param, 
                             tmb_model, 
                             opt_method = 'BFGS', 
                             use_parscale = TRUE){
  
  # Compile and generate TMB ADFun (see ?TMB::MakeADFun)
  obj_fun <- g3_tmb_adfun(tmb_model, param)
  
  # Sort out upper and lower
  if (opt_method == 'L-BFGS-B'){
    parhigh <- g3_tmb_upper(param)
    parlow <- g3_tmb_lower(param)
  }
  else{
    parhigh <- Inf
    parlow <- -Inf
  }
  
  ## Control list for optimisation
  opt_control <- list(
    trace = 2,
    maxit = 1000,
    reltol = .Machine$double.eps^2
  )
  if (use_parscale){
    opt_control <- c(opt_control, list(parscale = g3_tmb_parscale(param)))
  }
  
  ## Run optimiser
  fit.opt <- optim(g3_tmb_par(param),
                   obj_fun$fn,
                   obj_fun$gr,
                   method = opt_method,
                   lower = parlow,
                   upper = parhigh,
                   control = opt_control)
  
  ## Optimised parameters  
  p <- g3_tmb_relist(param,fit.opt$par)
  param$value[names(p)] <- p
  
  ## Add summary of input/output to an attribute of params
  attributes(param)$summary <- data.frame(opt_method = opt_method,
                                          max_iterations = opt_control$maxit,
                                          reltol = opt_control$reltol,
                                          function_calls = fit.opt$counts[1],
                                          gradient_calls = fit.opt$counts[2],
                                          convergence = ifelse(fit.opt$convergence == 0, TRUE, FALSE),
                                          score = fit.opt$value)
  
  
  return(param)
  
}

#' @export
g3_iterative_final <- function(lik_out_list){
  
  weights <- 
    lik_out_list %>% 
    dplyr::bind_rows(.id = 'group') %>% 
    dplyr::group_by(.data$comp) %>% 
    dplyr::filter(.data$value == min(.data$value)) %>% 
    dplyr::select(.data$comp, .data$df, .data$value) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(weight = ifelse(.data$value == 0, 0, .data$df/.data$value),
                  param_name = paste0(.data$comp, '_weight'))
  
  params <- 
    lik_out_list %>% 
    purrr::map(~attr(.,'param')) %>% 
    purrr::map(function(x){
      x$value[weights$param_name] <- 
        weights$weight
      x
    })
  return(params)
}

#' @export
tabulate_SS <- function(lik.out, grouping){
  
  group_list <- split(grouping, grouping$group)
  
  SS <- 
    lik.out %>% 
    dplyr::bind_rows(.id = '.group') %>% 
    dplyr::mutate(comp = gsub('(cdist|adist)_([A-Za-z]+)_(.+)', '\\3', .data$comp)) %>% 
    dplyr::select(.data$.group, .data$comp, .data$value) %>% 
    tidyr::pivot_wider(names_from = .data$comp, values_from = .data$value, names_sort = TRUE) %>% 
    dplyr::left_join(group_list %>% 
                       purrr::map(~tibble::tibble(.id = paste(.$comp, collapse = '.'))) %>% 
                       dplyr::bind_rows(.id = '.group'), by = '.group') %>% 
    dplyr::relocate(.data$.id) %>% 
    as.data.frame() 
  
  rownames(SS) <- SS$.group
  
  ## Normalise
  SS_norm <- SS
  for (group in names(group_list)){
    for (comp in group_list[[group]]$comp){
      SS_norm[,comp] <- SS_norm[,comp] / SS_norm[group, comp]
    }
  }
  
  return(list(SS = SS, SS_norm = SS_norm))
  
}

#' @export
g3_tmb_parscale <- function (parameters) {
  # Temporary location - should be in gadget3
  # Get all parameters we're thinking of optimising
  p <- parameters[
    !is.na(parameters[['parscale']]) &
      parameters$optimise, c('switch', 'value', 'parscale')]
  
  # Get the length of all values
  p$val_len <- vapply(p[['value']], length, integer(1))
  
  # Turn into a list with same dimensions as each value
  out <- structure(
    lapply(seq_len(nrow(p)), function (i) rep(p[i, 'parscale'], p[i, 'val_len'])),
    names = gadget3:::cpp_escape_varname(p$switch))
  
  # Unlist the result to condense list back to vector
  unlist(out)
}


