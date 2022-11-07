#' Perform multiple optimisation runs of a model, reweighting with each run
#'
#' @param gd Directory to store output
#' @param wgts Directory name within gd to store run outputs
#' @param model A G3 model, produced by g3_to_r()
#' @param params.in Initial parameters to use with the model
#' @param grouping List of component names to optmise together
#' @param use_parscale Logical indicating whether optim(control$parscale) should be used
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param shortcut If TRUE, weights for each component will be approximated and a final optimisation performed
#' @return Final set of parameters
#' @importFrom rlang .data
#' @importFrom stats optim
#' @export
g3_iterative <- function(gd, wgts = 'WGTS',
                         model, params.in, 
                         grouping = list(),
                         use_parscale = TRUE,
                         method = 'BFGS',
                         control = list(),
                         shortcut = FALSE){
  
  out_path <- file.path(gd, wgts)
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
  
  ## Compile the TMB model if the R version is passed in...
  if (inherits(model, 'g3_r')){
    model <- gadget3::g3_to_tmb(actions = attr(model, 'actions'))
  }
  
  ## ---------------------------------------------------------------------------
  ## INTERNAL FUNCTIONS
  ## ---------------------------------------------------------------------------
  
  ## Finds components that failed to run using the output of g3_lik_out
  find_failed <- function(dat){
    
    tmp <- sapply(dat, function(x){
      any(is.na(x$value))
    })
    
    if (any(tmp)) return(names(which(tmp)))
    else return(NULL)
    
  }
  
  ## ---------------------------------------------------------------------------
  ## APPROXIMATING WEIGHTS AND RUNNING THE OPTIMISATION
  ## ---------------------------------------------------------------------------
  
  ## Will write out approximated weights whether or not the shortcut is run
  approx_weights <- estimate_weights(model, params.in)
  
  if (shortcut){
    
    echo_message('##  RUNNING ITERATIVE SHORTCUT\n')
    
    ## First approximate weights and write to file
    approx_weights <- estimate_weights(model, params.in)
    write.g3.file(approx_weights, out_path, 'final.weights.shortcut')
    
    ## Merge estimated weights into parameter template
    params.in$value[match(approx_weights$comp, params.in$switch)] <- approx_weights$weight
    
    ## Run optimisation
    final_params <- g3_optim(model = model, 
                             params = params.in,
                             use_parscale = use_parscale,
                             method = method,
                             control = control)
    
    write.g3.param(final_params,
                   out_path,
                   'final.params.shortcut',
                   add_parscale = use_parscale)
    
  }else{
    
    ## -------------------------------------------------------------------------
    ## ITERATIVE REWEIGHTING
    ## -------------------------------------------------------------------------
    
    ## Compile and generate TMB ADFun (see ?TMB::MakeADFun)
    obj_fun <- gadget3::g3_tmb_adfun(model, tmb_param)
    
    ## -------------- Iterative re-weighting setup -------------------------------
    
    ## Run the R model to get the initial results
    init_weights <- g3_lik_out(model, params.in)
    
    if (is.na(attr(init_weights, 'nll'))){
      stop('The gadget model did not run')
    }
    
    write.g3.file(init_weights, out_path, 'lik.init')
    
    ## Setup the initial parameter files
    init_params <- g3_iterative_setup(init_weights[init_weights$weight > 0,], 
                                      grouping = grouping)
    
    ## Write groupings 
    write.g3.file(attr(init_params, 'grouping'), out_path, 'lik.groupings')
    
    ## Initial parameters (weights included in parameter template)
    for (i in names(init_params)){
      write.g3.param(init_params[[i]], 
                     out_path, 
                     paste0('params.init.stage1.', i),
                     add_parscale = use_parscale)
    }
    
    ## -------------- Run first stage of iterative re-weighting  -----------------
    
    echo_message('##  STAGE 1 OPTIMISATION\n')
    stage1_params <- parallel::mclapply(stats::setNames(names(init_params), 
                                                        names(init_params)),
                                        function(x){
                                          g3_optim(model = model,
                                                   params = init_params[[x]],
                                                   use_parscale = use_parscale,
                                                   method = method,
                                                   control = control,
                                                   print_status = TRUE,
                                                   print_id = x)
                                          },
                                        mc.cores = parallel::detectCores())
    
    
    ## Summary of optimisation settings and run details
    lapply(stage1_params, function(x) attr(x, 'summary')) %>% 
      dplyr::bind_rows(.id = 'group') %>% 
      write.g3.file(out_path, 'optim.summary.stage1')
    
    for (i in names(stage1_params)){
      attr(stage1_params[[i]], 'summary') <- NULL
      write.g3.param(stage1_params[[i]],
                     out_path,
                     paste0('params.out.stage1.', i),
                     add_parscale = use_parscale)
    }
    
    
    
    ## ------------ Update weights for second round of re-weighting --------------
    
    ## Update weights for second round of re-weighting
    lik_s1 <- parallel::mclapply(stage1_params, 
                                 function(x){ g3_lik_out(model, x) }, 
                                 mc.cores = parallel::detectCores())  
    
    ## Update the parameters
    int_params <- g3_update_weights(lik_s1, attr(init_params, 'grouping'))
    
    ## Identify any failed components
    failed_components <- find_failed(lik_s1)
    
    ## Update params and weights for failed components
    if (length(failed_components > 0)){
      
      ## Find the weights for the failed component
      bad_pars <- 
        attr(init_params, 'grouping') %>%
        dplyr::select(-comp) %>% 
        dplyr::rename(comp = param_name) %>% 
        dplyr::left_join(approx_weights, by = 'comp') %>% 
        dplyr::filter(group %in% failed_components) %>% 
        dplyr::select(comp, weight)   
      
      ## Merge the approximated weight into each component
      int_params <- 
        lapply(int_params, function(x){
          x$value[bad_pars$comp] <- bad_pars$weight
          return(x)
        })
      
      ## Now adjust the parameters
      for (i in seq_along(failed_components)){
        warning(paste0('## STAGE 1: optimisation for component ', i, ' failed, therefore corresponding weights are approximated using the shortcut method and optimised parameter values are taken from the initial values'))
        ## Identify the NA params
        na_params <- int_params[[i]][is.na(int_params[[i]]$value), 'switch']
        ## Merge values from previous param set
        int_params[[i]]$value[na_params] <- init_params[[i]]$value[na_params]
      }
    }
        
    for (i in names(int_params)){
      write.g3.param(int_params[[i]],
                     out_path,
                     paste0('params.init.stage2.', i),
                     add_parscale = use_parscale)
    }
    
    ## ----------- Second round of re-weighting ----------------------------------
    
    echo_message('\n##  STAGE 2 OPTIMISATION\n')
    stage2_params <- parallel::mclapply(stats::setNames(names(int_params), 
                                                        names(int_params)),
                                        function(x){
                                          g3_optim(model = model,
                                                   params = int_params[[x]],
                                                   use_parscale = use_parscale,
                                                   method = method,
                                                   control = control,
                                                   print_status = TRUE,
                                                   print_id = x)
                                          },
                                        mc.cores = parallel::detectCores())
    
    #save(stage2_params, file = file.path(out.dir, 'stage2_params.Rdata'))
    
    ## Summary of optimisation settings and run details
    lapply(stage2_params, function(x) attr(x, 'summary')) %>% 
      dplyr::bind_rows(.id = 'group') %>% 
      write.g3.file(out_path, 'optim.summary.stage2')
    
    for (i in names(stage2_params)){
      attr(stage2_params[[i]], 'summary') <- NULL
      write.g3.param(stage2_params[[i]],
                     out_path,
                     paste0('params.out.stage2.', i),
                     add_parscale = use_parscale)
    }
    
    ## ------------ Final parameter set ------------------------------------------
    
    final_lik <- 
      parallel::mclapply(stage2_params, 
                         function(x){ g3_lik_out(model, x) }, 
                         mc.cores = parallel::detectCores()) 
    
    final_score <- 
      final_lik %>% 
      dplyr::bind_rows(.id = 'group') %>% 
      dplyr::group_by(.data$group) %>% 
      dplyr::summarise(s = sum(.data$value*.data$weight))
    
    if (all(is.na(final_score$s))){
      warning('## STAGE 2: All optimisations failed')
      final_params <- NULL
    }
    else{
      
      if (any(is.na(final_score$s))){
        warning(paste0('## STAGE 2: optimisation failed for the following components: ', 
                       paste(final_score$group[is.na(final_score$s)], collapse = " ")))
      }  
    
      final_params <- 
        stage2_params[[final_score[which.min(final_score$s), 'group'][[1]]]]
      
      echo_message('\n##  Final parameters taken from component: ', final_score[which.min(final_score$s), 'group'][[1]])
      
      write.g3.param(final_params,
                     out_path,
                     'final.params',
                     add_parscale = use_parscale)
      
      ## Write the calculated and approximated weights to file
      approx_weights %>% 
        dplyr::select(comp, approx_weight = weight) %>% 
        dplyr::full_join(
          final_params %>% 
            dplyr::filter(grepl('_weight$', .$switch)) %>% 
            dplyr::select(comp = switch, weight = value) %>% 
            dplyr::mutate(weight = unlist(weight))
        ) %>% 
        write.g3.file(out_path, 'final.weights')
      
    }
  }
  save(final_params, file = file.path(out_path, 'final_params.Rdata'))
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
  
  ## Add groupings as an attributes
  attributes(out.params)$grouping <- ldf[,c('comp', 'group', 'param_name')]
  
  return(out.params)
}

#' @export
g3_lik_out <- function(model, param){
  
  if (inherits(model, 'g3_r')){
    model <- gadget3::g3_to_tmb(attr(model, 'actions'))
  }
  if (!inherits(param, 'data.frame')){
    stop('Error in param, expected a data.frame')
  }
  
  ## We only need reporting so build the objective function using type = 'Fun'
  adfun <- gadget3::g3_tmb_adfun(model, param, type = 'Fun')
  res <- adfun$report(adfun$par)
  
  lik.out <- 
    res[grep('dist_.+_obs__(wgt$|num$)', names(res), value = TRUE)] %>% 
    purrr::map(~sum(.>0)) %>% 
    purrr::map(~tibble::tibble(df = .)) %>% 
    dplyr::bind_rows(.id = 'comp') %>% 
    dplyr::mutate(comp = gsub('_obs__(wgt$|num$)', '', .data$comp)) %>% 
    dplyr::left_join(
      res[grep('nll_.dist_.+__(wgt$|num$)', names(res), value = TRUE)] %>% 
        purrr::map(sum) %>% 
        purrr::map(~tibble::tibble(value = .)) %>% 
        dplyr::bind_rows(.id = 'comp') %>% 
        dplyr::mutate(comp = gsub('nll_(.+)__(wgt$|num$)', '\\1', .data$comp)),
      by = 'comp') %>%
    dplyr::left_join(param %>% 
                       dplyr::select(comp = .data$switch, weight = .data$value) %>% 
                       dplyr::mutate(comp = gsub('_weight', '', .data$comp),
                                     weight = unlist(.data$weight)),
                     by = 'comp')
  
  attr(lik.out, 'param') <- param
  #attr(lik.out, 'actions') <- attr(model,'actions')
  #attr(lik.out, 'model_out') <- out
  attr(lik.out, 'nll') <- res$nll
  return(lik.out)
}

#' @export
g3_update_weights <- function(lik_out_list, grouping){
  
  ## Calculate new weights
  weights <- 
    lik_out_list %>% 
    dplyr::bind_rows(.id = 'group') %>% 
    dplyr::full_join(
      grouping %>% 
        dplyr::select(-comp) %>% 
        dplyr::rename(comp = param_name) %>% 
        dplyr::mutate(match = 1,
                      comp = gsub('_weight$', '', comp)),
      by = c('group', 'comp')
      ) %>% 
    dplyr::filter(!is.na(match)) %>% 
    dplyr::select(-match) %>% 
    dplyr::mutate(value = ifelse(.data$weight == 0, 0, .data$value),
                  weight = ifelse(.data$value == 0, 0, .data$df/.data$value),
                  param_name = paste0(.data$comp, '_weight')) 
  
  ## Merge into parameters
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
g3_iterative_final <- function(lik_out_list){
  
  weights <- 
    lik_out_list %>% 
    dplyr::bind_rows(.id = 'group') %>% 
    dplyr::mutate(value = ifelse(weight == 0, 0, value)) %>% 
    dplyr::group_by(.data$comp) %>% 
    dplyr::filter(.data$value == min(.data$value, na.rm = TRUE)) %>% 
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
