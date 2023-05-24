#' Perform multiple optimisation runs of a model, reweighting with each run
#' 
#' An implementation of the iterative reweigthing of likelihood
#' components in gadget. It analyzes a given gadget model and, after
#' a series of optimisations where each likelihood component is
#' heavily weigthed, suggests a weigthing for the components based on
#' the respective variance.  If one (or more) components, other than
#' understocking and penalty, are 0 then the gadget optimisation with
#' the final weights will not be completed.
#'
#' In Taylor et. al an objective reweighting scheme for likelihood
#' components is described for cod in Icelandic waters. The authors
#' nota that the issue of component weighting has been discussed for
#' some time, as the data sources have different natural scales (e.g
#' g vs. kg) that should not affect the outcome. A simple heuristic,
#' where the weights are the inverse of the initial sums of squares
#' for the respective component resulting in an initials score equal
#' to the number of components, is therfor often used. This has the
#' intutitive advantage of all components being normalised. There is
#' however a drawback to this since the component scores, given the
#' initial parametrisation, are most likely not equally far from
#' their respective optima resulting in sub-optimal weighting.  The
#' iterative reweighting heuristic tackles this problem by optimising
#' each component separately in order to determine the lowest
#' possible value for each component. This is then used to determine
#' the final weights.  The resoning for this approach is as follows:
#' Conceptually the likelihood components can be thought of as
#' residual sums of squares, and as such their variance can be
#' esimated by dividing the SS by the degrees of freedom. The optimal
#' weighting strategy is the inverse of the variance.  Here the
#' iteration starts with assigning the inverse SS as the initial
#' weight, that is the initial score of each component when
#' multiplied with the weight is 1. Then an optimisation run for each
#' component with the intial score for that component set to
#' 10000. After the optimisation run the inverse of the resulting SS
#' is multiplied by the effective number of datapoints and used as
#' the final weight for that particular component.  The effective
#' number of datapoints is used as a proxy for the degrees of freedom
#' is determined from the number of non-zero datapoints. This is
#' viewed as satisfactory proxy when the dataset is large, but for
#' smaller datasets this could be a gross overestimate. In
#' particular, if the surveyindices are weigthed on their own while
#' the yearly recruitment is esimated they could be overfitted. If
#' there are two surveys within the year Taylor et. al suggest that
#' the corresponding indices from each survey are weigthed
#' simultaneously in order to make sure that there are at least two
#' measurement for each yearly recruit, this is done through
#' component grouping which is implemented. 
#' 
#' Component grouping is oftern applied in cases where overfitting is likely, 
#' and this can also happen with catch composition data as well as survey indices. 
#' 
#' In addition to grouping, a maximum weight can be assigned survey indices vie the
#' cv_floor setting. The cv_floor parameter sets the minimum of the estimated component 
#' variance and thus the maximum of the inverser variance. 
#' 
#' @title Iterative re-weighting
#' @param gd Directory to store output
#' @param wgts Directory name within gd to store run outputs
#' @param model A G3 model, produced by g3_to_r() or g3_tmb_adfun()
#' @param params.in Initial parameters to use with the model
#' @param grouping List of component names to optmise together
#' @param use_parscale Logical indicating whether optim(control$parscale) should be used
#' @param method The optimisation method, see \code{\link[stats]{optim}}
#' @param control List of control options for optim, see \code{\link[stats]{optim}}
#' @param shortcut If TRUE, weights for each component will be approximated and a final optimisation performed
#' @param cv_floor Minimum value of survey components (adist_surveyindices) as 1/\code{cv_floor}, applied prior to second stage of iterations. 
#' @param resume_final Logical value. If TRUE the re-weighting procedure starts at the second stage.
#' @param serial_compile g3_tmb_adfun will be run in serial mode (i.e., not in parallel), potentially helping with memory issues
#' @param mc.cores number of cores used, defaults to the number of available cores
#' @return Final set of parameters
#' @details Weights are calculated using inverse-variance weighting (\eqn{1/\sigma^2}), and as \code{1/pmax(variance, cv_floor)}, hence the minimum value for survey components is 1/\code{cv_floor}. Use smaller \code{cv_floor} values to increase the weight of survey components. 
#' @importFrom rlang .data
#' @importFrom stats optim
#' @export
g3_iterative <- function(gd, wgts = 'WGTS',
                          model, params.in, 
                          grouping = list(),
                          use_parscale = TRUE,
                          method = 'BFGS',
                          control = list(),
                          shortcut = FALSE,
                          cv_floor = 0,
                          resume_final = FALSE,
                          serial_compile = FALSE,
                          mc.cores = parallel::detectCores()){
  
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
  
  ## Merges the 'summary' attributes from a list of optimised parameters
  collect_summary <- function(param_list){
    
    tmp <- lapply(names(param_list), function(x){
      return(
        cbind(data.frame(group = x),
              attr(param_list[[x]], 'summary'),
              stringsAsFactors = FALSE)
      )
    })
    return(do.call('rbind', tmp))
    
  }
  
  ## Internal function to setup and run g3_optim
  run_g3_optim <- function(model, params,
                           use_parscale, method, control,
                           serial_compile, mc.cores){
    
    ## Compiling the model prior to g3_optim?
    if (serial_compile){
      echo_message('##  COMPILING MODEL AND CREATING ADFUN IN SERIAL\n')
      
      objfns <- lapply(stats::setNames(names(params), names(params)), function(x){
        print(x)
        return(gadget3::g3_tmb_adfun(model, params[[x]]))
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
      out <- parallel::mclapply(stats::setNames(names(params), names(params)), function(x){
        
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
      }, mc.cores = mc.cores)
    }
    return(out)
  }
  
  ## ---------------------------------------------------------------------------
  ## APPROXIMATING WEIGHTS AND RUNNING THE OPTIMISATION
  ## ---------------------------------------------------------------------------
  
  ## Will write out approximated weights whether or not the shortcut is run
  approx_weights <- estimate_weights(model, params.in)
  
  if (shortcut){
    
    echo_message('##  RUNNING ITERATIVE SHORTCUT\n')
    
    ## First approximate weights and write to file
    write.g3.file(approx_weights, out_path, 'weights.final.shortcut')
    
    ## Merge estimated weights into parameter template
    params.in$value[match(approx_weights$comp, params.in$switch)] <- approx_weights$weight
    
    ## Run optimisation
    params_final <- g3_optim(model = model, 
                             params = params.in,
                             use_parscale = use_parscale,
                             method = method,
                             control = control)
    
    write.g3.param(params_final,
                   out_path,
                   'params.final.shortcut')
    
  }else{
    
    ## -------------------------------------------------------------------------
    ## ITERATIVE REWEIGHTING
    ## -------------------------------------------------------------------------
    
    ## Compile and generate TMB ADFun (see ?TMB::MakeADFun)
    obj_fun <- gadget3::g3_tmb_adfun(model, params.in)
    
    ## -------------- Iterative re-weighting setup -------------------------------
    
    ## Skip first stage if resume_final = TRUE
    if (!resume_final){
      
      ## Run the R model to get the initial results
      lik_init <- g3_lik_out(model, params.in)
      
      if (is.na(attr(lik_init, 'nll'))){
        stop('The gadget model did not run')
      }
      
      write.g3.file(lik_init, out_path, 'lik.init')
      
      ## Setup the initial parameter files
      params_in_s1 <- g3_iterative_setup(lik_init[lik_init$weight > 0,],
                                         grouping = grouping)
      
      ## Write groupings 
      write.g3.file(attr(params_in_s1, 'grouping'), out_path, 'lik.groupings')
      
      ## Initial parameters (weights included in parameter template)
      for (i in names(params_in_s1)){
        write.g3.param(params_in_s1[[i]], 
                       out_path, 
                       paste0('params.in.stage1.', i)
        )
      }
      
      save(params_in_s1, file = file.path(out_path, 'params_in_s1.Rdata'))
      
      ## -------------- Run first stage of iterative re-weighting  -----------------
      
      echo_message('##  STAGE 1 OPTIMISATION\n')
      
      params_out_s1 <- run_g3_optim(model, params_in_s1,
                                    use_parscale, method, control,
                                    serial_compile, mc.cores)
      
      ## Check whether NULLs were passed out
      params_out_s1 <- check_null_params(params_out_s1, params_in_s1)
      
      ## Save
      save(params_out_s1, file = file.path(out_path, 'params_out_s1.Rdata'))
      
      ## Summary of optimisation settings and run details
      collect_summary(params_out_s1) %>% write.g3.file(out_path, 'optim.summary.stage1')
      
      ## Optimised parameters
      for (i in names(params_out_s1)){
        attr(params_out_s1[[i]], 'summary') <- NULL
        write.g3.param(params_out_s1[[i]],
                       out_path,
                       paste0('params.out.stage1.', i)
        )
      }
    }
    else{ 
      
      echo_message('##  RESUMING ITERATIVE REWEIGHTING AT STAGE 2\n')
      
      ## Resuming from saved stage 1 parameters
      if (file.exists(file.path(out_path, 'params_out_s1.Rdata'))){
        load(file = file.path(out_path, 'params_out_s1.Rdata'))
        load(file = file.path(out_path, 'params_in_s1.Rdata'))
      }
      else stop(paste("'resume_final' was TRUE, however, the directory",
                      out_path, "does not contain a 'params_out_s1.Rdata' file"))
      
    }
    
    ## ------------ Update weights for second round of re-weighting --------------
    
    ## Update weights for second round of re-weighting
    lik_s1 <- parallel::mclapply(params_out_s1, 
                                 function(x){ g3_lik_out(model, x) }, 
                                 mc.cores =  mc.cores)
    
    ## Calculate and write SS
    ss_s1 <- tabulate_SS(lik_s1, attr(params_in_s1, 'grouping'))
    write.g3.file(ss_s1$SS, out_path, 'lik.comp.score.stage1')
    write.g3.file(ss_s1$SS_norm, out_path, 'lik.comp.score.norm.stage1')
    
    ## Update the parameters
    params_in_s2 <- g3_update_weights(lik_s1, 
                                      attr(params_in_s1, 'grouping'),
                                      cv_floor)
    
    ## Identify any failed components
    failed_components <- find_failed(lik_s1)
    
    ## Update params and weights for failed components
    if (length(failed_components > 0)){
      
      ## Find the weights for the failed component
      bad_pars <- 
        attr(params_in_s1, 'grouping') %>%
        dplyr::select(-.data$comp) %>% 
        dplyr::rename(comp = .data$param_name) %>% 
        dplyr::left_join(approx_weights, by = 'comp') %>% 
        dplyr::filter(.data$group %in% failed_components) %>% 
        dplyr::select(.data$comp, .data$weight)   
      
      ## Merge the approximated weight into each component
      params_in_s2 <- 
        lapply(params_in_s2, function(x){
          x$value[bad_pars$comp] <- bad_pars$weight
          return(x)
        })
      
      ## Now adjust the parameters
      for (i in failed_components){
        warning(paste0('## STAGE 1: optimisation for component ', i, ' failed, therefore corresponding weights are approximated using the shortcut method and optimised parameter values are taken from the initial values'))
        ## Identify the NA params
        na_params <- params_in_s2[[i]][is.na(params_in_s2[[i]]$value), 'switch']
        ## Merge values from previous param set
        params_in_s2[[i]]$value[na_params] <- params_in_s1[[i]]$value[na_params]
      }
    }
    
    for (i in names(params_in_s2)){
      write.g3.param(params_in_s2[[i]],
                     out_path,
                     paste0('params.in.stage2.', i)
      )
    }  
    
    save(params_in_s2, file = file.path(out_path, 'params_in_s2.Rdata'))
    
    ## ----------- Second round of re-weighting ----------------------------------
    
    echo_message('\n##  STAGE 2 OPTIMISATION\n')
    
    params_out_s2 <- run_g3_optim(model, params_in_s2,
                                  use_parscale, method, control,
                                  serial_compile, mc.cores)
    
    ## Check whether NULLs were passed out
    params_out_s2 <- check_null_params(params_out_s2, params_in_s2)
    
    ## Save
    save(params_out_s2, file = file.path(out_path, 'params_out_s2.Rdata'))
    
    ## Summary of optimisation settings and run details
    collect_summary(params_out_s2) %>% write.g3.file(out_path, 'optim.summary.stage2')
    
    for (i in names(params_out_s2)){
      attr(params_out_s2[[i]], 'summary') <- NULL
      write.g3.param(params_out_s2[[i]],
                     out_path,
                     paste0('params.out.stage2.', i)
      )
    }
    
    ## ------------ Final parameter set ------------------------------------------
    
    lik_final <- 
      parallel::mclapply(params_out_s2, 
                         function(x){ g3_lik_out(model, x) }, 
                         mc.cores =  mc.cores) 
    
    ## Calculate and write SS
    ss_s2 <- tabulate_SS(lik_final, attr(params_in_s1, 'grouping'))
    write.g3.file(ss_s2$SS, out_path, 'lik.comp.score.stage2')
    write.g3.file(ss_s2$SS_norm, out_path, 'lik.comp.score.norm.stage2')
    
    
    final_score <- 
      lik_final %>% 
      dplyr::bind_rows(.id = 'group') %>% 
      dplyr::group_by(.data$group) %>% 
      dplyr::summarise(s = sum(.data$value*.data$weight))
    
    if (all(is.na(final_score$s))){
      warning('## STAGE 2: All optimisations failed')
      params_final <- NULL
    }
    else{
      
      if (any(is.na(final_score$s))){
        warning(paste0('## STAGE 2: optimisation failed for the following components: ', 
                       paste(final_score$group[is.na(final_score$s)], collapse = " ")))
      }  
      
      params_final <- 
        params_out_s2[[final_score[which.min(final_score$s), 'group'][[1]]]]
      
      echo_message('\n##  Final parameters taken from component: ', final_score[which.min(final_score$s), 'group'][[1]])
      
      write.g3.param(params_final,
                     out_path,
                     'params.final')
      
      ## Write the calculated and approximated weights to file
      approx_weights %>% 
        dplyr::select(.data$comp, approx_weight = .data$weight) %>% 
        dplyr::full_join(
          params_final %>% 
            dplyr::filter(grepl('_weight$', .data$switch)) %>% 
            dplyr::select(comp = .data$switch, weight = .data$value) %>% 
            dplyr::mutate(weight = unlist(.data$weight))
          , by = 'comp') %>% 
        write.g3.file(out_path, 'weights.final')
      
    }
  }
  save(params_final, file = file.path(out_path, 'params_final.Rdata'))
  return(params_final)  
}

#' @title Initial parameters for iterative re-weighting
#' @param lik_out A likelihood summary dataframe. The output of g3_lik_out(model, param)
#' @param grouping A list describing how to group likelihood components for iterative re-weighting
#' @return A list of G3 parameter dataframes, one for each group of likelihood components
#' @details Initial weights are calculated by taking the inverse SS and then multiplying these by 10000 for the components that are to be optimised 
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
    split(ldf,ldf$group) %>% 
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

#' @title Likelihood summary for a model run
#' @param model A G3 model, produced by g3_to_r() or g3_tmb_adfun()
#' @param param A G3 parameter dataframe
#' @return A dataframe containing the likelihood component name, degrees of freedom, NLL value, and weight
#' @details Runs the model, sums the number of data points (df) and nll's for each component (value), and extracts the corresponding weights from the input parameters (weight)
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

#' @title Updates weights for second stage of iterative re-weighting
#' @param lik_out_list A list of likelihood summarys. The output of g3_lik_out(model, param) for each group of likelihood components 
#' @param grouping A list describing how to group likelihood components for iterative re-weighting
#' @param cv_floor Minimum value of survey components (adist_surveyindices) as 1/\code{cv_floor}, applied prior to second stage of iterations. 
#' @return A list of G3 parameter dataframes, one for each group of likelihood components
#' @details This function updates weights for the second stage of iterative re-weighting. If the cv_floor is > 0, it will be applied here.
#' @export
g3_update_weights <- function(lik_out_list, grouping, cv_floor){
  
  ## Calculate new weights
  weights <- 
    lik_out_list %>% 
    dplyr::bind_rows(.id = 'group') %>% 
    dplyr::full_join(
      grouping %>% 
        dplyr::select(-.data$comp) %>% 
        dplyr::rename(comp = .data$param_name) %>% 
        dplyr::mutate(match = 1,
                      comp = gsub('_weight$', '', .data$comp)),
      by = c('group', 'comp')
      ) %>% 
    dplyr::filter(!is.na(match)) %>% 
    dplyr::select(-match) %>% 
    dplyr::mutate(value = ifelse(.data$weight == 0, 0, .data$value),
                  variance = .data$value/.data$df)
  
  ## Apply CV floor
  weights <- 
    weights %>% 
    dplyr::mutate(variance = ifelse(grepl('_surveyindices_', .data$comp),
                             pmax(.data$variance, cv_floor),
                             .data$variance))
  
  ## Update weights
  weights <- 
    weights %>% 
    dplyr::mutate(weight = ifelse(.data$value == 0, 0, 1/.data$variance),
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

#' @title Calculates likelihood component scores, including scores normalised by the minimum value
#' @param lik_out Likelihood summary for a model run. Output of g3_lik_out(model, param)
#' @param grouping A list describing how to group likelihood components for iterative re-weighting
#' @return A list with a table of likelihood component scores (SS) and a table of normalised likelihood component scores (SS_norm)
#' @export
tabulate_SS <- function(lik_out, grouping){
  
  group_list <- split(grouping, grouping$group)
  
  SS <- 
    lik_out %>% 
    dplyr::bind_rows(.id = 'group') %>% 
    dplyr::mutate(comp = gsub('(cdist|adist)_([A-Za-z]+)_(.+)', '\\3', .data$comp)) %>% 
    dplyr::select(.data$group, .data$comp, .data$value) %>% 
    tidyr::pivot_wider(names_from = .data$comp, values_from = .data$value, names_sort = TRUE) %>% 
    dplyr::left_join(group_list %>% 
                       purrr::map(~tibble::tibble(.id = paste(.$comp, collapse = '.'))) %>% 
                       dplyr::bind_rows(.id = 'group'), by = 'group') %>% 
    dplyr::relocate(.data$.id) %>% 
    as.data.frame() 
  
  rownames(SS) <- SS$group
  
  ## Normalise
  SS_norm <- SS
  for (group in names(group_list)){
    for (comp in group_list[[group]]$comp){
      SS_norm[,comp] <- SS_norm[,comp] / SS_norm[group, comp]
    }
  }
  
  return(list(SS = SS, SS_norm = SS_norm))
  
}
