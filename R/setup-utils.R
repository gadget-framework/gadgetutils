
#' Parameterized initial conditions formula
#'
#' @param scalar formula substituted into initial conditions calculation
#' @param init formula substituted into initial conditions calculation
#' @param M formula substituted into initial conditions calculation
#' @param init_F formula substituted into initial conditions calculation
#' @param by_stock Controls how the 'scalar', 'init', and 'M' parameters are grouped
#' @param by_stock_f Controls how the 'init_F' parameter is grouped
#' @return g3a_initial_abund formula for use in g3a_initialconditions_normalparam()
#'
#' @export
g3a_initial_abund <- function(
    scalar = g3_parameterized('init', by_stock = by_stock),
    init = g3_parameterized('init.scalar', by_stock = by_stock, by_age = TRUE),
    M = g3_parameterized('M', by_stock = by_stock),
    init_F = g3_parameterized('init.F', by_stock = by_stock_f),
    by_stock = TRUE,
    by_stock_f = FALSE){
  
  gadget3:::f_substitute(
    ~scalar * init * exp(-1 * (M + init_F) * age),
    list(scalar = scalar,
         init = init,
         M = M,
         init_F = init_F)
  )
}

#' Generate parameterised g3a_initial_abund
#'
#' @param imm A g3 stock object for immature specimens
#' @param imm A g3 stock object for mature specimens
#' @param comp_id Part of stock name to use for parameters, e.g. 'species' will share parameters with both mature/immature
#' @param mature Generate actions for mature (TRUE) or immature (FALSE) stock
#' @param init_mode One of 0 (initialised at equilibrium), 1 (Initial parameter per age group (across stocks)), 2 (Initial parameter per age group per stock)
#' @param exp_init Logical, whether the initial parameters should be exponentiated
#' @param exp_init_f Logical, whether the parameter 'init.f' should be exponentiated
#' @param naturalmortality A g3 parameter for natural mortality
#' @return g3a_initial_abund formula for use in g3a_initialconditions_normalparam()
#'
#' @export
init_abund <- function(imm,
                       mat,
                       comp_id = 'species',
                       mature = TRUE,
                       init_mode = 2,
                       exp_init = FALSE,
                       exp_init_f = FALSE,
                       naturalmortality = g3_parameterized('M', by_stock = TRUE)){
  
  ## Checks
  stopifnot(gadget3:::g3_is_stock(imm))
  stopifnot(gadget3:::g3_is_stock(mat))
  
  g3a_initial_ageprop <- function(alpha, a50){
    gadget3:::f_substitute(
      ~bounded(-1*alpha*(age - a50),0,1),
      list(alpha = alpha, a50 = a50))
  }
  
  ## ---------------------------------------------------------------------------
  
  stock <- if (mature) mat else imm
  
  ## MODE 0: initialised at equilibrium (using carrying capacity B0) assuming constant natual M (M)
  if (init_mode == 0){
    
    ## Proportion mature
    prop_mat0 <- setup_g3_param('prop_mat0', by_stock = comp_id)
    
    if(!mature){
      prop_mat0 <- gadget3:::f_substitute(~1-prop_mat0, list(prop_mat0 = prop_mat0))
    }
    
    ## Equilibrium age distribution
    init_scalar <- gadget3:::f_substitute(~p0 * B0 * (1-exp(-1*M))/(1-exp(-1*maxage*M)),
                                          list(p0 = prop_mat0,
                                               B0 = setup_g3_param('B0', by_stock = comp_id),
                                               M = naturalmortality,
                                               maxage = gadget3:::g3_step(~stock_with(mat, mat__maxage))))
      
    out <- g3a_initial_abund(scalar = init_scalar,
                             init = 1,
                             M = naturalmortality,
                             init_F = setup_g3_param(name = 'init.F',
                                                     by_stock = comp_id,
                                                     exp_param = exp_init_f))
    
    
  }else{
    
    ## MODE 1: Initial parameter per age group (across stocks)  
    if (init_mode == 1){
      
      ## Proportion mature at age
      p_age <- g3a_initial_ageprop(setup_g3_param('mat_initial_alpha',
                                                  by_stock = comp_id),
                                   setup_g3_param('mat_initial_a50',
                                                  by_stock = comp_id))
      
      ## Invert for immature stock
      if(!mature){
        p_age <- gadget3:::f_substitute(~1-p_age, list(p_age = p_age))
      }
    
      out <- gadget3:::f_substitute(~scalar*init*exp(-1*(M+init_F)*(age-minage))*p_age,
                                    list(scalar = setup_g3_param('init.scalar', 
                                                                 by_stock = comp_id, 
                                                                 exp_param = exp_init),
                                         init = setup_g3_param(name = 'init', 
                                                               by_stock = comp_id, 
                                                               by_age = TRUE, 
                                                               exp_param = exp_init),
                                         M = naturalmortality,
                                         init_F = setup_g3_param(name = 'init.F',
                                                                 by_stock = comp_id,
                                                                 exp_param = exp_init_f),
                                         minage = gadget3:::g3_step(~stock_with(imm, imm__minage)),
                                         p_age = p_age))
      
      
    }
    else{
      
      ## MODE 2: Initial parameter per age group per stock
      out <- g3a_initial_abund(scalar = setup_g3_param(name = 'init.scalar',
                                                       by_stock = TRUE,
                                                       exp_param = exp_init),
                               init = setup_g3_param(name = 'init',
                                                     by_stock = TRUE,
                                                     by_age = TRUE,
                                                     exp_param = exp_init),
                               M = naturalmortality,
                               init_F = setup_g3_param(name = 'init.F',
                                                       by_stock = comp_id,
                                                       exp_param = exp_init_f))
      
    }
  }
  
  return(out)
  
}

#' Generate paramerterised stock renewal factor
#'
#' @param stock A g3 stock object
#' @param id Part of the stock name to use in parameter name
#' @param exponentiate Should the renewal parameters be exponentiated?
#' @return "scalar * renew" formula, parameterised for given stock
#' @export
stock_renewal <- function(stock, 
                          id = 'species', 
                          exponentiate = FALSE){
  
  ## String identifier for exponentiated parameters
  suffix <- ifelse(exponentiate, '_exp', '')
  
  g3_parameterized(name = paste0('rec', suffix),
                   by_stock = id,
                   by_year = TRUE,
                   exponentiate = exponentiate,
                   scale = 
                     g3_parameterized(name = paste0('rec.scalar', suffix),
                                      by_stock = id,
                                      exponentiate = exponentiate),
                   ifmissing = NaN)
  
}

#' Generate parameterised initial conditions 
#'
#' @param stock A g3 stock object
#' @param id Part of the stock name to use in parameter name
#' @param parametric Is the initial conditions stddev parameterised, or a table by age?
#' @return A formula suitable for g3a_initialconditions_normalparam()
#' @export
init_sd <- function(stock, id, parametric = FALSE){
  
  ## Helper from gadget3
  g3a_initial_sigma <- function(alpha, beta, gamma, mean_l){
    gadget3:::f_substitute(
      ~mean_l * ( alpha + beta/age + gamma * age),
      list(alpha = alpha,
           beta = beta,
           gamma = gamma,
           mean_l=mean_l)
    )
  }
  
  if (parametric){
    g3a_initial_sigma(
      g3_parameterized('initial_sigma_alpha', by_stock = id),
      g3_parameterized('initial_sigma_beta', by_stock = id),
      g3_parameterized('initial_sigma_gamma', by_stock = id),
      g3a_renewal_vonb(by_stock = id),
    )
  }
  else{
    g3_parameterized('init.sd', by_stock = TRUE, by_age = TRUE)
  }
}

#' Generate a standard set of model actions
#'
#' @param imm A g3 stock object for immature specimens
#' @param mat A g3 stock object for mature specimens
#' @param mlgg maxlengthgroupgrowth for growth of both mature/immature
#' @param mature Generate actions for mature (TRUE) or immature (FALSE) stock
#' @param comp_id Part of stock name to use for parameters, e.g. 'species' will share parameters with both mature/immature
#' @param init_mode 
#' @param parametric_sd Is the initial conditions stddev parameterised, or a table by age?
#' @param exp_params Which parameters should be exponentiated? exp_params is a vector of parameter names, possible parameters include: c('linf','k','bbin','recl','rec.sd','mat1','mat2','init','init.scalar','rec','rec.scalar','init.f','m','walpha','wbeta'). Note that is a scalar is exponentiated the annual values will be too, and vice versa.
#' @param tv_params Which parameters should be time-varying? tv_params is a vector of parameter names, possible time-varying parameters include: 'linf','k','walpha','beta','bbin','recl','rec.sd','mat1','mat2','m'
#' @return A list of g3 actions
#' @export
model_actions <- function(imm, 
                          mat, 
                          mlgg = 15,
                          mature = TRUE, 
                          comp_id = 'species', 
                          init_mode = 2, 
                          parametric_sd = FALSE,
                          exp_params = c(),
                          tv_params = c()){
  
  
  stock <- if(mature) mat else imm
  
  ## Stock specific variables
  if (mature) output_stock <- list() else output_stock <- list(mat)
  
  ## tv_params lookup to lower
  if (!is.null(tv_params)){ 
    param_list <- c('linf','k','walpha','wbeta','bbin','recl','rec.sd','mat1','mat2','m')
    tv_params <- casefold(tv_params)
    if (!all(tv_params %in% param_list)){
      stop(paste0("The following parameters are not currently available as time-varying: ", 
                  paste0(tv_params[!(tv_params %in% param_list)], collapse = ', ')))
    }  
  }
  
  ## Exponentiating parameters
  exp_init <- exp_rec <- FALSE
  
  if (!is.null(exp_params)){ 
    param_list <- c('linf','k','bbin','recl','rec.sd','mat1','mat2',
                    'init','init.scalar','rec','rec.scalar','init.f','m','walpha','wbeta')
    if (exp_params == 'all'){
      exp_params <- param_list
    }
    else{
      exp_params <- casefold(exp_params)
      
      ## If exponentiating init or rec, then the respective scalars should be too... and vice versa
      if (sum(c('init','init.scalar') %in% exp_params) > 0) exp_init <- TRUE
      if (sum(c('rec','rec.scalar') %in% exp_params) > 0) exp_rec <- TRUE
      
      if (!all(exp_params %in% param_list)){
        stop(paste0("The following parameters are not currently available to exponentiate: ", 
                    paste0(exp_params[!(exp_params %in% param_list)], collapse = ', ')))
      }
    }
  }
  
  ## Setup parameter references
  Linf <- setup_g3_param('Linf', comp_id, 'linf' %in% tv_params, 'linf' %in% exp_params)
  kk <- setup_g3_param('K', comp_id, 'k' %in% tv_params, 'k' %in% exp_params, scale = 0.001)
  walpha <- setup_g3_param('walpha', comp_id, 'walpha' %in% tv_params, 'walpha' %in% exp_params)
  wbeta <- setup_g3_param('wbeta', comp_id, 'wbeta' %in% tv_params, 'wbeta' %in% exp_params)
  bbin <- setup_g3_param('bbin', comp_id, 'bbin' %in% tv_params, 'bbin' %in% exp_params, scale = 10)
  recl <- setup_g3_param('recl', comp_id, 'recl' %in% tv_params, 'recl' %in% exp_params)
  recsd <- setup_g3_param('rec.sd', comp_id, 'rec.sd' %in% tv_params, 'rec.sd' %in% exp_params)
  mat_alpha <- setup_g3_param('mat1', comp_id, 'mat1' %in% tv_params, 'mat1' %in% exp_params, scale = 0.001)
  mat_l50 <- setup_g3_param('mat2', comp_id, 'mat2' %in% tv_params, 'mat2' %in% exp_params)
  natm <- setup_g3_param('M', TRUE, 'M' %in% tv_params, 'M' %in% exp_params)
  
  ## Create some variables
  initvonb <- g3a_renewal_vonb(Linf = Linf, K = kk, recl = recl)
  
  ## ---------------------------------------------------------------------------
  ## SETUP ACTIONS
  ## ---------------------------------------------------------------------------
  
  
  stock_actions <- list(
    ## INITIAL CONDITIONS
    g3a_initialconditions_normalparam(stock,
                                      # NB: area & age factor together (gadget2 just multiplied them)
                                      # initial abundance at age is 1e4 * q
                                      factor_f =
                                        init_abund(imm, mat, comp_id, mature, init_mode, 
                                                   exp_init,
                                                   'init.f' %in% exp_params,
                                                   natm),
                                      mean_f = initvonb,
                                      stddev_f = init_sd(stock, 
                                                         comp_id, 
                                                         parametric = parametric_sd),
                                      alpha_f = walpha,
                                      beta_f = wbeta),
    
    ## NATURAL MORALITY
    g3a_naturalmortality(stock, g3a_naturalmortality_exp(natm)),
    
    ## AGING
    g3a_age(stock, output_stocks = output_stock)
  )
  
  if (!mature){
    
    stock_actions <- c(stock_actions, list(
      
      ## RENEWAL
      g3a_renewal_normalparam(imm,
                              factor_f = stock_renewal(imm, id = list(imm, mat), exp_rec),
                              mean_f = initvonb,
                              stddev_f = recsd,
                              alpha_f = walpha,
                              beta_f = wbeta,
                              run_f = gadget3:::f_substitute(
                                ~cur_step == 1 && age == minage && cur_time > 0 && !cur_year_projection,
                                list(minage = gadget3:::g3_step(~stock_with(imm, imm__minage))))),
      
      ## GROWTH AND MATURATION
      g3a_growmature(imm,
                     impl_f = g3a_grow_impl_bbinom(
                       delta_len_f = g3a_grow_lengthvbsimple(Linf, kk),      
                       delta_wgt_f = g3a_grow_weightsimple(walpha, wbeta),   
                       beta_f = bbin,
                       maxlengthgroupgrowth = mlgg),
                     maturity_f = g3a_mature_continuous(
                       alpha = mat_alpha,
                       l50 = mat_l50),
                     output_stocks = list(mat),
                     transition_f = ~cur_time > 0),
      list()
    ))
  }
  else{
    
    stock_actions <- c(stock_actions, list(
      
      g3a_growmature(mat,
                     impl_f = g3a_grow_impl_bbinom(
                       delta_len_f = g3a_grow_lengthvbsimple(Linf, kk),       
                       delta_wgt_f = g3a_grow_weightsimple(walpha, wbeta),    
                       beta_f = bbin,
                       maxlengthgroupgrowth = mlgg)),
      list()
      
    ))
  }
  return(stock_actions)
}

## Helper for g3_parameterized that modifies the name of a parameter if it is being exponentiated
setup_g3_param <- function(name, 
                           by_stock,
                           by_age = FALSE,
                           tv_param = FALSE, 
                           exp_param = FALSE, ...){
  
  ## Suffix for exponentiated parameters
  suffix <- ifelse(exp_param, '_exp', '')
  
  g3_parameterized(paste0(name, suffix),
                   by_stock = by_stock,
                   by_age = by_age,
                   by_year = ifelse(tv_param, TRUE, FALSE),
                   exponentiate = ifelse(exp_param, TRUE, FALSE),
                   ...)
  
}
