#' Generate parameterised g3a_initial_abund
#'
#' @param imm A g3 stock object for immature specimens
#' @param imm A g3 stock object for mature specimens
#' @param comp_id Part of stock name to use for parameters, e.g. 'species' will share parameters with both mature/immature
#' @param mature Generate actions for mature (TRUE) or immature (FALSE) stock
#' @param init_mode One of 0 (initialised at equilibrium), 1 (Initial parameter per age group (across stocks)), 2 (Initial parameter per age group per stock)
#' @param bound_param Should parameters be normalised with g3 bounded() ?
#' @param exponentiate Should the initial abundance parameters be exponentiated?
#' @return g3a_initial_abund formula for use in g3a_initialconditions_normalparam()
#'
#' @export
init_abund <- function(imm,
                       mat,
                       comp_id = 'species',
                       mature = TRUE,
                       init_mode = 1,
                       bound_param = TRUE,
                       exponentiate = FALSE){
  
  stock <- if (mature) mat else imm
  
  ## Some defaults
  minage <- gadget3:::g3_step(~stock_with(imm, imm__minage))
  p_age <- 1
  
  ## MODE 0: initialised at equilibrium (using carrying capacity B0) assuming constant natual M (M)
  if (init_mode == 0){
    
    init <- 1
    m_table <- g3_stock_param(stock, comp_id, 'M', bound_param)
    prop_mat0 <- g3_stock_param(stock, comp_id, 'prop_mat0', bound_param)
    
    if(!mature){
      prop_mat0 <- gadget3:::f_substitute(~1-prop_mat0, list(prop_mat0 = prop_mat0))
    }
    
    ## Equilibrium age distribution
    init_scalar <- gadget3:::f_substitute(
      ~p0 * B0 * (1-exp(-1*M))/(1-exp(-1*maxage*M)),
      list(p0 = prop_mat0,
           B0 = g3_stock_param(stock, comp_id, 'B0', bound_param),
           M = m_table,
           maxage = gadget3:::g3_step(~stock_with(mat, mat__minage)))
    )
  }else{
    
    ## MODE 1: Initial parameter per age group (across stocks)  
    if (init_mode == 1){
      
      init <- g3_stock_table(list(imm, mat), comp_id, 'init', bound_param, exponentiate)
      init_scalar <- g3_stock_param(stock, comp_id, 'init.scalar', bound_param, exponentiate)
      m_table <- g3_stock_param(stock, 'full', 'M', bound_param)
      #m_table <- g3_stock_table(list(imm, mat), comp_id, 'M', bound_param)
    
      ## Proportion mature at age
      p_age <- 
        gadget3:::g3a_initial_ageprop(g3_stock_param(imm,
                                                    comp_id,
                                                    'mat_initial_alpha',
                                                    bound_param),
                                      g3_stock_param(imm,
                                                    comp_id, 
                                                    'mat_initial_a50',
                                                    bound_param))
      
      ## Invert for immature stock
      if(!mature){
        p_age <- gadget3:::f_substitute(~1-p_age, list(p_age = p_age))
      }
      
    }
    else{
      
      ## MODE 2: Initial parameter per age group per stock
      init <- g3_stock_table(stock, 'full', 'init', bound_param, exponentiate)
      init_scalar <- g3_stock_param(stock, 'full', 'init.scalar', bound_param, exponentiate)
      
      ## Minimum age taken from immature stock
      m_table <- g3_stock_param(stock, 'full', 'M', bound_param)
      
    }
  }
  
  ## Get the initial abundance
  gadget3:::g3a_initial_abund(
    scalar = init_scalar,
    init = init,
    M = m_table,
    init_F = g3_stock_param(stock, comp_id, 'init.F', bound_param),
    minage = minage,
    p_age = p_age)
  
}

#' Generate parameterised version of g3a_initial_vonb()
#'
#' @param stock A g3 stock object
#' @param id Part of the stock name to use in parameter name
#' @param bound_param Should this parameter be normalised with g3 bounded() ?
#' @return g3a_initial_vonb formula
#' @export
init_vonb <- function(stock, id = 'species', bound_param = TRUE){
  
  gadget3:::g3a_initial_vonb(
    recl = g3_stock_param(stock, id, 'recl', bound_param), 
    Linf = g3_stock_param(stock, id, 'Linf', bound_param),
    K = g3_stock_param(stock, id, 'K', bound_param))
  
}

#' Generate paramerterised stock renewal factor
#'
#' @param stock A g3 stock object
#' @param id Part of the stock name to use in parameter name
#' @param bound_param Should this parameter be normalised with g3 bounded() ?
#' @param exponentiate Should the renewal parameters be exponentiated?
#' @return "scalar * renew" formula, parameterised for given stock
#' @export
stock_renewal <- function(stock, 
                          id = 'species', 
                          bound_param = TRUE, 
                          exponentiate = FALSE){

  gadget3:::f_substitute(~scalar * renew,
                         list(scalar = g3_stock_param(stock,
                                                      id,
                                                      'rec.scalar', 
                                                      bound_param,
                                                      exponentiate),
                              renew = g3_year_table(stock,
                                                    id, 
                                                    'rec',
                                                    bound_param,
                                                    exponentiate)))
  
}

#' Generate parameterised initial conditions 
#'
#' @param stock A g3 stock object
#' @param id Part of the stock name to use in parameter name
#' @param parametric Is the initial conditions stddev parameterised, or a table by age?
#' @param bound_param Should this parameter be normalised with g3 bounded() ?
#' @return A formula suitable for g3a_initialconditions_normalparam()
#' @export
init_sd <- function(stock, id, parametric = TRUE, bound_param = TRUE){
  
  if (parametric){
    gadget3:::g3a_initial_sigma(
      g3_stock_param(stock, 'species', 'initial_sigma_alpha', FALSE),
      g3_stock_param(stock, 'species', 'initial_sigma_beta', FALSE),
      g3_stock_param(stock, 'species', 'initial_sigma_gamma', FALSE),
      init_vonb(stock, id, bound_param)
    )
  }
  else{
    g3_stock_table(stock, 'full', 'init.sd', bound_param = FALSE)
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
#' @param bound_param Should parameters be normalised with g3 bounded() ?
#' @param parametric_sd Is the initial conditions stddev parameterised, or a table by age?
#' @param exp_rec Should the recruitment parameters and scalar be exponentiated?
#' @param exp_init Should the initial abundance parameters and scalar be exponentiated?
#' @return A list of g3 actions
#' @export
model_actions <- function(imm, 
                          mat, 
                          mlgg = 15,
                          mature = TRUE, 
                          comp_id = 'species', 
                          init_mode = 1, 
                          bound_param = TRUE,
                          parametric_sd = FALSE,
                          exp_rec = FALSE,
                          exp_init = FALSE){
  
  stock <- if(mature) mat else imm

  ## Stock specific variables
  if (mature) output_stock <- list() else output_stock <- list(mat)
  
  ## ---------------------------------------------------------------------------
  ## SETUP ACTIONS
  stock_actions <- list(
    ## INITIAL CONDITIONS
    g3a_initialconditions_normalparam(stock,
                                      # NB: area & age factor together (gadget2 just multiplied them)
                                      # initial abundance at age is 1e4 * q
                                      factor_f =
                                        init_abund(imm, mat, comp_id, mature, init_mode, bound_param, exp_init),
                                      mean_f = init_vonb(stock, comp_id, bound_param),
                                      stddev_f = init_sd(stock, 
                                                         comp_id, 
                                                         parametric = parametric_sd, 
                                                         bound_param),
                                      alpha_f = g3_stock_param(stock, comp_id, 'walpha', bound_param),
                                      beta_f = g3_stock_param(stock, comp_id, 'wbeta', bound_param)),
    
    ## NATURAL MORALITY
    g3a_naturalmortality(stock, g3a_naturalmortality_exp({
      #bounded_table(stock, 'M', model_params, id = 'species')
      if (init_mode == 0){
        g3_stock_param(stock, comp_id, 'M', bound_param)
      }else{
        if (init_mode == 1){
          g3_stock_param(stock, 'full', 'M', bound_param)
          #g3_stock_table(list(imm, mat), comp_id, 'M', bound_param)
          }else{
            #g3_stock_table(stock, 'full', 'M', bound_param)
            g3_stock_param(stock, 'full', 'M', bound_param)
          }
      }
    })),
    
    ## AGING
    g3a_age(stock, output_stocks = output_stock)
  )
    
  if (!mature){
    
    stock_actions <- c(stock_actions, list(
      
      ## RENEWAL
      g3a_renewal_normalparam(imm,
                              factor_f = stock_renewal(imm, id = comp_id, bound_param, exp_rec),
                              mean_f = init_vonb(imm, comp_id, bound_param),
                              stddev_f = g3_stock_param(imm, comp_id, 'rec.sd', bound_param),
                              alpha_f = g3_stock_param(imm, comp_id, 'walpha', bound_param),
                              beta_f = g3_stock_param(imm, comp_id, 'wbeta', bound_param),
                              run_f = gadget3:::f_substitute(
                                ~cur_step == 1 && age == minage && cur_time > 0,
                                list(minage = gadget3:::g3_step(~stock_with(imm, imm__minage))))),
      
      ## GROWTH AND MATURATION
      g3a_growmature(imm,
                     impl_f = g3a_grow_impl_bbinom(
                       g3a_grow_lengthvbsimple(g3_stock_param(imm, comp_id, 'Linf', bound_param),
                                               {gadget3:::f_substitute(~(0.001 * K), 
                                                                       list(K = g3_stock_param(imm,
                                                                                               comp_id,
                                                                                               'K',
                                                                                               bound_param)))}),      
                       g3a_grow_weightsimple(g3_stock_param(imm, comp_id, 'walpha', bound_param),
                                             g3_stock_param(imm, comp_id, 'wbeta', bound_param)),   
                       beta_f = {gadget3:::f_substitute(~(10 * bbin), 
                                                        list(bbin = g3_stock_param(imm,
                                                                                   comp_id,
                                                                                   'bbin', 
                                                                                   bound_param)))},
                       maxlengthgroupgrowth = mlgg),
                     maturity_f = g3a_mature_continuous(
                       alpha = gadget3:::f_substitute(~(0.001 * exp(mat1)),
                                                      list(mat1 = g3_stock_param(imm,
                                                                                 comp_id,
                                                                                 'mat1',
                                                                                 bound_param))),
                       l50 = g3_stock_param(imm, comp_id, 'mat2', bound_param)
                     ),
                     output_stocks = list(mat),
                     transition_f = ~cur_time > 0),
      
      list()
      
    ))
  }
  else{
    
    stock_actions <- c(stock_actions, list(
      
      g3a_growmature(mat,
                     impl_f = g3a_grow_impl_bbinom(
                       g3a_grow_lengthvbsimple(g3_stock_param(mat, comp_id, 'Linf', bound_param),
                                               {gadget3:::f_substitute(~(0.001 * K), 
                                                                       list(K = g3_stock_param(mat,
                                                                                               comp_id,
                                                                                               'K',
                                                                                               bound_param)))}),       
                       g3a_grow_weightsimple(g3_stock_param(mat, comp_id, 'walpha', bound_param),
                                             g3_stock_param(mat, comp_id, 'wbeta', bound_param)),    
                       beta_f = {gadget3:::f_substitute(~(10 * bbin), 
                                                        list(bbin = g3_stock_param(mat,
                                                                                   comp_id,
                                                                                   'bbin', 
                                                                                   bound_param)))},
                       maxlengthgroupgrowth = mlgg)),
      list()
      
    ))
  }
  return(stock_actions)
}
