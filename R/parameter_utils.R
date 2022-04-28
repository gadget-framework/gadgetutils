#' Extract species name from stock
#'
#' @param stock A g3 stock object
#' @return Species name from stock
#' @export
stock_species <- function(stock){
  return(g3_stock_name(stock, 'species'))
}

#' Extract parts from stock name
#'
#' @param stock A g3 stock object
#' @param id Part to extract, or 'full' to return the stock name
#' @return Part from stock name
#' @export
#' @export
g3_stock_name <- function(stock, id = 'full'){
  stopifnot(gadget3:::g3_is_stock(stock))
  if(id[1] == 'full'){
    return(stock$name)
  } else if(setequal(id,intersect(id, names(stock$name_parts)))){
    return(paste(stock$name_parts[id],collapse = '_'))
  } else {
    stop(sprintf('id should be %s, supplied %s', 
                 paste(names(stock$name_parts),collapse = ', '),
                 id))
  }
}

#' Generate g3 code to reference a single parameter
#'
#' @param stock A g3 stock object
#' @param id Part of the stock name to use in parameter name
#' @param param_name Parameter name to append to the stock name
#' @param bound_param Should this parameter be normalised with g3 bounded() ?
#' @param exponentiate Should the parameter be exponentiated
#' @return formula to insert parameter into g3 model
#' @export
g3_stock_param <- function(stock, 
                           id = 'full', 
                           param_name, 
                           bound_param = FALSE,
                           exponentiate = FALSE){
  
  ## Checks
  stopifnot(gadget3:::g3_is_stock(stock))
  
  ## Add suffix to param name if exponentiating
  param_name <- paste0(param_name, ifelse(exponentiate, '_exp', ''))
  
  ## Parameter name
  stock_param <- paste(g3_stock_name(stock, id), param_name, sep = ".")
  
  ## Setup unbounded param references
  if (exponentiate) out <- gadget3:::f_substitute(~exp(g3_param(x)), list(x = stock_param))
  else out <- gadget3:::f_substitute(~g3_param(x), list(x = stock_param))
  
  ## Add lower and upper bounds if the parameter is to be bounded
  if (bound_param){
    
    out <- gadget3:::f_substitute(~bounded(x,
                                           g3_param(x_lower, optimise = FALSE),#, value = -9999),
                                           g3_param(x_upper, optimise = FALSE)),#, value = 9999)),
                                  list(x = out,
                                       x_lower = paste(stock_param, 'lower', sep = '.'),
                                       x_upper = paste(stock_param, 'upper', sep = '.')))
    
  }
  return(out)
}

#' Generate g3 code to reference a stock parameter, by stock age
#'
#' @param stock A g3 stock object or a list of g3 stock objects (if creating a table that includes the age ranges of multiple stocks)
#' @param id Part of the stock name to use in parameter name
#' @param param_name Parameter name to append to the stock name
#' @param bound_param Should this parameter be normalised with g3 bounded() ?
#' @param exponentiate Should the parameter be exponentiated
#' @return formula to insert parameter table into g3 model
#' @export
g3_stock_table <- function(stock, 
                           id = 'full', 
                           param_name, 
                           bound_param = FALSE,
                           exponentiate = FALSE){
  
  ## Check whether stock is a g3 stock
  is.stock <- gadget3:::g3_is_stock(stock)
  
  ## If not, assume list and check whether each element is a g3_stock
  if (!is.stock){
    lapply(stock, function(x){
      if (!gadget3:::g3_is_stock(x)){
        stop("The 'stock' argument should either be a g3_stock or a list of g3_stocks")
      }
    })
  }
  else stock <- list(stock)
  
  ## Check for unique id if using > 1 stock
  if (length(stock) > 1){
    tmp <- do.call('c', lapply(stock, function(x) g3_stock_name(x, id = id)))
    if (length(unique(tmp)) > 1){
      cat("A common id is required when creating a parameter table from multiple stocks, the following were found: \n", tmp, '\n')
      stop()
    }
  }
  
  ## Get min and max age formulas
  minages <- lapply(stock, function(x) gadget3:::g3_step(~stock_with(x, x__minage)))
  maxages <- lapply(stock, function(x) gadget3:::g3_step(~stock_with(x, x__maxage)))
  
  lage <- minages[[1]]
  uage <- maxages[[1]]
  
  ## Create 'min' and 'max' formulas for if there is more than one stock
  if (length(stock) > 1){
    
    for (i in 2:length(stock)){
      
      lage <- gadget3:::f_substitute(~min(x, y), list(x = lage,
                                                      y = minages[[i]]))
      uage <- gadget3:::f_substitute(~max(x, y), list(x = uage,
                                                      y = maxages[[i]]))
      
    }
  }
  
  ## Add suffix to param name if exponentiating
  param_name <- paste0(param_name, ifelse(exponentiate, '_exp', ''))
  
  ## Parameter name
  stock_param <- paste(g3_stock_name(stock[[1]], id = id), param_name, sep = '.')
  
  ## Setup unbounded parameter reference
  out <- gadget3:::f_substitute(~g3_param_table(x, 
                                                data.frame(age = seq(a0, a1))),
                                #  lower = -9999,
                                #  upper = 9999),
                                list(x = stock_param,
                                     a0 = lage,
                                     a1 = uage))
  
  ## Exponentiate the parameter
  if (exponentiate){
    out <- gadget3:::f_substitute(~exp(x), list(x = out))
  }
  
  ## Add lower and upper bounds if the parameter is to be bounded
  if (bound_param){
    
    out <- gadget3:::f_substitute(~bounded(x, 
                                           g3_param(x_lower, optimise = FALSE),#, value = -9999), 
                                           g3_param(x_upper, optimise = FALSE)),#, value = 9999)),
                                  list(x = out,
                                       x_lower = paste(stock_param, 'lower', sep = '.'),
                                       x_upper = paste(stock_param, 'upper', sep = '.')))
    
  }
  return(out)
}

#' Generate g3 code to reference a stock parameter, divided by model year
#'
#' @param stock A g3 stock object
#' @param id Part of the stock name to use in parameter name
#' @param param_name Parameter name to append to the stock name
#' @param bound_param Should this parameter be normalised with g3 bounded() ?
#' @param exponentiate Should the parameter be exponentiated
#' @return formula to insert parameter table into g3 model
#' @export
g3_year_table <- function(stock, 
                          id = 'full', 
                          param_name, 
                          bound_param = FALSE, 
                          exponentiate = FALSE,
                          random = FALSE){
  
  
  ## Checks
  stopifnot(gadget3:::g3_is_stock(stock))
  
  ## Add suffix to param name if exponentiating
  param_name <- paste0(param_name, ifelse(exponentiate, '_exp', ''))
  
  ## Parameter name
  stock_param <- paste(g3_stock_name(stock, id), param_name, sep = ".")
  
  ## Setup unbounded, unexponentiated param reference
  out <- gadget3:::f_substitute(~g3_param_table(x, 
                                                expand.grid(cur_year = seq(start_year, 
                                                                           end_year)),
                                                #lower = -9999,
                                                #upper = 9999,
                                                random = random),
                                list(x = stock_param,
                                     random = random))
  
  ## Exponentiate the parameter
  if (exponentiate){
    out <- gadget3:::f_substitute(~exp(x), list(x = out))
  } 
  
  ## Add lower and upper bounds if the parameter is to be bounded
  if (bound_param){
    
    out <- gadget3:::f_substitute(~bounded(x, 
                                           g3_param(x_lower, optimise = FALSE),#, value = -9999), 
                                           g3_param(x_upper, optimise = FALSE)),#, value = 9999)),
                                  list(x = out,
                                       x_lower = paste(stock_param, 'lower', sep = '.'),
                                       x_upper = paste(stock_param, 'upper', sep = '.'),
                                       random = random))  
    
  }
  return(out)
}

