#' @export
stock_species <- function(stock){
  return(g3_stock_name(stock, 'species'))
}

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

#' @export
g3_stock_param <- function(stock, id = 'full', param_name, bound_param = FALSE){
  
  ## Parameter name
  stock_param <- paste(g3_stock_name(stock, id), param_name, sep = ".")
  
  ## Unbounded parameter reference
  ## Note, using f_substitute so environment is attached (for subsequent f_substitute calls)
  if (!bound_param){
    out <- gadget3:::f_substitute(~g3_param(x),#, lower = -9999, upper = 9999), 
                                  list(x = stock_param))
  }
  else{
    ## Bounded parameter reference 
    ## Lower and upper bounds included as parameters that are not optimised
    out <- gadget3:::f_substitute(~bounded(g3_param(x),
                                           g3_param(x_lower, optimise = FALSE),#, value = -9999),
                                           g3_param(x_upper, optimise = FALSE)),#, value = 9999)),
                                  list(x = stock_param,
                                       x_lower = paste(stock_param, 'lower', sep = '.'),
                                       x_upper = paste(stock_param, 'upper', sep = '.')))
  }
  return(out)
}

#' @export
g3_stock_table <- function(stock, id = 'full', param_name, bound_param = FALSE){
  
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
  
  ## Parameter name
  stock_param <- paste(g3_stock_name(stock[[1]], id = id), param_name, sep = '.')
  
  
  if (!bound_param){
    ## Unbounded  
    out <- gadget3:::f_substitute(~g3_param_table(x, 
                                                  data.frame(age = seq(a0, a1))),
                                                #  lower = -9999,
                                                #  upper = 9999),
                                  list(x = stock_param,
                                       a0 = lage,
                                       a1 = uage))
  }
  else{
    ## Bounded
    out <- gadget3:::f_substitute(~bounded(g3_param_table(x,
                                                          data.frame(age = seq(a0, a1))), 
                                           g3_param(x_lower, optimise = FALSE),#, value = -9999), 
                                           g3_param(x_upper, optimise = FALSE)),#, value = 9999)),
                                  list(x = stock_param,
                                       a0 = lage,
                                       a1 = uage,
                                       x_lower = paste(stock_param, 'lower', sep = '.'),
                                       x_upper = paste(stock_param, 'upper', sep = '.')))
  }
  return(out)
}

#' @export
g3_year_table <- function(stock, id = 'full', param_name, bound_param = FALSE, random = FALSE){
  
  ## Checks
  stopifnot(gadget3:::g3_is_stock(stock))
  
  ## Parameter name
  stock_param <- paste(g3_stock_name(stock, id), param_name, sep = ".")
  
  if (!bound_param){
    ## Unbounded table
    out <- gadget3:::f_substitute(~g3_param_table(x, 
                                                  expand.grid(cur_year = seq(start_year, 
                                                                             end_year)),
                                                  #lower = -9999,
                                                  #upper = 9999,
                                                  random = random),
                                  list(x = stock_param,
                                       random = random))
  }
  else{
    ## Bounded
    out <- gadget3:::f_substitute(~bounded(g3_param_table(x,
                                                          expand.grid(cur_year = seq(start_year, 
                                                                                     end_year)),
                                                          random = random), 
                                           g3_param(x_lower, optimise = FALSE),#, value = -9999), 
                                           g3_param(x_upper, optimise = FALSE)),#, value = 9999)),
                                  list(x = stock_param,
                                       x_lower = paste(stock_param, 'lower', sep = '.'),
                                       x_upper = paste(stock_param, 'upper', sep = '.'),
                                       random = random))
  }
  return(out)
}

