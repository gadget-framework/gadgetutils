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
