# Functions that combine g3 fit objects

#' @title Combine individual elements from multiple gadget.fit objects
#' @description This function extracts a specific element from multiple gadget.fit objects
#' and combines them into a single dataframe. The output has an additional column 'id' 
#' that distinguishes each input gadget.fit object
#' @param fit_list A list of multiple gadget.fit objects. 
#' @param component A string specifying which list element to extract (ie, see names(fit))
#' @details
#' The 'id' column will be taken from names(fit_list). If fit_list has no names, the 'id'
#' column will be taken from the positions, e.g., 
#' fit_list = list(a = fit1, b = fit2) will be produces id's of 'a' and 'b', whereas 
#' fit_list = list(fit1, fit2) will produce id's of '1' and '2'
#' 
#' @returns A data frame that contains an 'id' column
#' @export
bind_fit_components <- function(fit_list, component){
  
  tmp <- lapply(fit_list, function(x, component){
    return(x[[component]])
  }, component = component)
  out <- dplyr::bind_rows(tmp, .id = 'id')
  return(out)
  
}
