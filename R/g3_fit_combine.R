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

#' @title Combine multiple gadget.fit objects
#' @description This function combines multiple gadget.fit objects into a single fit object.
#' An additional column 'id' is added to each element of the resulting list to identify 
#' the individual gadget.fits
#' @param fit_list A list of multiple gadget.fit objects. 
#' @details
#' The 'id' column will be taken from names(fit_list). If fit_list has no names, the 'id'
#' column will be taken from the positions, e.g., 
#' fit_list = list(a = fit1, b = fit2) will be produces id's of 'a' and 'b', whereas 
#' fit_list = list(fit1, fit2) will produce id's of '1' and '2'
#' 
#' @returns A gadget.fit object
#' @export
bind_fit <- function(fit_list){
  # initialise an empty list
  out <- fit_list[[1]]
  out[] <- NA
  for(i in 1:length(out)){
    if(sum(class(fit_list[[1]][[names(out)[i]]]) == "data.frame") == 1){
      tmp <- sapply(fit_list,function(x){x[names(out)[i]]})
      names(tmp) <- names(fit_list) # replace with readable model names
      out[[i]] <- dplyr::bind_rows(tmp, .id="id")
    } else { NULL }}

  return(out)
}
