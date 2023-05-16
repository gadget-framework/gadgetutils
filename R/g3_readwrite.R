#' @title A wrapper for write.g3.file that formats a gadget3 TMB parameter data.frame for writing to file
#' @param params A TMB parameter data.frame
#' @param gd Directory to write the file
#' @param file_name Name of the file
#' @export
write.g3.param <- function(params, gd, file_name){
  
  ## Modify columns
  params$value <- unlist(params$value)
  params$optimise <- as.numeric(params$optimise)
  params$random <- as.numeric(params$random)
  params[params$type == '', 'type'] <- '.'
  
  write.g3.file(params, gd, file_name)
  
}

#' @title Writes an object to a file, typically just used within gadgetutils g3_ functions 
#' @param obj An object to write to file
#' @param gd Directory to write the file
#' @param file_name Name of the file
#' @export
write.g3.file <- function(obj, gd, file_name){
  
  if (tibble::is_tibble(obj)) obj <- as.data.frame(obj)
  
  outfile <- file(file.path(gd, file_name), 'w')
  utils::capture.output(print(format(obj, 
                                     #scientific = TRUE, 
                                     drop0trailing = TRUE), 
                              path = gd, 
                              row.names = FALSE,
                              max = 999999), 
                        file = outfile)
  close(outfile)
  
}

#' @title Reads a file from write.g3.file to an object. Now depreciated, using RData files for read/write.
#' @param gd Directory to write the file
#' @param file_name Name of the file
#' @export
read.g3.param <- function(gd, file.name){
  .Defunct()
}


