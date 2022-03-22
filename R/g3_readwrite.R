#' @export
write.g3.param <- function(params, gd, file_name){
  
  params$value <- unlist(params$value)
  params$optimise <- as.numeric(params$optimise)
  params$random <- as.numeric(params$random)
  params[params$type == '', 'type'] <- '.'
  
  write.g3.file(params, gd, file_name)
  
}

#' @export
write.g3.file <- function(obj, gd, file_name){
  
  if (tibble::is_tibble(obj)) obj <- as.data.frame(obj)
  
  outfile <- file(file.path(gd, file_name), 'w')
  utils::capture.output(print(format(obj, 
                                     scientific = FALSE, 
                                     drop0trailing = TRUE), 
                              path = gd, 
                              row.names = FALSE,
                              max = 999999), 
                        file = outfile)
  close(outfile)
  
}

