#' @export
write.tmb.param <- function(params, dir, file.name){
  
  params$value <- unlist(params$value)
  params$optimise <- as.numeric(params$optimise)
  params$random <- as.numeric(params$random)
  
  write.g3.file(params, dir, file.name)
  
}

#' @export
write.g3.file <- function(obj, dir, file_name){
  
  if (tibble::is_tibble(obj)) obj <- as.data.frame(obj)
  
  outfile <- file(file.path(dir, file_name), 'w')
  utils::capture.output(print(format(obj, 
                                     scientific = FALSE, 
                                     drop0trailing = TRUE), 
                              path = dir, 
                              row.names = FALSE,
                              max = 999999), 
                        file = outfile)
  close(outfile)
  
}

