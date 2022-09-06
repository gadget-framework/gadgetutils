#' @export
write.g3.param <- function(params, gd, file_name, add_parscale = TRUE){
  
  ## Add parscale if used
  if (add_parscale) params <- g3_add_parscale(params)
  
  ## Modify columns
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

#' @export
read.g3.param <- function(gd, file.name){
  
  path <- file.path(gd, file.name)
  if (!file.exists(path)){
    stop(paste0('The file: ', path, ' does not exist'))
  }
  
  params <- utils::read.table(path, header = TRUE)
  
  ## Modify columns
  params[params$type == '.', 'type'] <- ''
  params$optimise <- as.logical(params$optimise)
  params$random <- as.logical(params$random)
  
  ## Convert value columns to list
  params$value <- as.list(params$value)
  names(params$value) <- params$switch
  
  row.names(params) <- params$switch
  return(params)
  
}

#' @export
g3_add_parscale <- function(parameters){
  par <- g3_tmb_parscale(parameters)
  out <- utils::relist(par, unclass(parameters$value[parameters$optimise]))
  parameters$parscale[match(names(out), parameters$switch)] <- out
  return(parameters)
}
