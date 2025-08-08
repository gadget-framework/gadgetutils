compare_fits <- function(fit1, fit2){
  
  out <- list()
  
  cat(" ------------- CHECKING NAMES ----------------\n\n")
  
  matvec <- all(sort(names(fit1)) == sort(names(fit2)))
  out$names <- matvec
  if (matvec) cat(" - Names matched\n\n")
  else{
    warning("Fit component names do not match")
    cat(" - Names do not match\n\n")
    # stop()
  }
  
  for (i in names(fit1)){
    
    cat(" ------------- COMPONENT ", i, "-------------\n\n")
    
    if (i == "params"){
      fit1$params$value <- unlist(fit1$params$value)
      fit2$params$value <- unlist(fit2$params$value)
    }
   
    tmp <- fit1[[i]]
    tmp2 <- fit2[[i]]
    
    if (is.null(tmp)){
      matvec <- is.null(tmp2)
      out[[paste(i, "_NULL")]] <- matvec
      if (matvec) cat(" - Components NULL in both fits\n\n")
      else{
        warning("Components not matching, one is null, the other is not")
        cat(" - Components not matching, one is null, the other is not")
      } 
      next
    }
    
    matvec <- all(dim(fit1[[i]]) == dim(fit2[[i]]))
    out[[paste0(i, "_dims")]] <- matvec
    if (matvec){
      cat(" - Dimensions match\n")
    }else{
      warning("Component dimensions do not match")
      cat(" - Component dimensions do not match\n")
    }
   
    out_tmp <- matrix(FALSE, nrow = 1, ncol = ncol(tmp))
    for (j in 1:ncol(tmp)){
      if (ncol(tmp2) < j) qq <- FALSE
      else qq <- all(tmp[[j]] == tmp2[[j]], na.rm = T)
      out_tmp[1,j] <- qq
      if (qq) cat(" - Column ", j, ", values match\n")
      else cat(" - Column ", j, ", values do not match\n")
    }
    out[[paste0(i, "_values")]] <- out_tmp
    cat("\n")
  }
  
  matvec <- any(!do.call("c", out))
  if (matvec)  cat("\nFITS DO NOT MATCH")
  else cat("\nFITS MATCH")
  
  out$FITS_MATCH <- !matvec
  
  return(out)
}
