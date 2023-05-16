#' Generate parameterised g3a_initial_abund
#'
#' @param fit A g3_fit object
#' @param retrofit A list of g3_fit objects that differ in their terminal year, i.e. the result of a retrospective analysis
#' @param vars Which variables would you like to calculate Mohn's Rho for? They must be columns in fit$res.by.year
#' @return A data.frame with mohn's rho for each variable by stock
#'
#' @export
mohnrho <- function(fit, retrofit, vars = c('total.biomass', 'F', 'recruitment')){
  
  ## Retro years
  retd <- do.call('rbind',
                  lapply(retrofit, function(x){
                    tmp <- x$res.by.year[,c('year', 'stock', vars)]
                    tmp <- tmp[tmp$year == max(tmp$year),]
                  }))
  
  ## Add column to maintain order
  retd$order <- 1:nrow(retd)
  
  ## Fit for corresponding years
  fitd <- merge(x = retd[, c('year', 'stock', 'order')], 
                y = fit$res.by.year[, c('year', 'stock', vars)], 
                by = c('year', 'stock'))
  
  ## Make sure the fit DF has the same order as the retro DF
  fitd <- fitd[fitd$order,]
  
  ## Calulate bias
  fitd[,vars] <- (retd[, vars] - fitd[, vars]) / fitd[, vars]
  
  ## And then mohn
  out <- stats::aggregate(fitd[,vars], list(stock = fitd$stock), FUN = mean)
  
  return(out)
  
}
