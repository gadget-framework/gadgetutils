#' Wrapper for plot.gadget.fit that saves all diagnostic graphs to a directory
#'
#' @param fit A gadget fit object (can be from gadget2 or gadget3)
#' @param path Directory path for saving figures, defaults to file.path(getwd(), 'figs')
#' @param quiet Logical indicating whether to print messages about the plotting process (set to \code{FALSE} to suppress the messages
#' @param ... Additional arguments passed to \link[ggplot2]{ggsave}
#' @return List of results for each optimisation run
#' @export
gadget_plots <- function(fit, path = NULL, quiet = FALSE, ...){
  
  if (is.null(path)){
    path <- file.path(getwd(), 'figs')
  }
  if (!dir.exists(path)) dir.create(path)
  
  ## Survey indices
  if(!quiet) message("Plotting survey indices")
  tmp <- plot(fit)
  ggsave("surveyIndices.png", 
         plot = tmp, 
         path = path,
         bg = "white",
         ...)
  
  ## Catch distributions
  if(!quiet) message("Plotting catch distributions")
  c1 <- plot(fit, data="catchdist.fleets")
  for (i in 1:length(c1)){
    
    ggsave(paste0("catchdistribution_", names(c1)[[i]], ".png"), 
           plot = print(c1[[i]]),
           path = path,
           bg = "white",
           ...)
  }
  
  ## Growth
  if(!quiet) message("Plotting growth")
  if (any(grepl('aldist', names(c1)))){
    c1 <- plot(fit, data="catchdist.fleets", type="growth")
    for (i in 1:length(c1)){
      ggsave(paste0("ageLength_", names(c1)[[i]], ".png"), 
             plot = print(c1[[i]]), 
             path=path,
             bg = "white",
             ...)
    }  
  }
  
  
  ## Parameters
  #  tmp <- plot(fit, data="params", height=par("din")[1]*2, width=par("din")[2])
  #  ggsave("parameters.png", path=path, plot = tmp, height=par("din")[1]*2, width=par("din")[2]*1.25, bg = "white")
  
  
  ## Residuals
  if(!quiet) message("Plotting residuals")
  tmp <- plot(fit, data = "catchdist.fleets", type = "resid")
  ggsave("residuals.png", 
         plot = tmp, 
         path=path,
         bg = "white",
         ...)
  
  ## ICES
  if(!quiet) message("Plotting annual output (ICES plot)")
  p1 <- plot(fit,data='res.by.year',type='total') + theme(legend.position = 'none') 
  p2 <- plot(fit,data='res.by.year',type='F') + theme(legend.position = 'none') 
  p3 <- plot(fit,data = 'res.by.year',type='catch') + theme(legend.position = 'none') 
  p4 <- plot(fit, data='res.by.year',type='rec')
  
  png(file=file.path(path, "annual_output.png"), width=960, height=960)
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2) 
  dev.off()
  
  ## Maturity ogives
  tmp <- plot(fit, data="stockdist")
  ggsave("maturity.png", plot = tmp[[1]], path = path, bg = "white", ...)
  
  ## Liklihood  
  #plot(fit, data="summary")
  #ggsave("liklihood.png", path=path)
  
  ## Liklihood  
  #plot(fit, data="summary", type="weighted")
  #ggsave("liklihood_weighted.png", path=path)
  
  ## Suitability
  tmp <- plot(fit, data="suitability")
  ggsave("suitability.png", plot = tmp, path=path, bg = "white", ...)
  
}
