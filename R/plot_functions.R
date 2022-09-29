
#' @export
plot_bio <- function(fit, plot_total = FALSE, panelrow = FALSE){
  
  if (length(fit) == 1) fit <- fit[[1]]
  if (inherits(fit, 'gadget.fit')){
    
    pl <-
      ggplot2::ggplot(fit$res.by.year, ggplot2::aes(.data$year,
                                                    .data$total.biomass/1e6,
                                                    col=.data$stock)) + 
      ggplot2::geom_line() +
      ggplot2::labs(y="Biomass (in '000 tonnes)", x='Year',col='Stock')
    
    if (plot_total){
      pl <- 
        pl + ggplot2::geom_line(ggplot2::aes(.data$year,
                                             .data$total.biomass/1e6),
                                data = 
                                  fit$res.by.year %>% 
                                  dplyr::group_by(year) %>% 
                                  dplyr::summarise(total.biomass = sum(total.biomass)) %>% 
                                  dplyr::mutate(stock = 'Total') 
                                  
        )
      }
  }
  else{
    
    ## Extract components from each fit and bind together
    dat <- bind_fit_components(fit, component = 'res.by.year')
    
    ## Plot
    pl <- 
      ggplot2::ggplot(dat, 
                      ggplot2::aes(.data$year,
                                   .data$total.biomass/1e6,
                                   col=.data$id)) + 
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~stock, ncol = if (panelrow) 1 else NULL) + 
      ggplot2::labs(y="Biomass (in '000 tonnes)", x='Year',col='Model')
    
    if (plot_total){
      
      pl <- 
        pl + ggplot2::geom_line(ggplot2::aes(.data$year,
                                             .data$total.biomass/1e6),
                                data = 
                                  dat %>% 
                                  dplyr::group_by(id, year) %>% 
                                  dplyr::summarise(total.biomass = sum(total.biomass)) %>% 
                                  dplyr::mutate(stock = 'Total'))
      
    }
  }
  return(pl)
}

#' @export
plot_z <- function(fit, panelrow = FALSE){
  
  if (length(fit) == 1) fit <- fit[[1]]
  if (inherits(fit, 'gadget.fit')){
    
    pl <-
      ggplot2::ggplot(fit$res.by.year, ggplot2::aes(.data$year,
                                                    .data$F,
                                                    col=.data$stock)) + 
      ggplot2::geom_line() +
      ggplot2::labs(y='F', x='Year', col='Stock')
    
  }
  else{
    
    ## Extract components from each fit and bind together
    dat <- bind_fit_components(fit, component = 'res.by.year')
    
    ## Plot
    pl <- 
      ggplot2::ggplot(dat, 
                      ggplot2::aes(.data$year,
                                   .data$F,
                                   col=.data$id)) + 
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~stock, ncol = if (panelrow) 1 else NULL) +
      ggplot2::labs(y='F', x='Year', col='Model')
    
  }
  
  return(pl)
}

#' @export
plot_rec <- function(fit, panelrow = FALSE, stock = NULL){
  
  if (length(fit) == 1) fit <- fit[[1]]
  if (inherits(fit, 'gadget.fit')){
    
    if (is.null(stock)) stock <- unique(fit$res.by.year$stock)
    pdat <- 
      fit$res.by.year %>% 
      dplyr::filter(stock %in% stock) %>% 
      tidyr::drop_na()
    
    pl <-
      ggplot2::ggplot(pdat, 
                      ggplot2::aes(.data$year,
                                   .data$recruitment/1e6,
                                   col=.data$stock)) + 
      ggplot2::geom_line() +
      ggplot2::labs(y='Recruitment (in millions)', x='Year', col='Stock')
    
  }
  else{
    
    ## Extract components from each fit and bind together
    dat <- bind_fit_components(fit, component = 'res.by.year')
    
    if (is.null(stock)) stock <- unique(dat$stock)
    pdat <- 
      dat %>% 
      dplyr::filter(stock %in% stock) %>% 
      tidyr::drop_na()
    
    ## Plot
    pl <- 
      ggplot2::ggplot(pdat %>% dplyr::filter(stock %in% stock), 
                      ggplot2::aes(.data$year,
                                   .data$recruitment/1e6,
                                   col=.data$id)) + 
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~stock, ncol = if (panelrow) 1 else NULL) +
      ggplot2::labs(y='Recruitment (in millions)', x='Year', col='Model')
    
  }
  
  return(pl)
}

bind_fit_components <- function(fit_list, component){
  
  tmp <- lapply(fit_list, function(x, component){
    return(x[[component]])
  }, component = component)  
  out <- dplyr::bind_rows(tmp, .id = 'id')
  return(out)
    
}
