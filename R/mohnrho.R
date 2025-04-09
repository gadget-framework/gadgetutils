#' Generate parameterised g3a_initial_abund
#'
#' @param fit A g3_fit object
#' @param retrofit A list of g3_fit objects that differ in their terminal year, i.e. the result of a retrospective analysis
#' @param vars Which variables would you like to calculate Mohn's Rho for? They must be columns in fit$res.by.year
#' @param lag vector of lags for each variable in terms of years
#' @param aggregate_biomass calculate rho for biomass over stocks? Only works if 'total.biomass' is in vars
#' @return A data.frame with mohn's rho for each variable by stock
#'
#' @export
mohnrho <- function(fit, retrofit, 
                    vars = c('total.biomass', 'F', 'recruitment'), 
                    lag = c(0,1,1), aggregate_biomass = TRUE){
  
  ## Final year from fit
  final_year <- fit$res.by.year$year |> max()
  ## Create a lag data.frame to subset correct years
  lag_df <- data.frame(name = vars, lag = lag)
  ## Turn off aggregate biomass of 'total.biomass' is not in vars
  if (!('total.biomass' %in% vars)) aggregate_biomass <- FALSE
  
  retd <- do.call('rbind',
                  lapply(retrofit, function(x, final_year, lag_df){
                    
                    tmp <- x$res.by.year[,c('year', 'stock', vars)] 
                    if (aggregate_biomass){
                      tmp <- 
                        tmp |> 
                        dplyr::bind_rows(
                          tmp |> 
                            dplyr::group_by(year) |> 
                            dplyr::summarise(total.biomass = sum(total.biomass)) |> 
                            dplyr::mutate(stock = 'All')
                        )
                    }
                    
                    return(
                      tmp |> 
                        dplyr::mutate(peel = final_year - max(.data$year)) |> 
                        tidyr::pivot_longer(cols = vars) |> 
                        dplyr::right_join(lag_df, by = c('name')) |> 
                        dplyr::group_by(stock, peel, name) |> 
                        dplyr::filter(year == max(year)-lag) |> 
                        dplyr::select(-lag)
                    )
                    
                  }, final_year = final_year, lag_df = lag_df))
  
  fitd <- fit$res.by.year[,c('year', 'stock', vars)] 
  if (aggregate_biomass){
    fitd <- 
      fitd |> 
      dplyr::bind_rows(
        fitd |> 
          dplyr::group_by(year) |> 
          dplyr::summarise(total.biomass = sum(total.biomass)) |> 
          dplyr::mutate(stock = 'All')
      )
  }
  
  fitd <- 
    fitd |> 
    tidyr::pivot_longer(cols = vars, values_to = 'value0') |> 
    dplyr::right_join(retd, by = c('year','stock','name')) 
  
  out <- 
    fitd |> 
    dplyr::mutate(bias = (value - value0)/value0) |> 
    dplyr::group_by(stock, name) |> 
    dplyr::summarise(rho = mean(bias)) |> 
    tidyr::pivot_wider(names_from = name, values_from = rho)
  
  return(out)
  
}