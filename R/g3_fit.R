#' Build output summary suitable for feeding into gadgetplots
#'
#' @param model A G3 model built using \code{\link[gadget3]{g3_to_r}} or \code{\link[gadget3]{g3_to_tmb}}. The model must include \code{\link[gadget3]{g3a_report_detail}}.
#' @param params The fitted parameters post-optimisation
#' @inheritParams g3f_stock.recruitment
#' @param steps Which steps to include in the annual output? Vector of int.
#' @inheritParams g3f_f.by.year
#' @return List of data frames
#' @export
g3_fit <- function(model, 
                   params, 
                   rec_step = 1, 
                   steps = 1,
                   f_by_age = NULL){
  
  ## Checks
  stopifnot(is.list(params))
  
  if (inherits(model, "g3_r")) {
    if (is.data.frame(params)) params <- params$value
      if ("report_detail" %in% names(params)) {
          params$report_detail <- 1L
          tmp <- attributes(model(params))
          data_env <- environment(model)
      } else {
          tmp <- NULL
      }
  } else if (inherits(model, "g3_cpp")) {
      if (is.data.frame(params) && "report_detail" %in% params$switch) {
          params['report_detail', 'value'] <- 1L
          obj_fun <- gadget3::g3_tmb_adfun(model, params, type = 'Fun')
          tmp <- obj_fun$report(gadget3::g3_tmb_par(params))
          data_env <- as.environment(attr(model, 'model_data'))
      } else {
          tmp <- NULL
      }
  } else {
      stop("Unknown model type: ", class(model))
  }

  # If model is missing actions, then add them and try again
  if (is.null(tmp)) {
    
    ## Incorporate reporting action
    re_actions <- c(attr(model, 'actions'),
                    list(gadget3::g3a_report_detail(attr(model, 'actions'))))
    
    ## Compile model
    model <- gadget3::g3_to_r(re_actions)
    
    ## Run model
    if (is.data.frame(params)) params <- params$value
    params$report_detail <- 1L
    tmp <- attributes(model(params))
    data_env <- environment(model)
  }

  ## Add run information back to reporting output, if not already there
  if (is.null(tmp$model_params)) tmp$model_params <- params
  if (is.null(tmp$model_data)) tmp$model_data <- data_env

  return(g3_fit_inner(tmp, rec_step = rec_step, steps = steps, f_by_age = f_by_age))
} 

# Generate fit report from reporting output of g3_fit()
g3_fit_inner <- function(tmp,
                         rec_step = 1,
                         steps = 1,
                         f_by_age = NULL) {
  data_env <- tmp$model_data

  ## Returning parameters as they are put in for now
  ## so model(fit$params) will work without any fiddling
  params <- tmp$model_params
  out_params <- params
  
  ## Stock, predator, and fleet names
  sp_names <- stock_pred_names(tmp)
  fleet_names <- get_fleet_names(tmp)
  
  ##############################################################################
  
  ## Catch distributions
  catchdist <- g3f_catchdistribution(tmp)
  
  ## Sparse distributions
  sparsedist <- g3f_sparsedata(tmp, data_env = data_env)
  
  ## Survey or other indices 
  sidat <- g3f_sidat(tmp)
  
  ## Suitability
  suitability <- g3f_suitability(tmp, stock_names = sp_names$stocks)
  
  ## Likelihood
  likelihood <- g3f_likelihood(tmp, data_env = data_env, params = params)
  
  ## Stock-recruitment
  stock.recruitment <- 
    g3f_stock.recruitment(tmp, rec_step = rec_step)
  
  ## ---------------------------------------------------------------------------
  
  ## Stock and fleet reports
  stock_reports <- g3f_stock_reports(tmp)
  fleet_reports <- g3f_fleet_reports(tmp, stock_names = sp_names$stocks)
  
  ## Stock summaries
  stock.full <- g3f_stock.full(tmp, stock_reports = stock_reports)
  stock.std <- g3f_stock.std(tmp, stock_reports = stock_reports)
  
  ## Predation - fleets & biological
  stock.prey <- g3f_stock.prey(tmp, 
                               fleet_reports = fleet_reports,
                               stock_std = stock.std,
                               fleet_names = fleet_names)
  
  predator.prey <- g3f_predator.prey(tmp,
                                     fleet_reports = fleet_reports,
                                     stock_reports = stock_reports,
                                     stock_full = stock.full,
                                     suitability = suitability)
  
  fleet.info <- g3f_fleet.info(tmp, predator_prey = predator.prey)
  
  ## Predation - fleets
  harv.suit <- g3f_harv.suit(tmp, 
                             predator_prey = predator.prey,
                             fleet_names = fleet_names,
                             stock_full = stock.full)
  
  f.by.year <- g3f_f.by.year(tmp, 
                             stock_prey = stock.prey, 
                             stock_full = stock.full,
                             f_by_age = f_by_age)
  
  ## Annual output
  res.by.year <- g3f_res.by.year(tmp, 
                                 steps = steps,
                                 stock_full = stock.full,
                                 harv_suit = harv.suit,
                                 f_by_year = f.by.year,
                                 stock_recruitment = stock.recruitment)
  
  out <- list(
    sidat = sidat,
    stockdist = catchdist$stockdist,
    catchdist.fleets = catchdist$catchdist.fleets,
    sparsedist = sparsedist,
    #suitability = predator.prey %>% 
    #  dplyr::select(.data$year,.data$step,area=.data$area, stock=.data$prey,
    #                fleet=.data$predator,.data$length,.data$suit) %>% 
    #  dplyr::ungroup(),
    suitability = suitability |> 
      add_missing_columns(
        list(stock = NA_character_, fleet = NA_character_,
             year = NA_real_, step = NA_real_, area = NA_character_,
             age = NA_integer_, length = NA_real_, predator_length = NA_character_)
        ),
    likelihood = likelihood,
    stock.prey = stock.prey,
    stock.std = stock.std,
    stock.full = stock.full,
    fleet.info = fleet.info,
    predator.prey = predator.prey,
    stock.recruitment = stock.recruitment,
    res.by.year = res.by.year,
    params = out_params,
    score = data.frame(nll = tmp$nll)
  )
  
  ## Add 'summary' attribute if it exists
  if ('summary' %in% names(attributes(out_params))){
    attributes(out)$summary <- attr(out_params, 'summary')
  }
  
  class(out) <- c('gadget.fit',class(out))
  return(out)
  
}


#' @title Extracts year and step columns
#' @description Gadget3 reports have a time column in the format 'YEAR-STEP'. This function creates a 'year' and 'step' column by splitting the existing 'time' column. This is useful for processing gagdet3 reports.
#' @param data dataframe of gadget3 reports. Must include a 'time' column.
#' @export
extract_year_step <- function(data){
  if (!("time" %in% names(data))) return(data)
  
  if (any(!grepl('-', data$time))){
    data[!grepl('-', data$time), 'time'][[1]] <- 
      paste0(data[!grepl('-', data$time), 'time'][[1]], '-0')
  }
  
  data <- 
    within(data, {
    step = gsub('^(\\d+)-(\\d+)$', "\\2", time) |> as.numeric();
    year = gsub('^(\\d+)-(\\d+)$', "\\1", time) |> as.numeric()
  })
  
  return(data[,names(data) != "time"])
  
}

## -----------------------------------------------------------------------------
## Internal functions used by g3_fit
## -----------------------------------------------------------------------------

split_length <- function(data){
  
  if (!('length' %in% names(data))) return(data)
  
  data <- 
    data |> 
    within( {
      upper = gsub('(.+):(.+)', '\\2', length) |> as.numeric();
      lower = gsub('(.+):(.+)', '\\1', length) |> as.numeric()
    } )
  
  return(data)
  
}

replace_inf <- function(data, group_col = NULL){
  
  if (!('upper' %in% names(data))) return(data)
  
  inner_fun <- function(x){
    ind <- is.infinite(x$upper)
    if (any(ind)){
      replacement <- abs(diff(sort(unique(x$upper[!ind]), decreasing = TRUE)[1:2]))
      if (!is.na(replacement)) x$upper[ind] <- x$lower[ind] + replacement
    } 
    return(x)
  }
  
  if (is.null(group_col)) return(inner_fun(data))
  else {
    return(
      do.call("rbind", lapply(split(data, data[[group_col]]), FUN = inner_fun))
    )
  }
  
  return(data)
}

## Need to check what happens if no predation
stock_pred_names <- function(reports, stock_names = NULL){
  
  inner_fun <- function(reports, re, pos){
    
    tmp_names <- regmatches(names(reports), regexec(re, names(reports)))
    tmp_names <- vapply(tmp_names, function(x){
      if (length(x) == pos) x[[pos]] else as.character(NA)
    }, as.character(1))
    return(
      unique(tmp_names[!is.na(tmp_names)])
    )
    
  }
  
  if (is.null(stock_names)) 
    stock_names <- inner_fun(reports, "^(?:dstart|detail)_(.+)__num$", 2)
  
  pred_names <- 
    inner_fun(reports,
              paste0("^detail",
                     "_(\\Q", paste(stock_names, collapse = "\\E|\\Q"), "\\E)",
                     "_(.+)",
                     "__cons"), 3)
  if (length(pred_names) == 0) pred_names <- NA_character_
  
  return(list(stocks = stock_names, preds = pred_names))
  
}

## Helper to add columns
add_missing_columns <- function(data, columns){
  
  for (i in names(columns)){
    if (!i %in% names(data))  data[[i]] <- rep(columns[[i]], length.out = nrow(data))
  }
  return(data)
}

## Converts all reports matching 're' into a single long data frame. 'labels'
## is a named vector of replacements applied to the report names via 're',
## e.g. c(stock = "\\2") adds a stock column. Columns missing from some reports
## (e.g. predator_length for fleets) are filled with NA
reports_to_df <- function(reports, re, responseName = "Freq", labels = NULL){
  
  report_names <- grep(re, names(reports), value = TRUE)
  if (length(report_names) == 0) return(NULL)
  
  out <- lapply(report_names, function(x){
    df <- as.data.frame.table(reports[[x]],
                              stringsAsFactors = FALSE,
                              responseName = responseName)
    for (i in names(labels)) df[[i]] <- gsub(re, labels[[i]], x)
    return(df)
  })
  
  return(rbind_fill(out))
}

## rbind a list of data frames, filling columns missing from some with NA
rbind_fill <- function(data_list){
  
  data_list <- data_list[!vapply(data_list, is.null, logical(1))]
  if (length(data_list) == 0) return(NULL)
  
  cols <- unique(unlist(lapply(data_list, names)))
  out <- do.call("rbind", lapply(data_list, function(x){
    x[setdiff(cols, names(x))] <- NA
    return(x[cols])
  }))
  row.names(out) <- NULL
  
  return(out)
}

## Applies FUN to 'cols' within each combination of the 'by' columns
## (base R equivalent of dplyr::group_by() |> dplyr::summarise()). 'by' columns
## not present in data are ignored
aggregate_by <- function(data, by, cols, FUN = sum, ...){
  
  by <- intersect(by, names(data))
  if (nrow(data) == 0) return(data[, c(by, cols), drop = FALSE])
  
  return(
    stats::aggregate(data[cols], by = data[by], FUN = FUN, ...)
  )
}

## Orders a data frame by the 'by' columns (those present) and resets row names
order_by <- function(data, by){
  
  by <- intersect(by, names(data))
  data <- data[do.call(order, unname(as.list(data[by]))), , drop = FALSE]
  row.names(data) <- NULL
  
  return(data)
}

## Length of each step as a proportion of the year, e.g. 0.25 for a 3 month step
get_step_size <- function(reports, step){
  return(reports$step_lengths[step] / sum(reports$step_lengths))
}

## Mortality rate (per year) from the numbers consumed in a step of length 'step_size' (proportion of year).
## Mortality is 0 if nothing was consumed (including when there is nothing to consume). Any other
## non-finite values (missing numbers or consumption, or consumption >= numbers available) are
## kept, with a warning, as they indicate a problem with the model or the reports
calc_mortality <- function(consumed, number, step_size){
  out <- -log(1 - consumed / number) / step_size
  out[consumed %in% 0] <- 0
  if (any(!is.finite(out))){
    warning(sum(!is.finite(out)), " mortality value(s) are NA, NaN or Inf, due to missing ",
            "numbers or consumption, or consumption greater than or equal to the numbers available")
  }
  return(out)
}

## Annual fishing mortality averaged over ages (Fbar) from stock.prey. The catch (C) and
## numbers (N) are pooled within the 'by' columns (and year). 'data' must only contain the
## ages to be averaged over. F in a step is -log(1 - C/N), and the annual F is the sum over steps.
##   method = "ices":   F at age and step, summed over steps, then the unweighted mean over ages
##   method = "gadget": F by step with C and N pooled over ages, summed over steps
calc_fbar <- function(data, by = NULL, method = c("ices", "gadget")){
  
  method <- match.arg(method)
  cols <- c("number_consumed", "number")
  
  if (method == "ices"){
    fbar <- aggregate_by(data, by = c(by, "year", "step", "age"), cols = cols)
    fbar$F <- calc_mortality(fbar$number_consumed, fbar$number, step_size = 1)
    fbar <- 
      aggregate_by(fbar, by = c(by, "year", "age"), cols = "F") |> 
      aggregate_by(by = c(by, "year"), cols = "F", FUN = mean)
  } else {
    fbar <- aggregate_by(data, by = c(by, "year", "step"), cols = cols)
    fbar$F <- calc_mortality(fbar$number_consumed, fbar$number, step_size = 1)
    fbar <- aggregate_by(fbar, by = c(by, "year"), cols = "F")
  }
  
  return(fbar[c(by, "year", "F")])
}

## Predators that are not stocks
get_fleet_names <- function(reports){
  sp_names <- stock_pred_names(reports)
  return(sp_names$preds[!(sp_names$preds %in% sp_names$stocks)])
}
