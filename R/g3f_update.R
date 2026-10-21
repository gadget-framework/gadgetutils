#' Adds recruitment at a given age to the annual output of a gadget.fit object
#'
#' @description Recruitment in \code{fit$stock.recruitment} (and the 'recruitment' column of \code{fit$res.by.year})
#' is the number of individuals entering a stock through renewal or spawning. This function instead
#' defines recruitment as the number of individuals of a given age, e.g. when they enter the fishery,
#' and adds it to \code{fit$res.by.year} as an additional column named \code{paste0("rec", age)}.
#' The numbers are those at the end of the step. As \code{fit$stock.std} contains the numbers at the
#' beginning of each step, the numbers at the end of a step are taken from the following step in
#' \code{fit$stock.std} (step 1 of the following year for the last step of the year, so there are no
#' numbers for the last step of the final year).
#' \code{fit$stock.recruitment} is not altered.
#'
#' @note Ageing occurs at the end of the last step of the year, so the numbers at the end of the last
#'   step are after ageing. If the recruitment step is the last step of the year, \code{age} should be
#'   one year older than the age of the recruits during that year, e.g. for recruits that are age 1
#'   during year y, use \code{age = 2} to get their numbers at the end of the last step of year y. The
#'   year in \code{fit$res.by.year} is year y. For all other steps, \code{age} is the age during the year.
#'
#' @param fit A gadget.fit object, see \code{\link{g3_fit}}
#' @param age The age of recruitment. Either an integer applied to all stocks, or a named list
#'   with an age for each stock, e.g. \code{list(imm1 = 1, imm2 = 5)}. Stocks not included in the list will be NA.
#' @param step The step at the end of which the recruitment numbers are taken. Either NULL, an integer applied to all stocks,
#'   or a named list with a step for each stock, e.g. \code{list(imm1 = 1, imm2 = 2)}. If NULL, or a stock
#'   is not included in the list, the numbers are taken from the same step(s) as \code{fit$res.by.year}
#'   (the 'steps' argument of \code{\link{g3_fit}}).
#' @return The gadget.fit object with additional column(s) in \code{fit$res.by.year}. As columns are
#'   named by age, stocks recruiting at different ages will have their recruitment in different columns.
#' @examples
#' \dontrun{
#' fit <- g3_fit(model, params)
#' fit <- g3_fit.update_rec(fit, age = list(imm1 = 1, imm2 = 5), step = list(imm1 = 1, imm2 = 2))
#' }
#' @export
g3_fit.update_rec <- function(fit, age, step = NULL){
  
  stopifnot(inherits(fit, "gadget.fit"))
  
  res <- fit$res.by.year
  std <- fit$stock.std
  
  ## Expand age and step to a named list with an element per stock
  stocks <- unique(std$stock)
  as_stock_list <- function(x, arg){
    if (is.null(x)) return(list())
    if (is.list(x)){
      if (is.null(names(x)) || any(names(x) == "")){
        stop("If '", arg, "' is a list, all elements must be named by stock")
      }
      if (length(setdiff(names(x), stocks)) > 0){
        warning("The following stocks in '", arg, "' are not in fit$stock.std: ",
                paste(setdiff(names(x), stocks), collapse = ", "))
      }
      return(x)
    }
    return(stats::setNames(rep(list(x), length(stocks)), stocks))
  }
  age_by_stock <- is.list(age)
  age <- as_stock_list(age, "age")
  step <- as_stock_list(step, "step")
  
  ## Numbers are not summed over ages or steps, as this would count individuals more than once
  for (arg in c("age", "step")){
    x <- get(arg)
    if (!all(vapply(x, function(y) is.numeric(y) && length(y) == 1, logical(1)))){
      stop("'", arg, "' must contain a single number for each stock")
    }
  }
  
  ## Join on all common columns (including 'id' for combined fits)
  by <- intersect(names(res), names(std))
  
  ## stock.std contains numbers at the start of the step, so the numbers at the end
  ## of a step are those at the start of the following step (step 1 of the following
  ## year for the last step)
  n_steps <- max(std$step)
  std$year <- ifelse(std$step == 1, std$year - 1, std$year)
  std$step <- ifelse(std$step == 1, n_steps, std$step - 1)
  
  for (s in intersect(names(age), stocks)){
    
    col <- paste0("rec", age[[s]])
    if (!col %in% names(res)) res[[col]] <- NA_real_
    
    rec <- std[std$stock == s & std$age == age[[s]], ]
    if (nrow(rec) == 0){
      if (age_by_stock) warning("Age ", age[[s]], " not found for stock '", s, "' in fit$stock.std")
      next
    }
    
    ## Specified step, otherwise the res.by.year step(s)
    s_by <- by
    if (!is.null(step[[s]])){
      rec <- rec[rec$step == step[[s]], ]
      s_by <- setdiff(by, "step")
    }
    
    ind <- res$stock == s
    res[[col]][ind] <-
      rec$number[match(do.call(paste, c(res[ind, s_by, drop = FALSE], sep = "\r")),
                       do.call(paste, c(rec[s_by], sep = "\r")))]
  }
  
  if (!age_by_stock && all(is.na(res[[paste0("rec", age[[1]])]]))){
    warning("Age ", age[[1]], " not found for any stock in fit$stock.std")
  }
  
  fit$res.by.year <- res
  return(fit)
}

#' Adds the average fishing mortality over an age range (Fbar) to the annual output of a gadget.fit object
#'
#' @description Calculates the annual fishing mortality averaged over a range of ages (Fbar) for one or more
#' groups of stocks, e.g. the immature and mature stocks of a species, and adds it to \code{fit$res.by.year}
#' as an additional column named \code{paste0("fbar", min(age), "_", max(age))}. The value of each group is
#' given to all stocks in the group, and is NA for stocks not included in any group. Stocks, and areas, within
#' a group are pooled. Only catches by fleets (\code{fit$stock.prey}) are included.
#'
#' In gadget3, catches are removed at the beginning of each step, before natural mortality, so the fishing
#' mortality in a step is \eqn{F = -log(1 - C/N)}, where \eqn{C} is the number caught and \eqn{N} is the number at the
#' beginning of the step. The annual F is the sum of F over the steps of the year. Two methods are available:
#' \describe{
#'   \item{"ices"}{(default) F is calculated for each age and step, with C and N summed over the stocks in the group. F at
#'     age is summed over the steps of the year, and Fbar is the unweighted mean of F at age over the age range.}
#'   \item{"gadget"}{(legacy) F is calculated for each step with C and N summed over the ages in the age range (and
#'     stocks in the group), and Fbar is the sum over the steps of the year. Fbar is therefore weighted by the numbers
#'     at age, so tends to be dominated by the most numerous (usually youngest) ages.}
#' }
#'
#' @param fit A gadget.fit object, see \code{\link{g3_fit}}
#' @param age The age range, i.e. a vector of the minimum and maximum age (or all ages, e.g. \code{3:6}). Either
#'   a vector applied to all groups, or a list with an age range for each group (named as \code{stocks}),
#'   e.g. \code{list(cod = c(4, 7), had = c(2, 5))}.
#' @param stocks The groups of stocks. Either NULL, in which case all stocks caught by fleets are pooled into a
#'   single group, a vector of stock names (a single group), or a named list of vectors of stock names,
#'   e.g. \code{list(cod = c("cod_imm", "cod_mat"), had = "had")}. A stock can only belong to one group.
#' @param method Either "ices" (default) or "gadget", see Description.
#' @return The gadget.fit object with additional column(s) in \code{fit$res.by.year}. As columns are named by age
#'   range, groups with different age ranges will have their Fbar in different columns.
#' @examples
#' \dontrun{
#' fit <- g3_fit(model, params)
#' ## All stocks pooled
#' fit <- g3_fit.update_fbar(fit, age = c(4, 7))
#' ## Two species
#' fit <- g3_fit.update_fbar(fit,
#'                           age = list(cod = c(4, 7), had = c(2, 5)),
#'                           stocks = list(cod = c("cod_imm", "cod_mat"), had = "had"))
#' }
#' @export
g3_fit.update_fbar <- function(fit, age, stocks = NULL, method = c("ices", "gadget")){
  
  stopifnot(inherits(fit, "gadget.fit"))
  method <- match.arg(method)
  
  sp <- fit$stock.prey
  if (is.null(sp) || nrow(sp) == 0) stop("There are no fleet catches in fit$stock.prey")
  
  ## Groups of stocks
  if (is.null(stocks)) stocks <- unique(sp$stock)
  if (!is.list(stocks)) stocks <- list(stocks)
  if (is.null(names(stocks))) names(stocks) <- seq_along(stocks)
  
  if (anyDuplicated(unlist(stocks))){
    stop("A stock can only belong to one group: ", 
         paste(unique(unlist(stocks)[duplicated(unlist(stocks))]), collapse = ", "))
  }
  if (length(setdiff(unlist(stocks), sp$stock)) > 0){
    warning("The following stocks have no fleet catches in fit$stock.prey: ",
            paste(setdiff(unlist(stocks), sp$stock), collapse = ", "))
  }
  
  ## Age range for each group
  if (is.list(age)){
    if (!setequal(names(age), names(stocks))){
      stop("If 'age' is a list, it must have an age range for each group in 'stocks', named as 'stocks'")
    }
  } else {
    age <- stats::setNames(rep(list(age), length(stocks)), names(stocks))
  }
  if (!all(vapply(age, function(x) is.numeric(x) && length(x) > 0, logical(1)))){
    stop("'age' must be numeric")
  }
  
  res <- fit$res.by.year
  id <- intersect("id", names(sp))   # combined fits
  
  for (g in names(stocks)){
    
    age_range <- range(age[[g]])
    col <- paste0("fbar", age_range[1], "_", age_range[2])
    if (!col %in% names(res)) res[[col]] <- NA_real_
    
    dat <- sp[sp$stock %in% stocks[[g]] & sp$age >= age_range[1] & sp$age <= age_range[2], ]
    
    missing_ages <- setdiff(seq(age_range[1], age_range[2]), dat$age)
    if (length(missing_ages) > 0){
      warning("The following ages are not in the stocks of group '", g, "': ", 
              paste(missing_ages, collapse = ", "))
    }
    if (nrow(dat) == 0) next
    
    fbar <- calc_fbar(dat, by = id, method = method)
    
    ind <- res$stock %in% stocks[[g]]
    res[[col]][ind] <- 
      fbar$F[match(do.call(paste, c(res[ind, c(id, "year"), drop = FALSE], sep = "\r")),
                   do.call(paste, c(fbar[c(id, "year")], sep = "\r")))]
  }
  
  fit$res.by.year <- res
  return(fit)
}
