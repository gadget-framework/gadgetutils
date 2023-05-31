#' @title Add gadget3 attributes to an existing data frame
#' @description Adds gadget3 attributes and missing columns to a a data frame ready to be passed to gadget3. Useful for importing data manually to gadget without the mfdb step.
#' @param x A data frame, tibble or data table.
#' @param params A list of parameters defining the grouping. See details in \code{mfdb_sample_count}.
#' @param warn Integer value to set the handling of warning messages. See the \code{warn} argument in \link[base]{options} for details. In brief, set to negative to ignore warnings, to 1 to print warnings as they occur and to 2 to stop executing and return an error instead of a warning. The last is meant to prevent passing information that may cause unintended behavior in gadget3 models and is recommended. 
#' @return Returns \code{x} with g3 attributes and missing columns as instructed by the \code{params} argument.
#' @export
# @examples
# read.csv('../ghl-gadget-data/data/out/Winter survey index.csv') %>%
#   dplyr::select(startyear, bm) %>%
#   rename("year" = "startyear", "weight" = "bm") %>%
#   add_g3_attributes(
#     params = list(
#       year = model_params$year_range,
#       step = model_params$timestep_fun,
#       length = mfdb_interval(
#         "all", c(40, stock_params$maxlength),
#         open_ended = c("upper"))
#     )
#   )

add_g3_attributes <- function(x, params, warn = 1) {
  
  ## Set warn
  
  orig_warn <- options()$warn
  options(warn = warn)
  
  ### Correct inconsistency in naming timestep
  
  if("timestep" %in% names(params) & !"step" %in% names(params)) {
    names(params)[names(params) == "timestep"] <- "step"
  }
  
  ## Column/attribute order
  
  param_order <- names(params)
  attrib_order <- c(names(attributes(x)), param_order)
  
  ## Assign x to new object for easier debugging
  
  out <- x
  
  ### Add attributes
  
  if("year" %in% names(params)) {
    attributes(out)$year <-
      stats::setNames(
        lapply(unique(params$year), function(k) k),
        lapply(unique(params$year), function(k) k)
      )
    params <- params[!names(params) %in% "year"] # Remove year from params to avoid duplicates
  } else {
    attributes(out)$year <-
      stats::setNames(
        lapply(unique(out$year), function(k) k),
        lapply(unique(out$year), function(k) k)
      )
  }
  
  if(any(sapply(params, function(k) inherits(k, "mfdb_group")))) {
    tmp_params <- sapply(params, function(k) inherits(k, "mfdb_group"))
    tmp_params <- names(tmp_params[tmp_params])
    attributes(out) <-
      c(attributes(out),
        stats::setNames(
          lapply(tmp_params, function(k) {
            params[[k]]
          }),
          tmp_params
        )
      )
  }
  
  if(any(sapply(params, function(k) inherits(k, "list")))) {
    tmp_params <- sapply(params, function(k) inherits(k, "list"))
    tmp_params <- names(tmp_params[tmp_params])
    attributes(out) <-
      c(attributes(out),
        stats::setNames(
          lapply(tmp_params, function(k) {
            params[[k]]
          }),
          tmp_params
        )
      )
  }
  
  if(any(sapply(params, function(k) inherits(k, "mfdb_interval")))) {
    tmp_params <- sapply(params, function(k) inherits(k, "mfdb_interval"))
    tmp_params <- names(tmp_params[tmp_params])
    attributes(out) <-
      c(attributes(out),
        stats::setNames(
          lapply(tmp_params, function(k) {
            tmp <- params[[k]]
            
            stats::setNames(lapply(seq_along(tmp)[-length(tmp)], function(i) {
              
              min_val <- unname(tmp[i])
              max_val <- unname(tmp[i+1])
              outout <- call("seq", min_val, max_val -1)
              
              attr(outout, "min") <- min_val
              attr(outout, "max") <- max_val
              
              if("lower" %in% attributes(tmp)$open_ended & i == 1) {
                attr(outout, "min_open_ended") <- TRUE
              }
              if("upper" %in% attributes(tmp)$open_ended & i == length(tmp) -1) {
                attr(outout, "max_open_ended") <- TRUE
              }
              
              outout
            }),
            names(tmp)[-length(tmp)]
            )
          }),
          tmp_params
        )
      )
  }
  
  attributes(out)$generator <- "add_g3_attributes"
  
  attributes(out) <- attributes(out)[c(attrib_order, setdiff(names(attributes(out)), attrib_order))]
  
  ### Add missing columns
  
  if(length(setdiff(names(params), colnames(out))) > 0) {
    
    outout <- cbind(out,
                    stats::setNames(lapply(setdiff(names(params), colnames(out)), function(k) {
                      rep(names(attributes(out)[names(attributes(out)) == k][[1]]), nrow(out))
                    }), setdiff(names(params), colnames(out)))
    )
    
    outout <- outout[c(intersect(names(outout), names(attributes(out))),
                       setdiff(names(outout), names(attributes(out))))]
    
    attributes(outout) <- c(attributes(outout),
                            attributes(out)[setdiff(names(attributes(out)), names(attributes(outout)))])
    
    out <- outout
  }
  
  ## Test whether attribute names match those of data
  
  test_attributes <- names(attributes(out)[names(attributes(out)) %in% names(out)])
  
  lapply(test_attributes, function(k) {
    if(!all(unique(out[[k]]) %in% names(attributes(out)[names(attributes(out)) == k][[1]]))) {
      warning(
        paste(
          "Unique column and attribute values for", k, "do not match:", 
          paste(setdiff(unique(out[[k]]), names(attributes(out)[names(attributes(out)) == k][[1]])), collapse = ", ")
        )
      )
    }
  })
  
  ## Set warn back
  
  options(warn = orig_warn)
  
  ## Return
  
  return(out)
  
}
