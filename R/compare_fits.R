## Compares two gadget.fit objects component by component, e.g. to check that changes to g3_fit()
## have not altered the output. Data frame columns are matched by name and rows are aligned by
## sorting on all common columns, so differences in column or row order are ignored. Numeric
## values are compared with a relative tolerance, and NA in one fit but not the other is a
## difference. Columns that differ in type only (integer vs double, tibble vs data.frame, or
## ages as "age1" vs 1) are compared after conversion, and reported as notes.
##
## Returns a list containing: 'components', a named logical vector (TRUE if the component
## matches); 'differences', a data frame of the differences (and notes) found; and 'FITS_MATCH'.
compare_fits <- function(fit1, fit2, tolerance = sqrt(.Machine$double.eps), verbose = TRUE){

  diffs <- data.frame(component = character(), column = character(),
                      issue = character(), note = logical())
  add_diff <- function(component, column, issue, note = FALSE){
    diffs[nrow(diffs) + 1, ] <<- list(component, column, issue, note)
  }

  ## Converts a column so columns that only differ in type can be compared
  normalise <- function(x){
    if (is.list(x)) x <- unlist(x)
    if (is.factor(x)) x <- as.character(x)
    if (is.integer(x)) x <- as.numeric(x)
    return(x)
  }
  as_numeric_age <- function(x) suppressWarnings(as.numeric(gsub("^age", "", x)))

  ## Compares two vectors, returns NULL if they match or a description of the difference
  compare_values <- function(x, y){
    if (length(x) != length(y)) return(paste("lengths differ:", length(x), "vs", length(y)))
    if (!identical(is.na(x), is.na(y))){
      return(paste0("NA in one fit only (", sum(is.na(x) != is.na(y)), " values)"))
    }
    if (is.numeric(x) && is.numeric(y)){
      ok <- !is.na(x)
      if (isTRUE(all.equal(x[ok], y[ok], tolerance = tolerance, check.attributes = FALSE))) return(NULL)
      rel <- abs(x[ok] - y[ok]) / pmax(abs(x[ok]), .Machine$double.xmin)
      return(paste0(sum(rel > tolerance), " values differ, max relative difference ", signif(max(rel), 3)))
    }
    if (identical(as.character(x), as.character(y))) return(NULL)
    return(paste0(sum(as.character(x) != as.character(y), na.rm = TRUE), " values differ"))
  }

  ## Component names
  for (i in setdiff(names(fit1), names(fit2))) add_diff(i, NA, "component in fit1 only")
  for (i in setdiff(names(fit2), names(fit1))) add_diff(i, NA, "component in fit2 only")

  for (i in intersect(names(fit1), names(fit2))){

    x <- fit1[[i]]
    y <- fit2[[i]]

    if (is.null(x) || is.null(y)){
      if (!(is.null(x) && is.null(y))) add_diff(i, NA, paste("NULL in", if (is.null(x)) "fit1" else "fit2", "only"))
      next
    }

    ## Non data frame components, e.g. params for R models
    if (!(is.data.frame(x) && is.data.frame(y))){
      issue <- compare_values(normalise(x), normalise(y))
      if (!is.null(issue)) add_diff(i, NA, issue)
      next
    }

    if (!identical(class(x), class(y))){
      add_diff(i, NA, paste0("class differs: ", class(x)[1], " vs ", class(y)[1]), note = TRUE)
    }
    x <- as.data.frame(lapply(x, normalise), stringsAsFactors = FALSE)
    y <- as.data.frame(lapply(y, normalise), stringsAsFactors = FALSE)

    for (j in setdiff(names(x), names(y))) add_diff(i, j, "column in fit1 only")
    for (j in setdiff(names(y), names(x))) add_diff(i, j, "column in fit2 only")
    cols <- intersect(names(x), names(y))

    if (!identical(cols, intersect(names(y), names(x)))) add_diff(i, NA, "column order differs", note = TRUE)
    if (nrow(x) != nrow(y)){
      add_diff(i, NA, paste("number of rows differs:", nrow(x), "vs", nrow(y)))
      next
    }

    ## Columns differing in type only, e.g. age as "age1" vs 1
    for (j in cols){
      if (is.character(x[[j]]) != is.character(y[[j]])){
        add_diff(i, j, paste0("type differs: ", class(x[[j]])[1], " vs ", class(y[[j]])[1]), note = TRUE)
        x[[j]] <- as_numeric_age(x[[j]])
        y[[j]] <- as_numeric_age(y[[j]])
      }
    }

    ## Align rows by sorting on the common columns (character columns first)
    key <- cols[order(!vapply(x[cols], is.character, logical(1)))]
    ox <- do.call(order, unname(as.list(x[key])))
    oy <- do.call(order, unname(as.list(y[key])))
    if (!identical(ox, oy)) add_diff(i, NA, "row order differs", note = TRUE)
    x <- x[ox, , drop = FALSE]
    y <- y[oy, , drop = FALSE]

    for (j in cols){
      issue <- compare_values(x[[j]], y[[j]])
      if (!is.null(issue)) add_diff(i, j, issue)
    }
  }

  ## Summary attribute, e.g. from g3_optim()
  if (!isTRUE(all.equal(attr(fit1, "summary"), attr(fit2, "summary")))){
    add_diff("attr(summary)", NA, "summary attributes differ")
  }

  components <- union(names(fit1), names(fit2))
  components <- stats::setNames(!components %in% diffs$component[!diffs$note], components)
  fits_match <- all(components) && !any(diffs$component[!diffs$note] == "attr(summary)")

  if (verbose){
    for (i in names(components)){
      cat(sprintf(" %-20s %s\n", i, if (components[[i]]) "match" else "DO NOT MATCH"))
      d <- diffs[diffs$component == i, ]
      for (k in seq_len(nrow(d))){
        cat(sprintf("     %s%s%s\n", if (d$note[k]) "(note) " else "- ",
                    if (is.na(d$column[k])) "" else paste0(d$column[k], ": "), d$issue[k]))
      }
    }
    cat(if (fits_match) "\nFITS MATCH\n" else "\nFITS DO NOT MATCH\n")
  }

  return(invisible(list(components = components, differences = diffs, FITS_MATCH = fits_match)))
}
