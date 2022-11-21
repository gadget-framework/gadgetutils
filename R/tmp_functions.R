g3_ageareadata <- function(lookup_name, df, value_field = 'value', areas = NULL) {
  stopifnot(is.null(areas) || !anyNA(as.integer(areas)))
  
  age_convert <- function(age, step){
    as.integer(age) * 100L + as.integer(step)
  }
  
  # What's the next power of 10?
  next_mult <- function (x) as.integer(10 ** ceiling(log10(x)))
  
  remove_multzero <- function (f) gadget3:::call_replace(f, "*" = function (x) {
    if (isTRUE(all.equal(x[[2]], 0))) return(0)
    if (isTRUE(all.equal(x[[3]], 0))) return(0)
    return(as.call(list(
      x[[1]],
      remove_multzero(x[[2]]),
      remove_multzero(x[[3]]))))
  })
  
  for (n in c('age', value_field)) {
    if (is.null(df[[n]])) stop("No ", n, " field in g3_timeareadata data.frame")
  }
  
  # All lengths the same, no point adding to the lookup
  times <- age_convert(df$age, df$step)  # NB: if step column missing, this will be NULL
  age_mult <- 100L
  
  # If have a non-identity area map, apply it to data first
  if (!(is.null(areas) || identical(names(areas), as.character(areas)))) {
    storage.mode(areas) <- "integer"  # Integer-ize without losing names
    df$area <- as.integer(areas[df$area])
  }
  
  # Count potential areas, 0, 1, many
  area_count <- if (is.null(df$area)) 0 else if (length(df$area) > 1 && any(df$area[[1]] != df$area)) 2 else 1
  area_mult <- if (area_count > 1) next_mult(max(times)) else 0L
  
  lookup <- gadget3:::g3_intlookup(lookup_name,
                                   keys = as.integer(times + (if (area_count > 1) area_mult * df$area else 0)),
                                   values = df[[value_field]])
  
  # Return formula that does the lookup
  out_f <- lookup('getdefault', remove_multzero(gadget3:::f_substitute(
    quote( area * area_mult + age * age_mult + cur_step * step_mult ),
    list(
      area_mult = area_mult,
      # Mult is zero ==> There is no step.
      step_mult = if (age_mult > 1L) 1L else 0L,
      age_mult = age_mult))), 0)
  
  if (area_count == 1) {
    # Wrap lookup with check that we're in the correct area
    out_f <- gadget3:::f_substitute(
      quote( if (area != our_area) 0 else out_f ),
      list(our_area = df$area[[1]], out_f = out_f))
  }
  
  return(out_f)
}
