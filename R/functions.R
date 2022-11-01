## -----------------------------------------------------------------------------
##
## Miscellaneous functions used in gadgetutils
##
## -----------------------------------------------------------------------------

#' @export
echo_message <- function(...){
  system(sprintf("echo '%s'", paste0(..., collapse = '')))
}

save_obj <- function(..., file) {
  x <- list(...)
  save(list = names(x), file = file, envir = list2env(x))
}