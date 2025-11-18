library(unittest)
library(gadgetutils)

library(unittest)

# Convert a string into a data.frame
table_string <- function (text, ...) {
  out <- read.table(
    text = text,
    blank.lines.skip = TRUE,
    header = TRUE,
    stringsAsFactors = FALSE,
    ...)
  rownames(out) <- out$switch
  return(out)
}

ok(ut_cmp_identical(g3_iterative_default_grouping(table_string('
switch					value
cdist_sumofsquares_comm_ldist_weight	1
cdist_sumofsquares_comm_aldist_weight	1
cdist_sumofsquares_comm_argle_weight	1
cdist_sumofsquares_comm_matp_weight	1
cdist_sumofsquares_fgn_ldist_weight	1
cdist_sumofsquares_fgn_aldist_weight	1
cdist_surveyindices_log_surv_si_weight	1
'), nll_dist_names = c("ldist", "aldist", "matp", "si")), list(
  # NB: argle is missing
  comm = c("comm_ldist", "comm_aldist", "comm_matp"),
  fgn = c("fgn_ldist", "fgn_aldist"),
  # NB: parsed the awkward surveyindices_log
  surv = c("surv_si")
)))

ok(ut_cmp_identical(suppressWarnings(g3_iterative_default_grouping(table_string('
switch					value
cdist_sumofsquares_comm_ldist_weight	1
cdist_sumofsquares_comm_aldist_weight	1
cdist_sumofsquares_comm_argle_weight	1
cdist_sumofsquares_comm_matp_weight	1
cdist_sumofsquares_fgn_ldist_weight	1
cdist_sumofsquares_fgn_aldist_weight	0
cdist_surveyindices_log_surv_si_weight	1
'), nll_dist_names = c("ldist", "aldist", "matp", "si"))), list(
  comm = c("comm_ldist", "comm_aldist", "comm_matp"),
  # NB: zero-weighted doesn't count
  fgn = c("fgn_ldist"),
  surv = c("surv_si")
)))
