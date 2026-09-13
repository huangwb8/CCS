# Endpoint-specific external candidate qualification contracts.

source(file.path("R", "ablation.R"))

reference <- data.frame(
  sample_id = paste0("r", seq_len(6L)),
  cohort = rep(c("r1", "r2", "r3"), each = 2L),
  cohort_key = rep(c("T/r1", "T/r2", "U/r3"), each = 2L),
  cancer_type = rep(c("A", "A", "B"), each = 2L),
  stringsAsFactors = FALSE
)
query <- data.frame(
  sample_id = paste0("q", seq_len(4L)),
  cohort = paste0("e", seq_len(4L)),
  cohort_key = paste0("T/e", seq_len(4L)),
  cancer_type = c("A", "B", "C", "A"),
  stringsAsFactors = FALSE
)
config <- list(anchors = list(min_reference_cohorts = 2L))
audit <- .ablation_endpoint_eligibility(reference, query, config)

stopifnot(
  nrow(audit) == 28L,
  all(audit$candidate_status == "candidate"),
  all(audit$qualification_status[audit$endpoint == "geometry"] == "estimable"),
  sum(audit$endpoint == "cancer_readout" &
    audit$qualification_status == "estimable") == 2L,
  all(audit$qualification_status[
    audit$cancer_type == "C" & audit$endpoint == "cancer_readout"
  ] == "not_estimable"),
  all(audit$qualification_reason[
    audit$cancer_type == "C" & audit$endpoint == "cancer_readout"
  ] == "no_reference_cancer_support"),
  all(audit$fit_target_overlap == 0L),
  all(nzchar(audit$sample_hash)),
  all(nzchar(audit$config_hash))
)

# Endpoint thresholds are independently configurable without changing the
# shared candidate set.
config$anchors$endpoint_min_reference_cohorts <- list(cancer_readout = 3L)
stricter <- .ablation_endpoint_eligibility(reference, query, config)
stopifnot(
  all(stricter$qualification_status[stricter$endpoint == "geometry"] == "estimable"),
  all(stricter$qualification_status[stricter$endpoint == "cancer_readout"] ==
    "not_estimable")
)

cat("endpoint eligibility contracts passed\n")
