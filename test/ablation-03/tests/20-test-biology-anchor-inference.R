# TDD RED test: cohort-level paired inference for biological-anchor utility.
# The helper is intentionally referenced before implementation in RED phase.

helper_path <- file.path("test", "ablation-03", "02.02.00. 生物锚点分析_functions.R")
if (!file.exists(helper_path)) {
  stop("Expected biology inference helper is missing.", call. = FALSE)
}
source(helper_path)

synthetic <- data.frame(
  anchor = rep("proliferation", 8L),
  representation = rep(c("Direct-GSClassifier", "Cohort-d1"), 4L),
  query_sample = rep(paste0("q", 1:4), each = 2L),
  query_cohort = rep(paste0("cohort", 1:4), each = 2L),
  utility = c(0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67)
)

result <- .biology_paired_contrast(
  synthetic,
  n_boot = 200L,
  seed = 20260903L
)

stopifnot(nrow(result) == 1L)
stopifnot(result$anchor == "proliferation")
stopifnot(result$cohort_count == 4L)
stopifnot(result$query_count == 4L)
stopifnot(is.finite(result$estimate))
stopifnot(result$ci_low <= result$estimate)
stopifnot(result$estimate <= result$ci_high)
stopifnot(is.finite(result$p_value), result$p_value >= 0, result$p_value <= 1)
stopifnot(is.finite(result$p_value_adj), result$p_value_adj >= 0, result$p_value_adj <= 1)
stopifnot(result$estimate > 0)
stopifnot(result$p_method == "exact_paired_sign_flip")

heterogeneous <- data.frame(
  anchor = rep("immune_tme", 8L),
  representation = rep(c("Direct-GSClassifier", "Cohort-d1"), 4L),
  query_sample = rep(paste0("q", 1:4), each = 2L),
  query_cohort = rep(c("cohort1", "cohort1", "cohort2", "cohort3"), each = 2L),
  utility = c(0.50, 0.70, 0.60, 0.80, 0.50, 0.40, 0.50, 0.60)
)
cohort_deltas <- .biology_cohort_deltas(heterogeneous)
heterogeneous_inference <- .biology_paired_contrast(
  heterogeneous, n_boot = 200L, seed = 20260903L
)
stopifnot(
  nrow(cohort_deltas) == 3L,
  cohort_deltas$query_count[cohort_deltas$query_cohort == "cohort1"] == 2L,
  isTRUE(all.equal(
    cohort_deltas$delta_d1_minus_direct,
    cohort_deltas$utility_d1 - cohort_deltas$utility_direct
  )),
  isTRUE(all.equal(
    mean(cohort_deltas$delta_d1_minus_direct),
    heterogeneous_inference$estimate
  ))
)
with_missing <- synthetic
with_missing$utility[with_missing$query_sample == "q4" &
  with_missing$representation == "Cohort-d1"] <- NA_real_
missing_inference <- .biology_paired_contrast(
  with_missing, n_boot = 200L, seed = 20260903L
)
stopifnot(missing_inference$cohort_count == 3L,
  missing_inference$query_count == 3L)

low_information <- .biology_paired_contrast(
  synthetic[synthetic$query_cohort %in% c("cohort1", "cohort2"), ],
  n_boot = 200L,
  seed = 20260903L,
  min_cohorts = 3L
)
stopifnot(
  low_information$inference_status == "not_estimable",
  is.finite(low_information$estimate),
  is.na(low_information$ci_low),
  is.na(low_information$ci_high),
  is.na(low_information$p_value),
  is.na(low_information$p_value_adj)
)
cat("biology anchor inference test passed\n")
