# TDD RED test: cohort-level paired inference for biological-anchor utility.
# The helper is intentionally referenced before implementation in RED phase.

helper_path <- file.path("test", "ablation-03", "03-ablation-biology_functions.R")
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
cat("biology anchor inference test passed\n")
