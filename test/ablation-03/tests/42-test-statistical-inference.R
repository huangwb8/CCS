# Small synthetic cases for the new conditional inference helpers.
source(file.path("test", "ablation-03", "targets", "statistical_inference.R"))

local_data <- data.frame(
  cohort = rep(c("a", "b"), each = 20L),
  delta = c(rep(1, 20L), rep(-1, 20L))
)
local_result <- .asi_local(local_data, "cohort", "delta", "paired_accuracy",
  min_n = 20L, n_boot = 100L, seed = 17L)
stopifnot(nrow(local_result) == 2L,
  identical(local_result$unit, rep("query_within_cohort", 2L)),
  all(local_result$ci_low == local_result$estimate),
  all(local_result$ci_high == local_result$estimate))

sparse <- data.frame(
  cohort_a = c("a", "a", "b"), cohort_b = c("b", "c", "c"),
  same_cancer_type = TRUE, delta_d1_minus_direct = c(0.1, -0.2, 0.3),
  direction = "reference_bank_to_external_targets"
)
structural <- .asi_structural(sparse, n_boot = 100L, seed = 19L)
stopifnot(nrow(structural) == 2L, all(is.na(structural$p_value)),
  all(structural$status == "not_estimable"),
  all(structural$n_pairs == 3L))

truth <- c(0, 0, 0)
prediction <- c(0.1, 0.2, 0.3)
metric <- .asi_decoder_metric(truth, prediction, "gene_pair")
stopifnot(is.na(metric[["balanced_accuracy"]]),
  is.finite(metric[["brier"]]))
stopifnot(is.na(.asi_decoder_metric(1, 1.1, "single_bin")[["spearman"]]))

cat("statistical inference helper tests passed\n")
