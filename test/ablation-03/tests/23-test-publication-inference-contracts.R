# Regression tests for publication-facing inferential contracts.

source(file.path("R", "ablation.R"))
source(file.path(
  "test", "ablation-03", "02-ablation03-representation_functions.R"
))

# Exact paired sign-flip tests use the exact tail proportion, without the
# Monte Carlo +1 correction. Four same-direction clusters have p = 2 / 16.
exact_data <- data.frame(
  cohort = paste0("c", 1:4),
  delta = rep(1, 4)
)
exact_result <- .ae_paired_inference(
  exact_data,
  delta_column = "delta",
  n_boot = 200L,
  seed = 20260909L
)
stopifnot(exact_result$p_method == "exact_paired_sign_flip")
stopifnot(exact_result$p_resamples == 16L)
stopifnot(isTRUE(all.equal(exact_result$p_value, 0.125)))

# Non-estimable returns retain the same reporting columns as complete results.
missing_result <- .ae_paired_inference(
  data.frame(other = 1),
  delta_column = "delta"
)
stopifnot(
  missing_result$p_method == "not_estimable",
  missing_result$p_resamples == 0L
)

# Cohort-level readout must use accuracy, not a balanced-accuracy field whose
# equality is only accidental when each cohort contains one class.
readout <- list(paired_by_cohort = data.frame(
  cohort = c("c1", "c2"),
  delta_accuracy = c(0.2, 0.4),
  delta_balanced_accuracy = c(-0.8, -0.6),
  sample_count_d1 = c(10L, 10L),
  sample_count_direct = c(10L, 10L)
))
readout_result <- .ae_readout_inference(
  readout,
  n_boot = 200L,
  seed = 20260909L,
  multiplicity_method = "none"
)
stopifnot(readout_result$endpoint == "cohort_accuracy")
stopifnot(isTRUE(all.equal(readout_result$estimate, 0.3)))
stopifnot(readout_result$unit == "query_cohort")

# Repeated fits of the same 100% training-cohort set are not independent bank
# designs and therefore cannot produce a design-level CI or p-value.
learning <- list(paired = data.frame(
  requested_fraction = rep(1, 10L),
  repeat_id = seq_len(10L),
  cohort_subset_hash_direct = rep("all-cohorts", 10L),
  cohort_subset_hash_d1 = rep("all-cohorts", 10L),
  delta_balanced_accuracy = seq(-0.10, -0.08, length.out = 10L)
))
learning_result <- .ae_learning_inference(
  learning,
  n_boot = 200L,
  seed = 20260909L
)
stopifnot(learning_result$status == "not_estimable")
stopifnot(learning_result$n_unique_design == 1L)
stopifnot(is.na(learning_result$ci_low), is.na(learning_result$p_value))

# Scaling slopes are explicitly per doubling of module count.
design <- data.frame(
  design_id = paste0("d", 1:6),
  parent_design_id = NA_character_,
  mean_cohort_depth = 1,
  stringsAsFactors = FALSE
)
metrics <- data.frame(
  design_id = design$design_id,
  design_family = "breadth",
  design_role = "sequence",
  repeat_id = rep(c("r1", "r2"), each = 3L),
  level = rep(1:3, 2L),
  pair_id = NA_character_,
  module_count = rep(c(1L, 2L, 4L), 2L),
  tissue_count = rep(1:3, 2L),
  metric_name = "normalized_effective_rank",
  metric_role = "primary_nonredundancy",
  query_coverage = "all",
  estimate = rep(c(0, 1, 2), 2L),
  status = "evaluated",
  reason = NA_character_,
  stringsAsFactors = FALSE
)
scaling_result <- .ablation_representation_scaling_summary(
  metrics,
  design,
  bootstrap = 200L,
  seed = 20260909L
)
scaling_summary <- scaling_result[
  scaling_result$aggregation == "bootstrap_summary",
  ,
  drop = FALSE
]
stopifnot(scaling_summary$component == "per_module_count_doubling")
stopifnot(isTRUE(all.equal(scaling_summary$estimate, 1)))

cat("publication inference contract tests passed\n")
