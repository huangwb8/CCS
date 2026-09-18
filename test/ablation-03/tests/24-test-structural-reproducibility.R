# Regression tests for structural-reproducibility contracts.

source(file.path(
  "test", "ablation-03",
  "07.00.00. 结构复现分析_functions.R"
))
source(file.path("R", "ablation.R"))

# Historical Undefined module labels must resolve before bank partitioning.
module_rows <- rbind(
  data.frame(
    module_id = paste0("A|r", seq_len(10L)),
    tissue = "A", cohort = paste0("r", seq_len(10L))
  ),
  data.frame(
    module_id = paste0("B|r", 11:22),
    tissue = "B", cohort = paste0("r", 11:22)
  ),
  data.frame(
    module_id = paste0("C|r", 23:150),
    tissue = "C", cohort = paste0("r", 23:150)
  ),
  data.frame(
    module_id = c("Undefined|e1", paste0("A|e", 2:10)),
    tissue = c("Undefined", rep("A", 9L)),
    cohort = paste0("e", seq_len(10L))
  ),
  data.frame(
    module_id = paste0("B|e", 11:22),
    tissue = "B", cohort = paste0("e", 11:22)
  ),
  data.frame(
    module_id = paste0("D|e", 23:43),
    tissue = "D", cohort = paste0("e", 23:43)
  )
)
module_rows$block_width <- 1L
module_rows$first_column <- seq_len(nrow(module_rows))
module_rows$last_column <- module_rows$first_column
resolution_audit <- data.frame(
  cohort = "e1",
  resolved_tissue = "A",
  stringsAsFactors = FALSE
)
external_cohorts <- c(
  paste0("A/e", seq_len(10L)),
  paste0("B/e", 11:22),
  paste0("D/e", 23:43)
)
module_table <- .asr_resolve_module_table(
  list(modules = module_rows),
  resolution_audit,
  external_cohorts
)
stopifnot(
  !any(module_table$tissue == "Undefined"),
  sum(module_table$bank_role == "reference") == 150L,
  sum(module_table$bank_role == "external") == 43L,
  length(intersect(
    module_table$module_id[module_table$bank_role == "reference"],
    module_table$module_id[module_table$bank_role == "external"]
  )) == 0L
)
matched_design <- .asr_matched_bank_design(
  module_table,
  n_repeats = 3L,
  seed = 20260912L
)
matched_counts <- table(matched_design$repeat_id, matched_design$bank_role)
stopifnot(all(matched_counts == 22L))
matched_tissue_counts <- table(
  matched_design$repeat_id,
  matched_design$tissue,
  matched_design$bank_role
)
stopifnot(all(
  matched_tissue_counts[, , "reference"] ==
    matched_tissue_counts[, , "external"]
))

set.seed(20260912L)
cohorts <- c("A/c1", "A/c2", "B/c3")
metadata <- data.frame(
  sample_id = paste0("s", seq_len(90L)),
  cohort_key = rep(cohorts, each = 30L),
  cancer_type = rep(c("A", "A", "B"), each = 30L),
  assay_type = "RNA",
  source_system = "fixture",
  stringsAsFactors = FALSE
)
scores <- do.call(rbind, lapply(c("anchor_1", "anchor_2"), function(anchor) {
  data.frame(
    sample_id = metadata$sample_id,
    cohort_key = metadata$cohort_key,
    anchor = anchor,
    score = rep(seq_len(30L), 3L) +
      ifelse(anchor == "anchor_2", rep(c(0, 2, 4), each = 30L), 0),
    gene_count = 10L,
    stringsAsFactors = FALSE
  )
}))
states <- .asr_assign_anchor_states(
  scores,
  tail_fraction = 0.25,
  min_entity_n = 5L
)
stopifnot(
  length(unique(states$entity)) == 4L,
  all(table(states$cohort_key, states$entity) == 7L)
)

base <- cbind(
  x = rep(seq_len(30L), 3L),
  y = rep(sin(seq_len(30L) / 5), 3L)
)
rownames(base) <- metadata$sample_id
direct <- .asr_build_cohort_geometries(base, metadata, states)
d1 <- .asr_build_cohort_geometries(base * 2, metadata, states)
pair_data <- .asr_compare_cohort_pairs(
  direct$geometries,
  d1$geometries,
  unique(metadata[, c("cohort_key", "cancer_type")]),
  min_shared_entities = 4L
)
stopifnot(
  nrow(pair_data) == 3L,
  all(pair_data$shared_entity_count == 4L),
  all(abs(pair_data$delta_d1_minus_direct) < 1e-12),
  sum(pair_data$same_cancer_type) == 1L
)

summary <- .asr_summarize_pairs(pair_data, n_boot = 200L, seed = 20260912L)
stopifnot(
  setequal(summary$scope, c("all_cohort_pairs", "same_cancer_type")),
  all(abs(summary$mean_delta) < 1e-12),
  all(is.na(summary$p_value)),
  all(summary$inference_status == "not_estimable"),
  all(is.na(summary$mean_delta_ci_low)),
  all(is.na(summary$mean_delta_ci_high))
)

summary_no_boot <- .asr_summarize_pairs(
  pair_data,
  n_boot = 0L,
  seed = 20260912L
)
stopifnot(
  all(is.na(summary_no_boot$mean_delta_ci_low)),
  all(is.na(summary_no_boot$mean_delta_ci_high)),
  all(summary_no_boot$valid_bootstrap == 0L)
)

# Direction metadata must survive the shared evaluator unchanged.
fit_ids <- paste0("f", seq_len(30L))
fit_d1 <- matrix(rnorm(120L), nrow = 30L, dimnames = list(fit_ids, NULL))
target_d1 <- cbind(base, base^2)
colnames(fit_d1) <- colnames(target_d1) <- paste0("d", seq_len(4L))
direction_base <- list(
  fit_sample_ids = fit_ids,
  target_sample_ids = rownames(target_d1),
  target_metadata = metadata,
  states = states,
  direct_geometry = direct
)
direction_result <- .asr_evaluate_direction(
  direction_base,
  fit_d1,
  target_d1,
  blocks = list(module_1 = 1:2, module_2 = 3:4),
  direction = "external_bank_to_reference_targets",
  bank_role = "external",
  target_role = "reference",
  min_shared_entities = 4L,
  n_boot = 0L,
  seed = 20260912L
)
stopifnot(
  all(direction_result$summary$direction ==
    "external_bank_to_reference_targets"),
  all(direction_result$summary$bank_role == "external"),
  all(direction_result$summary$target_role == "reference"),
  all(direction_result$pair_comparisons$direction ==
    "external_bank_to_reference_targets")
)

cat("structural reproducibility contract tests passed\n")
