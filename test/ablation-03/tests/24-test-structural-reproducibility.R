# Regression tests for structural-reproducibility contracts.

source(file.path(
  "test", "ablation-03",
  "04-ablation03-structural-reproducibility_functions.R"
))

set.seed(20260912L)
cohorts <- c("A/c1", "A/c2", "B/c3")
metadata <- data.frame(
  sample_id = paste0("s", seq_len(90L)),
  cohort_key = rep(cohorts, each = 30L),
  cancer_type = rep(c("A", "A", "B"), each = 30L),
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
  all(is.na(summary$p_value))
)

cat("structural reproducibility contract tests passed\n")
