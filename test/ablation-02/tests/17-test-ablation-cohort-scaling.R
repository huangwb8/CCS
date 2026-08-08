#!/usr/bin/env Rscript

# Purpose: Lock the two-dimensional cohort-evidence scaling contract.
# Input: Synthetic Direct features, unequal tissue banks, and external queries.
# Parameters: Breadth/depth repeats, matched bank sizes, and fixed readout settings.
# Output: Schema-v2 design, metrics, contrasts, diagnostics, and audit reasons.

luckyBase::Plus.library(c("xgboost", "digest"))
source(file.path("R", "ablation.R"))

stopifnot(exists(".ablation_cohort_bank_design", mode = "function"))
stopifnot(exists(".ablation_representation_scaling", mode = "function"))
stopifnot(exists(".ablation_representation_scaling_summary", mode = "function"))

set.seed(701)
reference_label <- rep(rep(c("A", "B"), each = 8), 4)
reference_signal <- ifelse(reference_label == "A", -1, 1)
reference_direct <- cbind(
  tsp1 = reference_signal + stats::rnorm(64, 0, 0.25),
  tsp2 = -reference_signal + stats::rnorm(64, 0, 0.25),
  tsp3 = stats::rnorm(64),
  tsp4 = stats::rnorm(64),
  tsp5 = stats::rnorm(64),
  tsp6 = stats::rnorm(64),
  single_bin = reference_signal + stats::rnorm(64, 0, 0.1),
  set_pair = stats::rnorm(64)
)
rownames(reference_direct) <- sprintf("RS%03d", seq_len(nrow(reference_direct)))
reference_metadata <- data.frame(
  sample_id = rownames(reference_direct),
  cohort = rep(paste0("R", seq_len(8)), each = 8),
  cancer_type = reference_label,
  tissue = rep(c("T1", "T2", "T3", "T1"), each = 16),
  platform_id = rep(c("P1", "P2"), 32),
  stringsAsFactors = FALSE
)

query_label <- rep(c("A", "B"), each = 12)
query_signal <- ifelse(query_label == "A", -1, 1)
query_direct <- cbind(
  tsp1 = query_signal + stats::rnorm(24, 0, 0.25),
  tsp2 = -query_signal + stats::rnorm(24, 0, 0.25),
  tsp3 = stats::rnorm(24),
  tsp4 = stats::rnorm(24),
  tsp5 = stats::rnorm(24),
  tsp6 = stats::rnorm(24),
  single_bin = query_signal + stats::rnorm(24, 0, 0.1),
  set_pair = stats::rnorm(24)
)
rownames(query_direct) <- sprintf("QS%03d", seq_len(nrow(query_direct)))
query_metadata <- data.frame(
  sample_id = rownames(query_direct),
  cohort = rep(c("Q1", "Q2"), each = 12),
  cancer_type = query_label,
  tissue = rep(c("T1", "T3"), each = 12),
  platform_id = rep(c("P1", "P2"), each = 12),
  d1_provenance = "external_frozen",
  stringsAsFactors = FALSE
)

module_ids <- paste0("M", seq_len(8))
module_tissues <- rep(c("T1", "T2", "T3"), c(3, 3, 2))
make_d1 <- function(signal, sample_id) {
  blocks <- lapply(seq_along(module_ids), function(i) {
    evidence <- stats::plogis(signal + stats::rnorm(length(signal), 0, 0.2))
    cbind(evidence, 1 - evidence)
  })
  result <- do.call(cbind, blocks)
  colnames(result) <- unlist(lapply(
    module_ids,
    function(module_id) paste(module_id, c("A", "B"), sep = "|")
  ))
  rownames(result) <- sample_id
  result
}
reference_d1 <- make_d1(reference_signal, rownames(reference_direct))
query_d1 <- make_d1(query_signal, rownames(query_direct))
selected_blocks <- split(seq_len(ncol(reference_d1)), rep(module_ids, each = 2))

modules <- data.frame(
  module_id = module_ids,
  tissue = module_tissues,
  cohort = paste0("C", seq_len(8)),
  stringsAsFactors = FALSE
)
prepared <- list(
  reference_direct = reference_direct,
  query_direct = query_direct,
  reference_d1 = reference_d1,
  query_d1 = query_d1,
  reference_metadata = reference_metadata,
  query_metadata = query_metadata,
  module_manifest = list(modules = modules),
  selected_blocks = selected_blocks,
  selected_module_ids = module_ids,
  feature_manifest = list(
    tsp_features = colnames(reference_direct)[seq_len(6)],
    features = colnames(reference_direct)
  ),
  cache_key = "synthetic-cohort-bank"
)
config <- list(
  anchors = list(technical = "platform_id"),
  geometry = list(
    k = c(3L, 5L),
    search = "exact",
    n_trees = 10L,
    search_k = -1L,
    distance_pairs = 100L
  ),
  scaling = list(
    module_counts = c(2L, 4L, 8L),
    sequences = 2L,
    direct_feature_type = "all",
    sensitivity_feature_type = "gene_pair",
    biology_anchors = character(),
    score_reference_samples = 40L,
    score_query_samples = 20L,
    lambda = 1,
    bootstrap = 100L
  ),
  validation = list(
    lambda = c(0.1, 1),
    inner_folds = 2L,
    nrounds = 6L,
    numCores = 1L
  )
)

# Bank design invariants are independent of the scoring engine.
design_a <- .ablation_cohort_bank_design(
  modules,
  config$scaling$module_counts,
  config$scaling$sequences,
  702L
)
design_b <- .ablation_cohort_bank_design(
  modules,
  config$scaling$module_counts,
  config$scaling$sequences,
  702L
)
stopifnot(identical(design_a$design_hash, design_b$design_hash))
stopifnot(identical(design_a$design$module_ids, design_b$design$module_ids))

breadth <- design_a$design[design_a$design$design_family == "breadth", , drop = FALSE]
stopifnot(all(breadth$min_cohort_depth == 1L))
stopifnot(all(breadth$max_cohort_depth == 1L))
stopifnot(all(breadth$tissue_count == breadth$level))

depth <- design_a$design[design_a$design$design_family == "depth", , drop = FALSE]
for (repeat_id in unique(depth$repeat_id)) {
  candidate <- depth[depth$repeat_id == repeat_id, , drop = FALSE]
  stopifnot(length(unique(candidate$tissue_count)) == 1L)
  stopifnot(length(unique(vapply(
    candidate$tissues,
    function(x) paste(x, collapse = ";"),
    character(1)
  ))) == 1L)
}

sequence_design <- design_a$design[
  design_a$design$design_family %in% c("breadth", "depth") &
    !is.na(design_a$design$parent_design_id),
  ,
  drop = FALSE
]
for (i in seq_len(nrow(sequence_design))) {
  parent <- design_a$design$module_ids[[match(
    sequence_design$parent_design_id[i],
    design_a$design$design_id
  )]]
  stopifnot(all(parent %in% sequence_design$module_ids[[i]]))
  stopifnot(length(parent) < length(sequence_design$module_ids[[i]]))
}

matched <- design_a$design[design_a$design$design_family == "matched", , drop = FALSE]
for (pair_id in unique(matched$pair_id)) {
  pair <- matched[matched$pair_id == pair_id, , drop = FALSE]
  stopifnot(nrow(pair) == 2L)
  stopifnot(length(unique(pair$module_count)) == 1L)
  stopifnot(
    pair$tissue_count[pair$design_role == "breadth_heavy"] >
      pair$tissue_count[pair$design_role == "depth_heavy"]
  )
}
stopifnot(any(
  design_a$exclusions$reason == "no_tissue_diversity_contrast_at_this_size"
))

# Wider probability blocks receive the same total squared distance weight.
equal_weight_reference <- matrix(stats::rnorm(80), nrow = 20)
equal_weight_query <- matrix(stats::rnorm(40), nrow = 10)
colnames(equal_weight_reference) <- colnames(equal_weight_query) <- paste0("V", 1:4)
balanced <- .ablation_module_balanced_transform(
  equal_weight_reference,
  equal_weight_query,
  list(narrow = 1L, wide = 2:4)
)
stopifnot(all.equal(
  sum(balanced$weights[1]^2),
  sum(balanced$weights[2:4]^2)
))

scaling <- .ablation_representation_scaling(
  prepared = prepared,
  config = config,
  label_column = "cancer_type",
  seed = 702L,
  verbose = FALSE
)

stopifnot(scaling$status == "complete")
stopifnot(identical(scaling$schema_version, 2L))
stopifnot(identical(scaling$design_hash, design_a$design_hash))
stopifnot(all(c("design", "metrics", "contrasts", "diagnostics", "reasons") %in%
  names(scaling)))
stopifnot(all(c(
  "d1_feature_count",
  "query_covered_count",
  "query_covered_fraction",
  "query_coverage_hash"
) %in% colnames(scaling$design)))
stopifnot(all(scaling$design$d1_feature_count == 2L * scaling$design$module_count))
stopifnot(all(c(
  "primary_nonredundancy",
  "primary_technical",
  "primary_stability",
  "diagnostic_mechanism",
  "diagnostic_lineage"
) %in% unique(scaling$metrics$metric_role)))
stopifnot(all(
  scaling$metrics$metric_role[
    grepl("^lineage_", scaling$metrics$metric_name)
  ] == "diagnostic_lineage"
))
stopifnot(any(scaling$reasons$status == "not_evaluated"))
stopifnot(any(
  scaling$reasons$reason == "no_independent_biology_anchor_configured"
))
stopifnot(all(c(
  "breadth_slope",
  "depth_slope",
  "matched_size",
  "marginal_gain"
) %in%
  unique(scaling$contrasts$contrast_type)))
stopifnot(any(scaling$contrasts$aggregation == "bootstrap_summary"))

direct_contracts <- scaling$diagnostics$direct_contracts
stopifnot(identical(
  direct_contracts$group,
  c("Direct-GSClassifier", "Direct-GSClassifier-TSP")
))
stopifnot(identical(direct_contracts$feature_count, c(8L, 6L)))
stopifnot(length(unique(direct_contracts$feature_hash)) == 2L)
stopifnot(identical(scaling$diagnostics$score_reference_count, 40L))
stopifnot(identical(scaling$diagnostics$score_query_count, 20L))
stopifnot(length(unique(scaling$metrics$design_id)) == nrow(scaling$design))
stopifnot(identical(
  scaling$test_sample_hash,
  digest::digest(sort(query_metadata$sample_id), algo = "md5")
))

message("17-test-ablation-cohort-scaling: all tests passed")
