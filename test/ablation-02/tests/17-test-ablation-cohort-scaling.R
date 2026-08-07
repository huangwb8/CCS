#!/usr/bin/env Rscript

# Purpose: Lock the representation-specific cohort-bank scaling contract.
# Input: Synthetic TSP features, four frozen-module blocks, and external queries.
# Parameters: Two nested bank sizes and two tissue-balanced module sequences.
# Output: Paired pointwise contrasts plus sequence-level trend summaries.

luckyBase::Plus.library(c("xgboost", "digest"))
source(file.path("R", "ablation.R"))

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
  tsp6 = stats::rnorm(64)
)
rownames(reference_direct) <- sprintf("RS%03d", seq_len(nrow(reference_direct)))
reference_metadata <- data.frame(
  sample_id = rownames(reference_direct),
  cohort = rep(paste0("R", seq_len(8)), each = 8),
  cancer_type = reference_label,
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
  tsp6 = stats::rnorm(24)
)
rownames(query_direct) <- sprintf("QS%03d", seq_len(nrow(query_direct)))
query_metadata <- data.frame(
  sample_id = rownames(query_direct),
  cohort = rep(c("Q1", "Q2"), each = 12),
  cancer_type = query_label,
  d1_provenance = "external_frozen",
  stringsAsFactors = FALSE
)

module_ids <- c("T1|C1", "T1|C2", "T2|C3", "T2|C4")
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

prepared <- list(
  reference_direct = reference_direct,
  query_direct = query_direct,
  reference_d1 = reference_d1,
  query_d1 = query_d1,
  reference_metadata = reference_metadata,
  query_metadata = query_metadata,
  module_manifest = list(
    modules = data.frame(
      module_id = module_ids,
      tissue = c("T1", "T1", "T2", "T2"),
      cohort = paste0("C", seq_len(4)),
      stringsAsFactors = FALSE
    )
  ),
  selected_blocks = selected_blocks,
  selected_module_ids = module_ids,
  feature_manifest = list(
    tsp_features = colnames(reference_direct),
    features = colnames(reference_direct)
  )
)
config <- list(
  scaling = list(
    module_counts = c(2L, 4L),
    sequences = 2L,
    direct_feature_type = "gene_pair",
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

scaling <- .ablation_representation_scaling(
  prepared = prepared,
  config = config,
  label_column = "cancer_type",
  seed = 702L,
  verbose = FALSE
)

stopifnot(scaling$status == "complete")
stopifnot(identical(scaling$direct_group, "Direct-GSClassifier-TSP"))
stopifnot(identical(scaling$module_counts, c(2L, 4L)))
stopifnot(nrow(scaling$metrics) == 4L)
stopifnot(nrow(scaling$pointwise) == 12L)
stopifnot(nrow(scaling$trend) == 6L)
stopifnot(all(scaling$metrics$direct_feature_count == 6L))
stopifnot(all(scaling$metrics$d1_feature_count == 2L * scaling$metrics$module_count))
stopifnot(length(unique(scaling$metrics$test_sample_hash)) == 1L)

direct_by_sequence <- split(
  scaling$metrics$balanced_accuracy_direct,
  scaling$metrics$sequence_id
)
stopifnot(all(vapply(direct_by_sequence, function(x) length(unique(x)) == 1L, logical(1))))
stopifnot(all(is.finite(scaling$trend$estimate)))

message("17-test-ablation-cohort-scaling: all tests passed")
