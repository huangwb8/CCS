#!/usr/bin/env Rscript

# Purpose: Prove that optimized frozen-bank encoding is identical on real CCS data.
# Input: 100 non-duplicate reference samples and the first three frozen modules.
# Parameters: Sequential encoding to keep this contract test deterministic and light.
# Output: Exact feature-contract and probability-parity assertions.

source(file.path("test", "ablation-02", "ablation-test-data.R"))
source(file.path("R", "ablation.R"))

module_manifest <- .ablation_module_manifest(resCCS)
feature_manifest <- .ablation_frozen_feature_manifest(resCCS, module_manifest)
stopifnot(length(feature_manifest$features) == 529L)
stopifnot(length(feature_manifest$tsp_features) == 496L)
stopifnot(identical(
  as.integer(table(feature_manifest$feature_manifest$feature_type)[
    c("single_bin", "gene_pair", "set_pair")
  ]),
  c(32L, 496L, 1L)
))

flattened <- .ablation_flatten_expression(reference_data)
d1 <- as.matrix(resCCS@Data$Probability$d1)
sample_ids <- intersect(rownames(d1), colnames(flattened$expr))
sample_ids <- sample_ids[seq_len(100L)]
direct <- .ablation_gsclassifier_matrix(
  resCCS,
  flattened$expr[, sample_ids, drop = FALSE],
  feature_manifest
)
module_ids <- module_manifest$modules$module_id[seq_len(3L)]

encoded <- .ablation_encode_d1_from_direct(
  object = resCCS,
  direct = direct,
  module_manifest = module_manifest,
  module_ids = module_ids,
  numCores = 1L,
  verbose = FALSE
)
expected_columns <- unlist(module_manifest$blocks[module_ids], use.names = FALSE)
expected <- d1[sample_ids, expected_columns, drop = FALSE]

stopifnot(identical(dim(encoded), dim(expected)))
stopifnot(identical(rownames(encoded), rownames(expected)))
stopifnot(identical(colnames(encoded), colnames(expected)))
stopifnot(max(abs(encoded - expected)) == 0)

message("test-ablation-real-encoder: all tests passed")
