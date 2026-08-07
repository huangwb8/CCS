#!/usr/bin/env Rscript

# Purpose: Lock feature-type-specific information-loss metrics for d1 decoding.
# Input: Low-rank synthetic d1 and reconstructable Direct-GSClassifier features.
# Parameters: Rank two ridge decoder with a fixed penalty.
# Output: Per-feature and feature-type summaries with the planned metric families.

luckyBase::Plus.library(c("irlba", "digest"))
source(file.path("R", "ablation.R"))

stopifnot(exists(".ablation_decode_direct_features", mode = "function"))

set.seed(801)
reference_d1 <- matrix(stats::rnorm(400), nrow = 200, ncol = 2)
query_d1 <- matrix(stats::rnorm(160), nrow = 80, ncol = 2)
colnames(reference_d1) <- colnames(query_d1) <- c("M1|1", "M1|2")
rownames(reference_d1) <- sprintf("R%03d", seq_len(nrow(reference_d1)))
rownames(query_d1) <- sprintf("Q%03d", seq_len(nrow(query_d1)))

make_direct <- function(z) {
  cbind(
    GENE1 = as.integer(cut(z[, 1], breaks = c(-Inf, -0.5, 0.5, Inf))) - 1L,
    `GENE1:GENE2` = as.integer(z[, 1] > z[, 2]),
    s1s2 = z[, 1] - z[, 2]
  )
}
reference_direct <- make_direct(reference_d1)
query_direct <- make_direct(query_d1)
feature_manifest <- data.frame(
  feature = colnames(reference_direct),
  feature_type = c("single_bin", "gene_pair", "set_pair"),
  stringsAsFactors = FALSE
)

decoded <- .ablation_decode_direct_features(
  reference_d1 = reference_d1,
  query_d1 = query_d1,
  reference_direct = reference_direct,
  query_direct = query_direct,
  feature_manifest = feature_manifest,
  blocks = list(M1 = 1:2),
  rank = 2L,
  lambda = 0.01
)

stopifnot(decoded$status == "complete")
stopifnot(nrow(decoded$per_feature) == 3L)
stopifnot(nrow(decoded$summary) == 3L)
stopifnot(all(c(
  "balanced_accuracy",
  "brier",
  "spearman",
  "mae",
  "rmse"
) %in% colnames(decoded$per_feature)))
gene_pair <- decoded$per_feature[
  decoded$per_feature$feature_type == "gene_pair",
  ,
  drop = FALSE
]
single_bin <- decoded$per_feature[
  decoded$per_feature$feature_type == "single_bin",
  ,
  drop = FALSE
]
set_pair <- decoded$per_feature[
  decoded$per_feature$feature_type == "set_pair",
  ,
  drop = FALSE
]
stopifnot(gene_pair$balanced_accuracy > 0.8)
stopifnot(gene_pair$brier < 0.2)
stopifnot(single_bin$spearman > 0.8)
stopifnot(single_bin$mae < 0.5)
stopifnot(set_pair$spearman > 0.9)
stopifnot(set_pair$rmse < 0.5)

message("13-test-ablation-decoder: all tests passed")
