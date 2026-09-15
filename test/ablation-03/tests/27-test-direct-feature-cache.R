#!/usr/bin/env Rscript

# Verify the persistent Direct-GSClassifier cache helper: matching inputs hit,
# changed expression content misses, and corrupted files are ignored.
source(file.path("R", "ccs.R"))
source(file.path("R", "ablation.R"))

object <- methods::new(
  "CCS",
  Repeat = list(
    method = "GSClassifier",
    geneSet = list(A = c("g1", "g2")),
    geneAnnotation = data.frame(ENSEMBL = c("g1", "g2", "g3")),
    geneid = "ensembl"
  ),
  Data = list(Probability = list(d1 = matrix(numeric(), 0L, 0L)))
)
expr <- matrix(
  c(1, 2, 3, 4, 5, 6),
  nrow = 3L,
  dimnames = list(c("g1", "g2", "g3"), c("s1", "s2"))
)
feature_manifest <- list(
  features = c("g1", "g1:g2"),
  feature_type = c("single_bin", "gene_pair"),
  break_vec = c(0, 1)
)
sample_ids <- colnames(expr)
key <- .ablation_direct_feature_cache_key(
  object, expr, feature_manifest, sample_ids
)
cache_path <- file.path(tempdir(), "ccs-direct-feature-cache-test.rds")
direct <- matrix(
  c(1, 2, 3, 4),
  nrow = 2L,
  dimnames = list(sample_ids, feature_manifest$features)
)
.ablation_atomic_save_rds(
  list(
    schema_version = 1L,
    key = key,
    sample_ids = sample_ids,
    feature_manifest = feature_manifest,
    direct = direct
  ),
  cache_path
)
stopifnot(identical(
  .ablation_read_direct_feature_cache(
    cache_path, key, sample_ids, feature_manifest
  ),
  direct
))

expr_changed <- expr
expr_changed[1L, 1L] <- 99
key_changed <- .ablation_direct_feature_cache_key(
  object, expr_changed, feature_manifest, sample_ids
)
stopifnot(is.null(
  .ablation_read_direct_feature_cache(
    cache_path, key_changed, sample_ids, feature_manifest
  )
))

writeLines("corrupted", cache_path)
stopifnot(is.null(
  .ablation_read_direct_feature_cache(
    cache_path, key, sample_ids, feature_manifest
  )
))
cat("direct feature cache test passed\n")
