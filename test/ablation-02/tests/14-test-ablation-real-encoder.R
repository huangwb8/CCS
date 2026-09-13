#!/usr/bin/env Rscript

# Purpose: Verify that representation ablation consumes precomputed d1 only.
# Input: Real CCS data with a bounded reference/query slice.
# Output: Query rows are the intersection with object@Data$Probability$d1;
# missing rows are reported and never sent through a model encoder.

source(file.path("test", "ablation-02", "01-ablation-test-data.R"))
source(file.path("R", "ablation.R"))

module_manifest <- .ablation_module_manifest(resCCS)
module_ids <- module_manifest$modules$module_id[seq_len(3L)]
config <- .ablation_resolve_representation_config(
  seed = 1401L,
  params = list(
    comparison = list(module_ids = module_ids),
    provenance = list(max_reference_samples = 80L, max_query_samples = 40L),
    geometry = list(k = c(1L, 3L), search = "exact", geometry_samples = 80L),
    validation = list(enabled = FALSE),
    controls = list(null_rp = FALSE, null_perm = TRUE),
    output = list(cover = TRUE)
  )
)

prepared <- .ablation_prepare_representation_input(
  object = resCCS,
  data = data_all,
  metadata = ablation_metadata[
    ablation_metadata$sample_id %in%
      rownames(resCCS@Data$Probability$d1),
    ,
    drop = FALSE
  ],
  config = config,
  output.dir = tempfile("ablation-real-input-"),
  seed = 1401L,
  verbose = FALSE
)

d1_ids <- rownames(resCCS@Data$Probability$d1)
stopifnot(all(rownames(prepared$query_d1) %in% d1_ids))
stopifnot(all(!prepared$excluded_query_d1_ids %in% d1_ids))
stopifnot(identical(
  sort(rownames(prepared$query_d1)),
  sort(intersect(rownames(prepared$query_direct), d1_ids))
))

message("14-test-ablation-real-encoder: precomputed d1 contract passed")
