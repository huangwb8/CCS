#!/usr/bin/env Rscript

# Purpose: Recompute only the two-dimensional cohort-evidence scaling extension.
# Input: Stage-01 objects, the verified external-d1 cache, and existing ablation outputs.
# Parameters: Shared experiment params with breadth/depth repeats and matched banks.
# Output: Schema-v2 cohort_scaling.rds plus synchronized manifest/result references.

source(file.path("test", "ablation-02", "00.Environment.R"))
source(file.path("test", "ablation-02", "01-ablation-test-data.R"))
source(file.path("test", "ablation-02", "02-ablation-experiment_functions.R"))
source(file.path("R", "ablation.R"))

seed <- 20260805L
output_dir <- file.path(
  "test", "ablation-02", "tmp", "ablation-experiment"
)
params <- .ae_ablation_params(filtered_cohorts, n_cores)
config <- .ablation_resolve_representation_config(seed, params)
analysis <- .ablation_prepare_representation_analysis(
  object = resCCS,
  data = data_all,
  metadata = ablation_metadata,
  config = config,
  output.dir = output_dir,
  seed = seed,
  verbose = TRUE
)
prepared <- analysis$prepared
cohort_scaling <- .ablation_representation_scaling(
  prepared = prepared,
  config = config,
  label_column = analysis$anchor,
  seed = seed + 35000L,
  verbose = TRUE,
  cache_path = file.path(output_dir, "cohort-scaling-fit-cache.rds")
)

saveRDS(
  cohort_scaling,
  file.path(output_dir, "cohort_scaling.rds")
)

manifest_path <- file.path(output_dir, "manifest.rds")
manifest <- readRDS(manifest_path)
feature_types <- prepared$feature_manifest$feature_manifest$feature_type
manifest$version <- max(5L, manifest$version)
manifest$scaling_created <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
manifest$gene_signature_count <- sum(feature_types == "single_bin")
manifest$scaling_direct_feature_count <-
  cohort_scaling$diagnostics$direct_contracts$feature_count[
    cohort_scaling$diagnostics$direct_contracts$contract_role == "main"
  ][1]
manifest$scaling_schema_version <- cohort_scaling$schema_version
manifest$config <- config
manifest$config_hash <- digest::digest(config, algo = "md5")
saveRDS(manifest, manifest_path)

result_path <- file.path(output_dir, "ablation-result.rds")
ablation_result <- readRDS(result_path)
ablation_result$manifest <- manifest
ablation_result$cohort_scaling <- cohort_scaling
saveRDS(ablation_result, result_path)

luckyBase::LuckyVerbose(
  "02-ablation-cohort-scaling: complete; output = ",
  normalizePath(
    file.path(output_dir, "cohort_scaling.rds"),
    winslash = "/",
    mustWork = TRUE
  )
)
