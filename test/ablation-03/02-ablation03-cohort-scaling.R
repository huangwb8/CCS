#!/usr/bin/env Rscript

# Purpose: Recompute only the two-dimensional cohort-evidence scaling extension.
# Input: Stage-01 objects and existing ablation outputs. External d1 is read
# directly from the full precomputed resCCS object.
# Parameters: Shared experiment params with breadth/depth repeats and matched banks.
# Output: Schema-v2 cohort_scaling.rds plus synchronized manifest/result references.

env_candidates <- c(
  file.path(getwd(), "00.Environment.R"),
  file.path(getwd(), "test", "ablation-03", "00.Environment.R")
)
env_path <- env_candidates[file.exists(env_candidates)][1L]
if (is.na(env_path)) stop("ablation-03: run from the project directory or repository root.", call. = FALSE)
source(env_path)
source(.ablation03_path("01-ablation03-test-data.R"))
source(.ablation03_path("02-ablation03-experiment_functions.R"))
source(.ablation03_repo_path("R", "ablation.R"))

seed <- 20260805L
output_dir <- .ablation03_path("tmp", "ablation-experiment")
params <- .ae_ablation_params(filtered_cohorts, n_cores)
config <- .ablation_resolve_representation_config(seed, params)
analysis <- .ablation_prepare_representation_analysis(
  object = resCCS,
  d1_source = resCCS_full,
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
