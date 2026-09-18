# Optional scaling-only rerun; other calculations are unchanged.
bootstrap <- c(file.path("..", "templates", "workflow_helpers.R"),
  file.path("test", "ablation-03", "templates", "workflow_helpers.R"))
bootstrap <- bootstrap[file.exists(bootstrap)][1L]
if (is.na(bootstrap)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(bootstrap, local = TRUE)
bundle <- .wf_read("01-representations", "representation-inputs.rds")
source(.ablation03_repo_path("R", "ablation.R"))
seed <- bundle$seed
config <- bundle$config
analysis <- bundle$analysis
prepared <- analysis$prepared
output_dir <- .wf_output("ablation-experiment")
.wf_validate(output_dir)
previous_receipt <- readRDS(file.path(output_dir, "stage-receipt.rds"))
previous_manifest <- readRDS(file.path(output_dir, "manifest.rds"))
if (!identical(previous_manifest$input_key, prepared$input_key) ||
    !identical(previous_manifest$config, config)) {
  stop("Scaling-only update requires the same inputs/config as stage 02; rerun stage 02.", call. = FALSE)
}
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
.ae_write_stage_receipt(output_dir,
  c(names(previous_receipt$hashes), .wf_path("tools", "recompute-cohort-scaling.R")), character())

luckyBase::LuckyVerbose(
  "02-ablation-cohort-scaling: complete; output = ",
  normalizePath(
    file.path(output_dir, "cohort_scaling.rds"),
    winslash = "/",
    mustWork = TRUE
  )
)
