# Run the original representation experiment on persisted inputs only.
bootstrap <- c(file.path("scripts", "helpers", "workflow_helpers.R"),
  file.path("test", "ablation-03", "scripts", "helpers", "workflow_helpers.R"))
bootstrap <- bootstrap[file.exists(bootstrap)][1L]
if (is.na(bootstrap)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(bootstrap, local = TRUE)
bundle <- .wf_read("01-representations", "representation-inputs.rds")
source(.ablation03_path("02.01.00. 表示分析_functions.R"))
bundle$config <- .ae_apply_runtime_config(bundle$config, bundle$analysis)
output_dir <- .wf_output("ablation-experiment")
stage_parameters <- list(
  seed = bundle$seed,
  config_hash = digest::digest(bundle$config, algo = "md5")
)
if (.wf_cache_hit("ablation-experiment", stage_parameters)) quit(save = "no", status = 0L)
.wf_begin("ablation-experiment", stage_parameters)
ablation_result <- .ablation_run_prepared_representation(
  analysis = bundle$analysis, config = bundle$config,
  output.dir = output_dir, seed = bundle$seed, verbose = TRUE
)

.ablation_atomic_save_rds(
  ablation_result,
  file.path(output_dir, "ablation-result.rds")
)

.ablation_atomic_save_rds(
  .wf_read("01-data", "data-profile.rds"),
  file.path(output_dir, "data-profile.rds")
)
.wf_receipt("ablation-experiment", "02.01.00. 表示分析",
  inputs = c(.wf_output("01-representations", "stage-receipt.rds"),
    .wf_output("01-data", "stage-receipt.rds"),
    .ablation03_ccs_description,
    .ablation03_path("02.01.00. 表示分析_functions.R")),
  outputs = file.path(output_dir, c("manifest.rds", "native_geometry.rds",
    "retrieval.rds", "anchor_retrieval.rds", "sample-contract.rds", "readout.rds",
    "learning_curve.rds", "cohort_scaling.rds", "tradeoffs.rds", "ablation-result.rds",
    "endpoint_eligibility.rds", "data-profile.rds")))
