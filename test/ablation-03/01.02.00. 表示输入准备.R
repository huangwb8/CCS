# Prepare Direct/d1 with the original functions, parameters and stage-specific seeds.
bootstrap <- c(file.path("scripts", "helpers", "workflow_helpers.R"),
  file.path("test", "ablation-03", "scripts", "helpers", "workflow_helpers.R"))
bootstrap <- bootstrap[file.exists(bootstrap)][1L]
if (is.na(bootstrap)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(bootstrap, local = TRUE)
inputs <- .wf_read("01-data", "inputs.rds")
source(.ablation03_repo_path("R", "ablation.R"))
source(.ablation03_path("02.01.00. 表示分析_functions.R"))
ablation_params <- .ae_ablation_params(inputs$filtered_cohorts, inputs$n_cores)
stage_parameters <- list(
  params = ablation_params,
  representation_seed = 20260805L,
  structural_seed = 20260912L
)
if (.wf_cache_hit("01-representations", stage_parameters)) quit(save = "no", status = 0L)
.wf_begin("01-representations", stage_parameters)
output_dir <- .wf_output("01-representations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(.wf_output("ablation-experiment"), recursive = TRUE, showWarnings = FALSE)
for (stage in c("representation", "structural")) {
  seed <- if (stage == "representation") 20260805L else 20260912L
  config <- .ablation_resolve_representation_config(seed, ablation_params)
  analysis <- .ablation_prepare_representation_analysis(
    object = inputs$resCCS_ablation, data = inputs$data_all,
    metadata = inputs$ablation_metadata, config = config,
    output.dir = .wf_output("ablation-experiment"), seed = seed, verbose = TRUE)
  saveRDS(list(analysis = analysis, config = config, params = ablation_params,
    seed = seed, structural = if (stage == "structural") list(
      full_d1 = inputs$full_d1,
      full_module_manifest = .ablation_module_manifest(inputs$resCCS_full),
      tissue_resolution_audit = inputs$tissue_resolution_audit,
      ablation_metadata = inputs$ablation_metadata) else NULL),
    file.path(output_dir, paste0(stage, "-inputs.rds")))
  if (stage == "representation") {
    saveRDS(list(reference = analysis$prepared$reference_metadata,
      query = analysis$prepared$query_metadata), file.path(output_dir, "sample-contract.rds"))
  }
}
.wf_receipt("01-representations", "01.02.00. 表示输入准备",
  inputs = c(.wf_output("01-data", "stage-receipt.rds"),
    .ablation03_repo_path("R", "ablation.R"),
    .ablation03_path("02.01.00. 表示分析_functions.R"),
    list.files(inputs$resCCS_ablation@Repeat$model.dir, pattern = "modelFit.rds$",
      recursive = TRUE, full.names = TRUE)),
  outputs = list.files(
    output_dir,
    pattern = "^(representation-inputs|structural-inputs|sample-contract)\\.rds$",
    full.names = TRUE
  ))
