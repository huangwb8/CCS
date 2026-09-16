# Run the original representation experiment on persisted inputs only.
bootstrap <- c("00-workflow_functions.R",
  "test/ablation-03/00-workflow_functions.R")
bootstrap <- bootstrap[file.exists(bootstrap)][1L]
if (is.na(bootstrap)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(bootstrap, local = TRUE)
bundle <- .wf_read("01-representations", "representation-inputs.rds")
source(.ablation03_repo_path("R", "ablation.R"))
output_dir <- .wf_output("ablation-experiment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ablation_result <- .ablation_run_prepared_representation(
  analysis = bundle$analysis, config = bundle$config,
  output.dir = output_dir, seed = bundle$seed, verbose = TRUE
)

saveRDS(
  ablation_result,
  file.path(output_dir, "ablation-result.rds")
)

saveRDS(.wf_read("01-data", "data-profile.rds"), file.path(output_dir, "data-profile.rds"))
.wf_receipt("ablation-experiment", "02-ablation03-representation.R",
  inputs = c(.wf_output("01-representations", "stage-receipt.rds"),
    .wf_output("01-data", "stage-receipt.rds")),
  outputs = file.path(output_dir, c("manifest.rds", "native_geometry.rds",
    "retrieval.rds", "anchor_retrieval.rds", "sample-contract.rds", "readout.rds",
    "learning_curve.rds", "cohort_scaling.rds", "tradeoffs.rds", "ablation-result.rds",
    "endpoint_eligibility.rds", "data-profile.rds")))
