#!/usr/bin/env Rscript

# Purpose: Run the complete d1 vs Direct-GSClassifier ablation on real CCS data.
# Input: Objects prepared by 01-ablation03-test-data.R and functions from R/ablation.R.
# Parameters: External cohorts, paired learning curves, 2D bank scaling, and controls.
# Output: Full unfiltered result objects under the shared tmp/ablation-experiment/ cache.

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

output_dir <- .ablation03_path("tmp", "ablation-experiment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ablation_params <- .ae_ablation_params(filtered_cohorts, n_cores)

ablation_result <- ablation(
  object = resCCS_ablation,
  data = data_all,
  metadata = ablation_metadata,
  experiment = "representation",
  output.dir = output_dir,
  params = ablation_params,
  seed = 20260805L,
  verbose = TRUE
)

saveRDS(
  ablation_result,
  file.path(output_dir, "ablation-result.rds")
)

luckyBase::LuckyVerbose(
  "02-ablation-experiment: complete; evidence level = ",
  ablation_result$evidence_level,
  "; output = ",
  normalizePath(output_dir, winslash = "/", mustWork = TRUE)
)
