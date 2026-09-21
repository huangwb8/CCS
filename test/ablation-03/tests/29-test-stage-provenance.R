# Stage provenance belongs to targets, not numbered-script receipts.
stopifnot(!file.exists("test/ablation-03/scripts/helpers/checkpoint_helpers.R"))
stopifnot(!file.exists("test/ablation-03/scripts/helpers/stage_receipt_helpers.R"))
stopifnot(!file.exists("test/ablation-03/scripts/helpers/run_identity_helpers.R"))
target_lines <- readLines("test/ablation-03/_targets.R", warn = FALSE, encoding = "UTF-8")
stopifnot(
  any(grepl("targets::tar_target", target_lines, fixed = TRUE)),
  any(grepl("representation_analysis", target_lines, fixed = TRUE)),
  any(grepl("biology_analysis", target_lines, fixed = TRUE)),
  any(grepl("structural_analysis", target_lines, fixed = TRUE))
)
cat("targets stage provenance boundary passed\n")
