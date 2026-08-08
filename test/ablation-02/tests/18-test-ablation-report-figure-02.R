#!/usr/bin/env Rscript

# Purpose: Verify Figure 2 explicitly displays the two sample-level overall means.
# Input: The ablation experiment R Markdown report source.
# Output: Source-contract assertions for the Figure 2 reference lines.

report_path <- file.path("test", "ablation-02", "02-ablation-experiment.Rmd")
report_source <- paste(readLines(report_path, warn = FALSE), collapse = "\n")

stopifnot(grepl("retrieval_reference", report_source, fixed = TRUE))
stopifnot(grepl("retrieval_mean_delta", report_source, fixed = TRUE))
stopifnot(grepl("retrieval_mrr_delta", report_source, fixed = TRUE))
stopifnot(grepl('linetype = "dashed"', report_source, fixed = TRUE))

message("18-test-ablation-report-figure-02: all tests passed")
