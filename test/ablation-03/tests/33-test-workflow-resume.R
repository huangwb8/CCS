# The old numbered-script resume protocol was removed. targets is now the only
# owner of cache invalidation and recovery; this test guards the new boundary.
target_script <- "test/ablation-03/_targets.R"
stopifnot(file.exists(target_script))
lines <- readLines(target_script, warn = FALSE, encoding = "UTF-8")
stopifnot(
  any(grepl("tar_option_set", lines, fixed = TRUE)),
  any(grepl("controller = controller", lines, fixed = TRUE)),
  any(grepl("representation_inputs", lines, fixed = TRUE)),
  any(grepl("biology_analysis", lines, fixed = TRUE)),
  any(grepl("structural_analysis", lines, fixed = TRUE))
)
stopifnot(!file.exists("test/ablation-03/run-ablation-03.R"))
stopifnot(!any(grepl("system2", lines, fixed = TRUE)))

cat("targets workflow boundary tests passed\n")
