# Verify that ablation-03 has one targets entry point and no legacy runner.
stopifnot(file.exists("test/ablation-03/_targets.R"))
stopifnot(file.exists("test/ablation-03/scripts/run-targets-renv.ps1"))
stopifnot(!file.exists("test/ablation-03/run-ablation-03.R"))
stopifnot(!file.exists("test/ablation-03/scripts/run-fresh-analysis.ps1"))

target_lines <- readLines("test/ablation-03/_targets.R", warn = FALSE, encoding = "UTF-8")
launcher_lines <- readLines(
  "test/ablation-03/scripts/run-targets-renv.ps1",
  warn = FALSE,
  encoding = "UTF-8"
)
stopifnot(
  any(grepl("tar_make", launcher_lines, fixed = TRUE)),
  any(grepl("tar_option_set", target_lines, fixed = TRUE)),
  any(grepl("CCS_ABLATION_INPUT_RDS", launcher_lines, fixed = TRUE)),
  any(grepl(".local\\bin", launcher_lines, fixed = TRUE)),
  !any(grepl(
    "CCS_ABLATION_PROFILE|[-]Profile|profile[ ]*=",
    launcher_lines,
    ignore.case = TRUE
  )),
  !any(grepl("system2", target_lines, fixed = TRUE))
)

cat("targets run boundary tests passed\n")
