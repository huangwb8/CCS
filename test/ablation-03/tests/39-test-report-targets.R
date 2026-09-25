# Verify that every ablation-03 Rmd is rendered by the formal targets graph.
project_dir <- file.path(getwd(), "test", "ablation-03")
old_dir <- setwd(project_dir)
on.exit(setwd(old_dir), add = TRUE)
source(file.path("renv", "activate.R"))

manifest <- targets::tar_manifest(script = "_targets.R")
stopifnot(
  'input_rds' %in% manifest$name,
  grepl('input_rds', manifest$command[match('data_preparation', manifest$name)], fixed = TRUE)
)
expected <- c(
  "data_overview_report",
  "representation_report",
  "biology_report",
  "structural_report"
)
source_targets <- paste0(expected, "_sources")
stopifnot(
  all(c(expected, source_targets) %in% manifest$name),
  all(grepl(
    ".ablation03_render_report",
    manifest$command[match(expected, manifest$name)],
    fixed = TRUE
  )),
  all(vapply(seq_along(expected), function(i) {
    grepl(source_targets[i], manifest$command[match(expected[i], manifest$name)], fixed = TRUE)
  }, logical(1)))
)
inference_targets <- c("statistical_inference", "learning_query_inference",
  "geometry_sensitivity")
stopifnot(
  all(inference_targets %in% manifest$name),
  all(grepl("statistical_inference",
    manifest$command[match(c("representation_report", "biology_report",
      "structural_report"), manifest$name)], fixed = TRUE)),
  all(vapply(c("learning_query_inference", "geometry_sensitivity"),
    function(target) grepl(target,
      manifest$command[match("representation_report", manifest$name)], fixed = TRUE),
    logical(1)))
)

target_lines <- readLines("_targets.R", warn = FALSE, encoding = "UTF-8")
launcher_lines <- readLines(
  file.path("scripts", "run-targets-renv.ps1"),
  warn = FALSE,
  encoding = "UTF-8"
)
stopifnot(
  sum(grepl('format = "file"', target_lines, fixed = TRUE)) >=
    length(c(expected, source_targets)),
  !any(grepl("rmarkdown::render", launcher_lines, fixed = TRUE))
)

cat("targets report wiring tests passed\n")
