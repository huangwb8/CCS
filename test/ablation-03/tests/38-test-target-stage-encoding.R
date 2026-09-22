# R 4.3.1 `sys.source()` has no `encoding` argument. The targets stage
# loader must use `source()` while preserving its isolated environment.
stopifnot(
  !"encoding" %in% names(formals(base::sys.source)),
  "encoding" %in% names(formals(base::source))
)

target_lines <- readLines(
  "test/ablation-03/targets/functions.R",
  warn = FALSE,
  encoding = "UTF-8"
)
stopifnot(
  any(grepl(
    'source(stage_path, local = stage_env, encoding = "UTF-8")',
    target_lines,
    fixed = TRUE
  )),
  !any(grepl("sys.source(stage_path", target_lines, fixed = TRUE))
)

stage_file <- tempfile("ablation-03-stage-", fileext = ".R")
on.exit(unlink(stage_file, force = TRUE), add = TRUE)
writeLines(
  "stage_only_value <- 1L",
  stage_file,
  useBytes = TRUE
)

stage_env <- new.env(parent = globalenv())
source(stage_file, local = stage_env, encoding = "UTF-8")
stopifnot(
  identical(stage_env$stage_only_value, 1L),
  !exists("stage_only_value", envir = .GlobalEnv, inherits = FALSE)
)

cat("targets stage UTF-8 loader passed\n")
