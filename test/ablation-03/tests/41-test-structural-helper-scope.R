if (.Platform$OS.type == "windows") {
  stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")))
}
project_dir <- file.path(getwd(), "test", "ablation-03")
setwd(project_dir)
script_file <- list.files(project_dir, pattern = "^02[.]03[.]00.*[.]R$", full.names = TRUE)
script_file <- script_file[!grepl("_functions[.]R$", script_file)]
functions_file <- list.files(project_dir, pattern = "^02[.]03[.]00.*_functions[.]R$", full.names = TRUE)
stopifnot(length(script_file) == 1L, length(functions_file) == 1L)
script_lines <- readLines(script_file, warn = FALSE, encoding = "UTF-8")
stopifnot(
  any(grepl("_functions.R", script_lines, fixed = TRUE)),
  any(grepl("local = TRUE", script_lines, fixed = TRUE)),
  any(grepl(".ablation_scale_train_apply <- getFromNamespace", script_lines, fixed = TRUE)),
  any(grepl(".ablation_module_balanced_transform <- getFromNamespace", script_lines, fixed = TRUE))
)

source(file.path(project_dir, "renv", "activate.R"))
stage_env <- new.env(parent = globalenv())
source(functions_file, local = stage_env, encoding = "UTF-8")
stage_env$.ablation_scale_train_apply <- getFromNamespace(".ablation_scale_train_apply", "CCS")
stage_env$.ablation_module_balanced_transform <- getFromNamespace(
  ".ablation_module_balanced_transform", "CCS"
)
for (name in c(".asr_prepare_direction_base", ".asr_evaluate_direction")) {
  helper_env <- environment(get(name, envir = stage_env))
  stopifnot(
    identical(helper_env, stage_env),
    is.function(get(".ablation_scale_train_apply", envir = helper_env)),
    is.function(get(".ablation_module_balanced_transform", envir = helper_env))
  )
}
cat("isolated structural helper scope passed\n")
