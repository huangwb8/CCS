#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

`%||%` <- function(x, y) if (is.null(x) || !nzchar(x)) y else x
file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else {
  sys.frames()[[1L]]$ofile %||% "scripts/renv-ablation03.R"
}
project <- normalizePath(
  file.path(dirname(normalizePath(script_file, mustWork = FALSE)), ".."),
  winslash = "/",
  mustWork = TRUE
)
command <- value_after("--command", "check")

activate <- function() {
  activation <- file.path(project, "renv", "activate.R")
  if (!file.exists(activation)) stop("ablation-03 renv is not initialized. Run --command init first.", call. = FALSE)
  source(activation, local = .GlobalEnv)
  if (!requireNamespace("renv", quietly = TRUE)) stop("The project renv library is unavailable after activation.", call. = FALSE)
  invisible(TRUE)
}

ensure_global_renv <- function() {
  if (requireNamespace("renv", quietly = TRUE)) return(invisible(TRUE))
  repos <- getOption("repos")
  if (is.null(repos) || identical(unname(repos[["CRAN"]]), "@CRAN@")) options(repos = c(CRAN = Sys.getenv("CRAN_REPOS", "https://cloud.r-project.org")))
  install.packages("renv")
  if (!requireNamespace("renv", quietly = TRUE)) stop("Unable to install renv into the active R library.", call. = FALSE)
  invisible(TRUE)
}

required_packages <- c("targets", "digest", "rmarkdown", "knitr", "yaml")
project_lock <- file.path(project, "renv.lock")
setwd(project)

if (identical(command, "init")) {
  ensure_global_renv()
  if (!file.exists(file.path(project, "renv", "activate.R"))) renv::init(project = project, bare = TRUE, restart = FALSE) else renv::activate(project = project)
  message("Initialized ablation-03 renv project at ", project)
  quit(status = 0L)
}

activate()
if (identical(command, "install-core")) {
  renv::install(required_packages)
  message("Installed core packages: ", paste(required_packages, collapse = ", "))
} else if (identical(command, "install-ccs")) {
  ccs_library <- Sys.getenv("CCS_ABLATION_CCS_LIBRARY", unset = "")
  if (!nzchar(ccs_library)) {
    candidates <- c(
      file.path(R.home(".."), "library"),
      file.path(R.home(), "library"),
      "C:/R/R-4.3.1/library"
    )
    candidates <- normalizePath(candidates[dir.exists(candidates)], winslash = "/", mustWork = FALSE)
    ccs_library <- candidates[vapply(candidates, function(x) dir.exists(file.path(x, "CCS")), logical(1))][1L]
  }
  if (is.na(ccs_library) || !nzchar(ccs_library) || !dir.exists(file.path(ccs_library, "CCS"))) {
    stop(
      "Set CCS_ABLATION_CCS_LIBRARY to a library containing the checked CCS 0.8.3 package.",
      call. = FALSE
    )
  }
  renv::hydrate(
    packages = c("CCS", "GSClassifier", "luckyBase"),
    sources = ccs_library,
    prompt = FALSE,
    project = project
  )
  message("Hydrated CCS and companion packages from ", ccs_library)
} else if (identical(command, "snapshot")) {
  renv::snapshot(project = project, prompt = FALSE, type = "all", force = TRUE)
  message("Wrote ", project_lock)
} else if (identical(command, "restore")) {
  if (!file.exists(project_lock)) stop("Missing renv.lock; run --command snapshot first.", call. = FALSE)
  renv::restore(project = project, prompt = FALSE)
} else if (identical(command, "status")) {
  if (!file.exists(project_lock)) stop("Missing renv.lock.", call. = FALSE)
  renv::status(project = project)
} else if (identical(command, "check")) {
  if (!file.exists(project_lock)) stop("Missing renv.lock.", call. = FALSE)
  package_version <- tryCatch(as.character(utils::packageVersion("CCS")), error = function(e) NA_character_)
  if (is.na(package_version) || package_version != "0.8.3") {
    found <- if (is.na(package_version)) "not installed" else package_version
    stop("ablation-03 requires installed CCS 0.8.3; found ", found, call. = FALSE)
  }
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  message("renv project: ", project)
  message("renv lock md5: ", unname(tools::md5sum(project_lock)))
  message("CCS version: ", package_version)
  message("CCS path: ", find.package("CCS"))
} else {
  stop("Unknown command: ", command, ". Use init, install-core, install-ccs, snapshot, restore, status or check.", call. = FALSE)
}
