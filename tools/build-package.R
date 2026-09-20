#!/usr/bin/env Rscript

# Build CCS from a clean, .Rbuildignore-aware staging tree. R 4.3.1 copies an
# external package path before applying .Rbuildignore, so invoking R CMD build
# on the repository itself can traverse large ignored reparse-point folders.

if (identical(.Platform$OS.type, "windows")) {
  invisible(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"))
}

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(NULL)
  args[[index + 1L]]
}
usage <- paste(
  "Usage: Rscript --vanilla tools/build-package.R",
  "--output-dir <directory>"
)
if (any(args %in% c("-h", "--help"))) {
  cat(usage, "\n")
  quit(save = "no", status = 0L)
}
output_dir <- value_after("--output-dir")
if (is.null(output_dir) || !nzchar(output_dir)) {
  stop("--output-dir is required.\n", usage, call. = FALSE)
}

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- all_args[grepl("^--file=", all_args)]
if (length(file_arg) == 0L) stop("Cannot resolve build script path.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
source_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Cannot create --output-dir: ", output_dir, call. = FALSE)

ignore_file <- file.path(source_root, ".Rbuildignore")
if (!file.exists(ignore_file)) stop("Missing .Rbuildignore: ", source_root, call. = FALSE)
patterns <- trimws(readLines(ignore_file, warn = FALSE, encoding = "UTF-8"))
patterns <- patterns[nzchar(patterns) & !startsWith(patterns, "#")]
ignored <- function(relative) {
  relative <- gsub("\\\\", "/", relative)
  relative == ".git" || startsWith(relative, ".git/") ||
    any(vapply(patterns, grepl, logical(1), x = relative, perl = TRUE))
}

staging_parent <- tempfile("ccs-build-stage-")
staging_root <- file.path(staging_parent, "CCS")
dir.create(staging_root, recursive = TRUE)
on.exit(unlink(staging_parent, recursive = TRUE, force = TRUE), add = TRUE)

copy_entry <- function(relative) {
  if (ignored(relative)) return(invisible(FALSE))
  source <- file.path(source_root, relative)
  destination <- file.path(staging_root, relative)
  if (dir.exists(source)) {
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    children <- list.files(source, all.files = TRUE, no.. = TRUE)
    for (child in children) copy_entry(file.path(relative, child))
    return(invisible(TRUE))
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, destination, overwrite = FALSE, copy.date = TRUE)) {
    stop("Cannot stage package file: ", relative, call. = FALSE)
  }
  invisible(TRUE)
}

entries <- list.files(source_root, all.files = TRUE, no.. = TRUE)
for (entry in entries) copy_entry(entry)

forbidden <- c(".bensz-api", ".ccs-cache", "renv", "test", "tmp")
present <- forbidden[file.exists(file.path(staging_root, forbidden))]
if (length(present) > 0L) {
  stop("Ignored package directories reached staging: ", paste(present, collapse = ", "), call. = FALSE)
}

build_output <- tempfile("ccs-build-output-")
dir.create(build_output)
on.exit(unlink(build_output, recursive = TRUE, force = TRUE), add = TRUE)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(build_output)
r_binary <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
status <- system2(r_binary, c("CMD", "build", shQuote(staging_root)))
if (!identical(status, 0L)) stop("R CMD build failed with status ", status, ".", call. = FALSE)

version <- read.dcf(file.path(source_root, "DESCRIPTION"), fields = "Version")[[1L]]
tarball_name <- paste0("CCS_", version, ".tar.gz")
built_tarball <- file.path(build_output, tarball_name)
tarball <- file.path(output_dir, tarball_name)
if (!file.exists(built_tarball)) {
  stop("Expected tarball was not created: ", built_tarball, call. = FALSE)
}
if (!file.copy(built_tarball, tarball, overwrite = TRUE, copy.date = TRUE)) {
  stop("Cannot copy package tarball to --output-dir: ", tarball, call. = FALSE)
}
cat(normalizePath(tarball, winslash = "/", mustWork = TRUE), "\n")
