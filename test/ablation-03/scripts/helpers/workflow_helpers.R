# Compatibility helpers for scientific stage code executed by targets.
#
# targets owns dependency tracking, cache invalidation, recovery, locking and
# parallel scheduling. These helpers only resolve paths and provide ordinary
# file I/O needed by the scientific calculations retained in numbered R files.

env_path <- c("00.Environment.R", "test/ablation-03/00.Environment.R")
env_path <- env_path[file.exists(env_path)][1L]
if (is.na(env_path)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(env_path, local = TRUE)

.ablation03_ccs_version <- as.character(utils::packageVersion("CCS"))
if (!identical(.ablation03_ccs_version, "0.8.3")) {
  stop("ablation-03 requires installed CCS 0.8.3; found ", .ablation03_ccs_version, ".", call. = FALSE)
}
.ablation03_ccs_namespace <- asNamespace("CCS")
.ablation03_ccs_symbols <- grep("^\\.ablation_", ls(.ablation03_ccs_namespace, all.names = TRUE), value = TRUE)
for (.ablation03_symbol in .ablation03_ccs_symbols) {
  assign(.ablation03_symbol, get(.ablation03_symbol, envir = .ablation03_ccs_namespace, inherits = FALSE), envir = environment())
}
rm(.ablation03_symbol)
.ablation03_ccs_description <- file.path(find.package("CCS"), "DESCRIPTION")

.wf_root <- .ablation03_dir
.wf_path <- function(...) .ablation03_path(...)
.wf_output <- function(...) file.path(.ablation03_cache_root, ...)
.wf_mode <- "targets"
.wf_targets_mode <- TRUE

.wf_atomic_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(value, temporary, version = 3)
  if (!file.rename(temporary, path)) stop("Cannot finalize target artifact: ", path, call. = FALSE)
  invisible(path)
}

.wf_read_rds_safe <- function(path, validator = NULL) {
  if (!file.exists(path)) return(NULL)
  value <- tryCatch(readRDS(path), error = function(error) NULL)
  if (is.null(value)) return(NULL)
  if (!is.null(validator) && !isTRUE(validator(value))) return(NULL)
  value
}

.wf_cache_hit <- function(stage, parameters = list()) FALSE

.wf_begin <- function(stage, parameters = list()) {
  dir.create(.wf_output(stage), recursive = TRUE, showWarnings = FALSE)
  invisible(.wf_output(stage))
}

.wf_validate <- function(stage_dir) {
  if (!dir.exists(stage_dir)) stop("Missing target artifact directory: ", stage_dir, call. = FALSE)
  invisible(TRUE)
}

.wf_read <- function(stage, filename) {
  path <- file.path(.wf_output(stage), filename)
  if (!file.exists(path)) stop("Missing target artifact: ", path, call. = FALSE)
  readRDS(path)
}

# Stage receipts and project checkpoint markers belonged to the removed
# numbered-script runner. Keep a no-op call boundary while scientific scripts
# are reduced to pure target commands.
.wf_receipt <- function(stage, unit_stem, inputs = character(), outputs = character()) invisible(NULL)
.wf_fail <- function(stage, message) invisible(NULL)
