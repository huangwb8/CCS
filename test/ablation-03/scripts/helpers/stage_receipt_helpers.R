# Generic stage provenance helpers; no scientific computation.

.ae_atomic_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(value, temporary, version = 3)
  invisible(readRDS(temporary))
  backup <- paste0(path, ".bak-", Sys.getpid())
  if (file.exists(path) && !file.rename(path, backup)) {
    stop("Cannot stage provenance file for replacement: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    stop("Cannot commit provenance file: ", path, call. = FALSE)
  }
  if (file.exists(backup)) unlink(backup, force = TRUE)
  invisible(path)
}

.ae_write_stage_receipt <- function(stage_dir, inputs, outputs) {
  paths <- unique(normalizePath(c(inputs, outputs), winslash = "/", mustWork = TRUE))
  receipt <- list(
    schema_version = 2L,
    hashes = tools::md5sum(paths)
  )
  .ae_atomic_save_rds(receipt, file.path(stage_dir, "stage-receipt.rds"))
  invisible(receipt)
}

.ae_validate_stage_receipt <- function(stage_dir) {
  path <- file.path(stage_dir, "stage-receipt.rds")
  if (!file.exists(path)) stop("Missing stage receipt; rerun analysis: ", stage_dir, call. = FALSE)
  receipt <- tryCatch(readRDS(path), error = function(error) NULL)
  if (!is.list(receipt) || !identical(receipt$schema_version, 2L) ||
      is.null(receipt$hashes) || is.null(names(receipt$hashes))) {
    stop("Invalid stage receipt; rerun analysis: ", stage_dir, call. = FALSE)
  }
  current <- tools::md5sum(names(receipt$hashes))
  if (anyNA(current) || !identical(current, receipt$hashes)) {
    stop("Stale or mixed analysis products; rerun stage: ", stage_dir, call. = FALSE)
  }
  invisible(TRUE)
}

