# Storage configuration only; this file never activates renv or runs analysis.
ccs_renv_paths <- function(environment_root) {
  environment_root <- normalizePath(environment_root, winslash = "/", mustWork = FALSE)
  Sys.setenv(
    RENV_PATHS_LIBRARY = file.path(environment_root, "library"),
    RENV_PATHS_CACHE = file.path(environment_root, "cache"),
    RENV_PATHS_SOURCE = file.path(environment_root, "source"),
    RENV_PATHS_BINARY = file.path(environment_root, "binary")
  )
  # Isolated tests may explicitly supply their own sandbox before activation.
  if (!nzchar(Sys.getenv("RENV_PATHS_SANDBOX"))) {
    Sys.setenv(RENV_PATHS_SANDBOX = file.path(environment_root, "sandbox"))
  }
  invisible(environment_root)
}
