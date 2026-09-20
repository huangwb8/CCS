# Regression test: loading the analysis environment must not open a graphics device.

env_path <- file.path("test", "ablation-03", "00.Environment.R")
if (!file.exists(env_path)) {
  stop("Expected ablation-03 environment file is missing.", call. = FALSE)
}

old_cache_root <- Sys.getenv("CCS_ABLATION_CACHE_ROOT", unset = NA_character_)
Sys.setenv(CCS_ABLATION_CACHE_ROOT = tempfile("ablation-env-test-"))
on.exit({
  if (is.na(old_cache_root)) Sys.unsetenv("CCS_ABLATION_CACHE_ROOT") else {
    Sys.setenv(CCS_ABLATION_CACHE_ROOT = old_cache_root)
  }
}, add = TRUE)

device_requests <- 0L
original_device <- getOption("device")
on.exit(options(device = original_device), add = TRUE)
options(device = function(...) {
  device_requests <<- device_requests + 1L
  stop("Unexpected default graphics device request.", call. = FALSE)
})

source(env_path)

stopifnot(device_requests == 0L)
stopifnot(
  is.character(getOption("bensz.base_family")),
  nzchar(getOption("bensz.base_family"))
)

grDevices::pdf(file = NULL)
on.exit(grDevices::dev.off(), add = TRUE)
invisible(bensz_setup_plot_font(prefer = "sans"))
stopifnot(identical(graphics::par("family"), "sans"))
cat("environment graphics-device contract test passed\n")
