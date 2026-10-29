# Isolated execution of the same benchmark and render entries, with smaller inputs.
root <- normalizePath(".", winslash = "/")
run_root <- file.path(root, "tmp/tests", paste0("normCCS-", format(Sys.time(), "%Y%m%d-%H%M%S")))
dir.create(run_root, recursive = TRUE)
products <- file.path(run_root, "products")
reports <- file.path(run_root, "reports")
dir.create(reports)
subject_files <- c("R/phenotype.R", "R/subtype_norm_dp.R",
                   "test/normCCS-5-more/benchmark.R", "test/normCCS-5-more/fixtures.R",
                   "test/normCCS-5-more/report.Rmd", "test/normCCS-5-more/render-report.R",
                   "test/normCCS-5-more/benchmark-config.json", "test/normCCS-5-more/00.Environment.R",
                   "test/normCCS-5-more/reference/subtype-norm-0.8.4.R",
                   "test/normCCS-5-more/reference/quality-single-cut.csv")
subject_before <- tools::md5sum(subject_files)
formal_snapshot <- function() {
  report <- "test/normCCS-5-more/report.html"
  paths <- sort(c(list.files("test/normCCS-5-more/products", recursive = TRUE, full.names = TRUE),
                  list.files("test/normCCS-5-more/reports", recursive = TRUE, full.names = TRUE),
                  if (file.exists(report)) report))
  list(paths = paths, hashes = tools::md5sum(paths), metadata = file.info(paths)[c("size", "mtime")])
}
formal_before <- formal_snapshot()
raw_metadata <- function() {
  paths <- list.files("raw", recursive = TRUE, full.names = TRUE)
  digest::digest(if (length(paths)) file.info(paths)[c("size", "mtime")] else character())
}
raw_before <- raw_metadata()
cfg <- jsonlite::fromJSON("test/normCCS-5-more/benchmark-config.json")
cfg$quality_m <- 6L
cfg$quality_c <- 3L
cfg$quality_seeds <- 1L
cfg$quality_k <- c(3L, 5L)
cfg$timing_m <- 6L
cfg$timing_c <- 3L
cfg$timing_k <- 5L
cfg$timing_reps <- 1L
cfg$scaling_m <- 6L
cfg$scaling_c <- 3L
cfg$scaling_k <- 5L
cfg$support_m <- 6L
cfg$support_c <- 3L
cfg$support_k <- 5L
config_file <- file.path(run_root, "config.json")
jsonlite::write_json(cfg, config_file, auto_unbox = TRUE, pretty = TRUE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
status <- system2(rscript, c("--vanilla", "test/normCCS-5-more/benchmark.R",
                             shQuote(products), shQuote(paste0("--config=", config_file))),
                  stdout = file.path(run_root, "benchmark.log"), stderr = file.path(run_root, "benchmark-errors.log"))
stopifnot(status == 0L)
quality <- read.csv(file.path(products, "quality.csv"))
timing <- read.csv(file.path(products, "timings.csv"))
scaling <- read.csv(file.path(products, "scaling.csv"))
stopifnot(nrow(quality) == 4L, nrow(timing) == 4L, nrow(scaling) == 2L,
          all(is.finite(quality$gap)), all(abs(quality$gap) <= 1),
          all(scaling$median_auc >= 0 & scaling$median_auc <= 1))
old_products <- Sys.getenv("BENSZ_NORM_PRODUCTS", unset = NA_character_)
Sys.setenv(BENSZ_NORM_PRODUCTS = products)
report_path <- file.path(reports, "report.html")
status <- system2(rscript, c("--vanilla", "test/normCCS-5-more/render-report.R", shQuote(report_path)),
                  stdout = file.path(run_root, "render.log"), stderr = file.path(run_root, "render-errors.log"))
if (is.na(old_products)) Sys.unsetenv("BENSZ_NORM_PRODUCTS") else Sys.setenv(BENSZ_NORM_PRODUCTS = old_products)
stopifnot(status == 0L, file.exists(report_path), file.info(report_path)$size > 20000L,
          file.exists(file.path(reports, "figures/speedup.pdf")),
          identical(unname(tools::md5sum(subject_files)), unname(subject_before)),
          identical(formal_snapshot(), formal_before),
          identical(raw_metadata(), raw_before))
record <- list(project_state = "existing", workflow_mode = "simple", test_style = "synthetic_fixture",
               formal_entrypoint = c("benchmark.R", "render-report.R"),
               run_root = substring(run_root, nchar(root) + 2L),
               subject_identity = as.list(subject_before), preflight = "PASS",
               lightweight_execution = "PASS", full_data_execution = "NOT_RUN",
               command_exit_status = status, raw_metadata_unchanged = TRUE,
               formal_products_unchanged = TRUE, code_unchanged = TRUE)
jsonlite::write_json(record, file.path(run_root, "verification.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Isolated benchmark and report passed:", record$run_root, "\n")
