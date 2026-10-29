# Same benchmark and report entries, with isolated small synthetic inputs.
root <- normalizePath(".", winslash = "/")
run_root <- file.path(root, "tmp/tests", paste0("norm-objectives-", format(Sys.time(), "%Y%m%d-%H%M%S")))
dir.create(run_root, recursive = TRUE)
products <- file.path(run_root, "products")
reports <- file.path(run_root, "reports")
dir.create(products)
dir.create(reports)
subject <- c("DESCRIPTION", "R/phenotype.R", "R/subtype_norm_dp.R", "R/subtype_norm_objective.R",
             "test/normCCS-5-more/benchmark-objectives.R", "test/normCCS-5-more/objective-report.Rmd",
             "test/normCCS-5-more/render-report.R", "test/normCCS-5-more/fixtures.R",
             "test/normCCS-5-more/objective-benchmark-config.json",
             "test/normCCS-5-more/00.Environment.R",
             "test/normCCS-5-more/reference/subtype-norm-0.8.5.R",
             "test/normCCS-5-more/reference/subtype-norm-dp-0.8.5.R")
identity <- tools::md5sum(subject)
snapshot <- function() {
  paths <- sort(c(list.files("test/normCCS-5-more/products", recursive = TRUE, full.names = TRUE),
                  list.files("test/normCCS-5-more/reports", recursive = TRUE, full.names = TRUE),
                  list.files("test/normCCS-5-more", pattern = "\\.html$", full.names = TRUE)))
  list(paths = paths, hashes = tools::md5sum(paths), metadata = file.info(paths)[c("size", "mtime")])
}
raw_snapshot <- function() {
  paths <- list.files("raw", recursive = TRUE, full.names = TRUE)
  if (length(paths)) file.info(paths)[c("size", "mtime")] else character()
}
before <- snapshot()
raw_before <- raw_snapshot()
cfg <- jsonlite::fromJSON("test/normCCS-5-more/objective-benchmark-config.json")
cfg$quality_m <- 6L
cfg$quality_c <- 3L
cfg$quality_seeds <- 1L
cfg$quality_k <- 4L
cfg$timing_m <- 6L
cfg$timing_c <- 3L
cfg$timing_k <- 4L
cfg$timing_reps <- 1L
cfg$scaling_m <- 6L
cfg$scaling_c <- 3L
cfg$scaling_k <- 4L
config <- file.path(run_root, "config.json")
jsonlite::write_json(cfg, config, auto_unbox = TRUE, pretty = TRUE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
status <- system2(rscript, c("--vanilla", "test/normCCS-5-more/benchmark-objectives.R", shQuote(products),
                             shQuote(paste0("--config=", config))),
                  stdout = file.path(run_root, "benchmark.log"), stderr = file.path(run_root, "benchmark-errors.log"))
stopifnot(status == 0L)
q <- read.csv(file.path(products, "objective-quality.csv"))
t <- read.csv(file.path(products, "objective-timing.csv"))
s <- read.csv(file.path(products, "objective-scaling.csv"))
stopifnot(nrow(q) == 10L, nrow(t) == 8L, nrow(s) == 10L,
          all(is.finite(q$gap)), all(q$gap >= -1e-10),
          all(q$objective_score >= q$initial_score - 1e-12),
          all(s$cohort_mean >= 0 & s$cohort_mean <= 1))
old <- Sys.getenv("BENSZ_NORM_PRODUCTS", unset = NA_character_)
Sys.setenv(BENSZ_NORM_PRODUCTS = products)
report <- file.path(reports, "objective-report.html")
status <- system2(rscript, c("--vanilla", "test/normCCS-5-more/render-report.R", shQuote(report), "--objectives"),
                  stdout = file.path(run_root, "render.log"), stderr = file.path(run_root, "render-errors.log"))
if (is.na(old)) Sys.unsetenv("BENSZ_NORM_PRODUCTS") else Sys.setenv(BENSZ_NORM_PRODUCTS = old)
stopifnot(status == 0L, file.info(report)$size > 20000,
          file.exists(file.path(reports, "figures/objective-gap.pdf")),
          identical(snapshot(), before), identical(raw_snapshot(), raw_before),
          identical(tools::md5sum(subject), identity))
jsonlite::write_json(list(project_state = "existing", workflow_mode = "simple", test_style = "synthetic_fixture",
                          formal_entrypoint = c("benchmark-objectives.R", "render-report.R --objectives"),
                          subject_identity = as.list(identity), preflight = "PASS", lightweight_execution = "PASS",
                          full_data_execution = "NOT_RUN", command_exit_status = status,
                          formal_products_unchanged = TRUE, raw_metadata_unchanged = TRUE),
                     file.path(run_root, "verification.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Isolated objective benchmark and report passed:", substring(run_root, nchar(root) + 2L), "\n")
