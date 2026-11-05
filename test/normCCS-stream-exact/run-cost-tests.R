# Run from the repository root with R 4.3.1; set R temporary directories to D:
# before starting R. Uses the current source adapter and installed native kernel.
repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
test_lib <- Sys.getenv("CCS_TEST_LIB",
                       "D:/cache/ccs/normCCS-stream-exact-088/environment/library")
.libPaths(c(test_lib, .libPaths()))
stopifnot(requireNamespace("CCS", quietly = TRUE))
env <- new.env(parent = asNamespace("CCS"))
sys.source(file.path(repo, "R/norm_exact.R"), envir = env)

cost <- env$norm_exact_search_cost(107L, 8L, 21L, "Raw", 16L)
stopifnot(cost$n_candidates == 24370067800,
          abs(cost$seconds - 409.035) < 1e-9, cost$explosion)
stopifnot(abs(env$norm_exact_search_cost(107L, 6L, 21L,
                                      "OptimizeRank", 16L)$seconds - 7.570) < 1e-9)
stopifnot(env$norm_exact_search_cost(107L, 8L, 21L, "Raw", 8L)$seconds ==
            2 * cost$seconds,
          env$norm_exact_search_cost(107L, 8L, 42L, "Raw", 16L)$seconds ==
            2 * cost$seconds,
          env$norm_exact_search_cost(107L, 8L, 21L, "Raw", 32L)$seconds ==
            cost$seconds)
# Large K alone is not sufficient, and a smaller K can still explode.
stopifnot(!env$norm_exact_search_cost(8L, 8L, 3L, "Raw", 1L)$explosion,
          env$norm_exact_search_cost(512L, 5L, 21L, "Raw", 16L)$explosion,
          env$norm_exact_search_cost(107L, 5L, 21000L, "Raw", 16L)$explosion)
stopifnot(length(capture.output(env$norm_report_exact_cost(cost, FALSE))) == 0L)
text <- paste(capture.output(env$norm_report_exact_cost(cost, TRUE)), collapse = "\n")
stopifnot(grepl("24,370,067,800 combinations", text, fixed = TRUE),
          grepl("6.8 minutes", text, fixed = TRUE),
          grepl("WARNING: combinatorial explosion", text, fixed = TRUE),
          grepl("not a 95% interval", text, fixed = TRUE))
stopifnot(env$norm_format_search_time(0.1) == "less than 1 second",
          grepl("hours", env$norm_format_search_time(7200)),
          grepl("days", env$norm_format_search_time(172800)),
          grepl("years", env$norm_format_search_time(1e9)),
          grepl("beyond", env$norm_format_search_time(Inf)))

# Stop at entry to preprocessing: an explosive search must log BEFORE any work.
env$norm_prepare_counts <- function(...) stop("PREPROCESSING_SENTINEL")
df <- data.frame(Cohort = rep(paste0("C", 1:20), each = 2))
ranked <- data.frame(Subtype = seq_len(107))
run <- function(k, verbose) {
  captured_error <- NULL
  output <- capture.output(tryCatch(
    env$subtype_norm_exact(df, ranked, k, "Raw", 0, verbose,
                          num_cores = 16L),
    error = function(e) captured_error <<- conditionMessage(e)))
  list(output = paste(output, collapse = "\n"), error = captured_error)
}
large <- run(8L, TRUE)
stopifnot(large$error == "PREPROCESSING_SENTINEL",
          grepl("WARNING: combinatorial explosion", large$output, fixed = TRUE))
quiet <- run(8L, FALSE)
stopifnot(quiet$error == "PREPROCESSING_SENTINEL", quiet$output == "")
unsupported <- run(9L, TRUE)
stopifnot(grepl("supports 2 <= K", unsupported$error, fixed = TRUE),
          grepl("301,579,589,025 combinations", unsupported$output, fixed = TRUE),
          grepl("unsupported extrapolation", unsupported$output, fixed = TRUE))

# Compare real small-search results with the previous installed adapter.
rm("norm_prepare_counts", envir = env)
fixture <- expand.grid(Cohort = c("A", "B"), Subtype = 1:6, response = 0:1)
fixture$SampleIDs <- paste0("synthetic_", seq_len(nrow(fixture)))
fixture$tumor_type <- "Lung"
ranked <- get("norm_rank_subtypes", asNamespace("CCS"))(fixture, "RR")
args <- list(df = fixture, ranked = ranked, k = 4L, method = "Raw",
             started = proc.time()[["elapsed"]], verbose = FALSE, num_cores = 1L)
old <- do.call(get("subtype_norm_exact", asNamespace("CCS")), args)
new <- do.call(env$subtype_norm_exact, args)
stopifnot(identical(old$Raw, new$Raw), identical(old$Summary, new$Summary),
          new$SearchDiagnostics$n_candidates_evaluated == choose(5, 3),
          is.finite(new$SearchDiagnostics$estimated_search_seconds))
cat("Search-cost checks passed: calibration, scaling, warning timing, quiet mode,",
    "unsupported K and unchanged small-search results.\n")
