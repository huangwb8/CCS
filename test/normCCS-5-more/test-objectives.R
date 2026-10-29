# Rscript test/normCCS-5-more/test-objectives.R [baseline|objectives|all]
source("test/normCCS-5-more/00.Environment.R", encoding = "UTF-8")
devtools::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
section <- if (length(args)) args[1L] else "all"
if (!section %in% c("all", "baseline", "objectives")) stop("Unknown objective test section.")
norm <- getFromNamespace("subtypeNorm", "CCS")
objectives <- c("median", "cohort_median", "mean", "lower_tail", "mean_lower_tail")

# Independent quantile integration, not the package scorer.
reference_objective <- function(x, cohorts, objective) {
  if (objective == "median") return(stats::median(x, na.rm = TRUE))
  x <- x[cohorts != "All"]
  nc <- length(x)
  weights <- pmax(0, pmin(seq_len(nc) / nc, 0.25) - (seq_len(nc) - 1) / nc) / 0.25
  lower <- sum(sort(x) * weights)
  switch(objective, cohort_median = stats::median(x), mean = mean(x),
         lower_tail = lower, mean_lower_tail = (mean(x) + lower) / 2)
}

if (section %in% c("all", "baseline")) {
  frozen <- readRDS("test/normCCS-5-more/reference/baseline-0.8.5.rds")
  strip_runtime <- function(x) {
    if (!is.null(x$SearchDiagnostics)) x$SearchDiagnostics$elapsed_seconds <- NULL
    x
  }
  for (i in seq_len(nrow(frozen$grid))) {
    x <- frozen$grid[i, ]
    cat("0.8.5 regression", i, "/", nrow(frozen$grid), "\n")
    arguments <- list(frozen$df, x$k, x$method, x$rank, 1L, FALSE, x$search)
    current <- suppressWarnings(do.call(norm, arguments))
    explicit <- suppressWarnings(do.call(norm, c(arguments, list(norm_objective = "median"))))
    testthat::expect_equal(strip_runtime(current), strip_runtime(frozen$results[[i]]), tolerance = 1e-12)
    testthat::expect_equal(strip_runtime(explicit), strip_runtime(current), tolerance = 1e-12)
  }
}

if (section %in% c("all", "objectives")) {
  df <- make_norm_fixture(7L, 5L, 18L, sparse = TRUE)
  for (method in c("Raw", "OptimizeRank")) {
    for (objective in objectives) {
      cat("Objective integration", method, objective, "\n")
      exhaustive <- suppressWarnings(norm(df, 4L, method, "NRR", 1L, FALSE, "bruteforce", objective))
      dp <- norm(df, 4L, method, "NRR", 1L, FALSE, "dp", objective)
      expected <- vapply(split(exhaustive$Raw, exhaustive$Raw$Info), function(x) {
        reference_objective(x$ROCAUC, x$Cohort, objective)
      }, numeric(1))
      get_best <- function(x) if (objective == "median") x$Summary$medianROCAUC[1L] else x$Summary$objectiveScore[1L]
      testthat::expect_equal(get_best(exhaustive), max(expected), tolerance = 1e-12)
      testthat::expect_lte(get_best(dp), max(expected) + 1e-10)
      if (objective != "median") {
        expected_summary <- exhaustive$Summary
        testthat::expect_equal(expected_summary$objectiveScore,
                               unname(expected[expected_summary$Info]), tolerance = 1e-12)
        expected_mean <- vapply(split(exhaustive$Raw, exhaustive$Raw$Info), function(x) mean(x$ROCAUC[x$Cohort != "All"]), numeric(1))
        expected_order <- order(-expected, -expected_mean, names(expected))
        testthat::expect_identical(expected_summary$Info, names(expected)[expected_order])
        d <- dp$SearchDiagnostics
        testthat::expect_identical(d$objective, objective)
        testthat::expect_identical(d$objective_scope, "cohorts_only")
        testthat::expect_equal(d$n_real_cohorts, 5L)
        testthat::expect_false(d$exact_global_optimum)
        testthat::expect_true(d$dp_subproblems_exact)
        testthat::expect_null(d$best_initial_median_auc)
        testthat::expect_equal(d$best_final_objective_score, get_best(dp))
        testthat::expect_gte(get_best(dp), d$best_initial_objective_score - 1e-12)
        for (refinement in d$refinements) {
          testthat::expect_gte(refinement$final_objective_score, refinement$initial_objective_score - 1e-12)
          testthat::expect_null(refinement$initial_auc)
        }
        first <- dp$Raw[dp$Raw$Info == dp$Summary$Info[1L], ]
        testthat::expect_equal(get_best(dp), reference_objective(first$ROCAUC, first$Cohort, objective), tolerance = 1e-10)
        testthat::expect_equal(dp$Summary$medianROCAUC[1L], median(first$ROCAUC))
      }
    }
  }
  # No RNG/options/backend/working-directory changes on any new DP objective.
  for (objective in objectives[-1L]) {
    set.seed(42)
    before <- list(rng = .Random.seed, options = options(), wd = getwd(), backend = foreach::getDoParName())
    a <- norm(df, 5L, "Raw", "RR", 1L, FALSE, "dp", objective)
    b <- norm(df, 5L, "Raw", "RR", 2L, FALSE, "dp", objective)
    after <- list(rng = .Random.seed, options = options(), wd = getwd(), backend = foreach::getDoParName())
    testthat::expect_identical(before, after)
    testthat::expect_equal(a$Raw, b$Raw)
    testthat::expect_equal(a$Summary, b$Summary)
    full <- norm(df, 7L, "Raw", "RR", 1L, FALSE, "dp", objective)
    testthat::expect_equal(nrow(full$Summary), 1L)
  }
}
cat("Objective tests passed:", section, "\n")
