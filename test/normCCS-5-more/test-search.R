# Formal component regression entry: Rscript test/normCCS-5-more/test-search.R
source("test/normCCS-5-more/00.Environment.R", encoding = "UTF-8")
devtools::load_all(quiet = TRUE)
legacy_env <- new.env(parent = asNamespace("CCS"))
sys.source("test/normCCS-5-more/reference/subtype-norm-0.8.4.R", legacy_env)
norm <- getFromNamespace("subtypeNorm", "CCS")
prepare <- getFromNamespace("norm_prepare_counts", "CCS")
batch <- getFromNamespace("norm_auc_batch", "CCS")
solve <- getFromNamespace("norm_solve_dp", "CCS")
info <- getFromNamespace("norm_info", "CCS")
marker_eval <- getFromNamespace("get_perform_markers", "CCS")
test_args <- commandArgs(trailingOnly = TRUE)
section <- if (length(test_args)) test_args[1L] else "all"
legacy_k <- if (length(test_args) > 1L) as.integer(test_args[2L]) else 2:4

if (section %in% c("all", "legacy")) {
  df <- make_norm_fixture(6L, 3L)
  for (k in legacy_k) for (method in c("Raw", "OptimizeRank")) {
    for (rank in c("RR", "NRR", "RR+NRR", "NRR+RR")) {
      cat("Legacy regression", k, method, rank, "\n")
      args <- list(df = df, n_norm_subtype = k, norm_roc_method = method,
                   norm_rank_marker = rank, numCores = 1L, verbose = FALSE)
      old <- suppressWarnings(do.call(legacy_env$subtypeNorm, args))
      default <- suppressWarnings(do.call(norm, args))
      explicit <- suppressWarnings(do.call(norm, c(args, list(norm_search_method = "bruteforce"))))
      testthat::expect_equal(default, old, tolerance = 1e-12)
      testthat::expect_equal(explicit, old, tolerance = 1e-12)
    }
  }
}

if (section %in% c("all", "kernels")) testthat::test_that("all partitions agree with pROC, including OptimizeRank and automatic direction", {
  for (seed in c(1:5, 18L, 1001L, 1018L)) {
    base_seed <- seed %% 1000L
    df <- make_norm_fixture(7L, 3L, base_seed, sparse = TRUE)
    ranked <- if (base_seed == 18L) data.frame(Subtype = sort(unique(df$Subtype))) else
      reference_rank(df, c("RR", "NRR")[(base_seed - 1L) %% 2L + 1L])
    if (seed > 1000L) df$response <- factor(df$response, levels = c(1, 0))
    model <- prepare(df, ranked)
    for (k in 2:5) {
      cuts <- enumerate_cuts(model$m, k)
      for (method in c("Raw", "OptimizeRank")) {
        auc <- batch(cuts, model, method)
        for (i in seq_len(nrow(cuts))) {
          group <- k - 1L - findInterval(match(df$Subtype, ranked$Subtype) - 1L, cuts[i, ])
          for (c in seq_along(model$rows)) {
            ids <- model$rows[[c]]
            expected <- suppressWarnings(marker_eval(df$response[ids], group[ids], k, method)$ROCAUC)
            testthat::expect_equal(unname(auc[i, c]), as.numeric(expected), tolerance = 1e-10,
                                   info = paste(seed, k, method, i, c))
          }
          testthat::expect_match(info(cuts[i, ], ranked), "^0=")
        }
      }
      for (weight in list(rep(1 / 4, 4), c(0, 0.1, 0.3, 0.6), c(1, -1, 0, 0))) {
        cost <- apply(model$loss * weight, c(2, 3), sum)
        dp <- solve(cost, k)
        brute_loss <- apply(cuts, 1, function(x) sum(cost[cbind(c(1, x + 1L), c(x, model$m))]))
        testthat::expect_equal(dp$loss, min(brute_loss), tolerance = 1e-12)
      }
    }
  }
})

if (section %in% c("all", "state")) testthat::test_that("DP is deterministic, safe, bounded, and keeps caller state", {
  df <- make_norm_fixture(12L, 4L, 99, sparse = TRUE)
  set.seed(42)
  rng <- .Random.seed
  old_options <- options()
  backend <- foreach::getDoParName()
  old_wd <- getwd()
  a <- norm(df, 5, "Raw", "NRR", 1, FALSE, "dp")
  b <- norm(df, 5, "Raw", "NRR", 2, FALSE, "dp")
  testthat::expect_equal(a$Raw, b$Raw)
  testthat::expect_equal(a$Summary, b$Summary)
  testthat::expect_identical(.Random.seed, rng)
  testthat::expect_identical(options(), old_options)
  testthat::expect_identical(foreach::getDoParName(), backend)
  testthat::expect_identical(getwd(), old_wd)
  testthat::expect_false(a$SearchDiagnostics$exact_global_optimum)
  testthat::expect_equal(a$SearchDiagnostics$workers, 1L)
  testthat::expect_gte(a$Summary$medianROCAUC[1], a$SearchDiagnostics$best_initial_median_auc - 1e-12)
  full <- norm(df, 12, "Raw", "RR", 1, FALSE, "dp")
  testthat::expect_equal(nrow(full$Summary), 1L)
  for (k in list(1, 2.5, NA_real_, 13)) testthat::expect_error(norm(df, k, verbose = FALSE))
  testthat::expect_error(norm(df, 5, numCores = 0, verbose = FALSE), "numCores")
  df$response[1] <- NA
  testthat::expect_error(norm(df, 5, verbose = FALSE), "complete")
})

cat("All normCCS component tests passed.\n")
