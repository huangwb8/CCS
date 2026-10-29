test_that("interval DP and legacy scoring agree on a small synthetic problem", {
  grid <- expand.grid(Cohort = c("A", "B", "C"), Subtype = letters[1:6],
                      response = 0:1, replicate = 1:3, stringsAsFactors = FALSE)
  grid$SampleIDs <- paste0("synthetic_", seq_len(nrow(grid)))
  grid$tumor_type <- "Lung"
  ranked <- data.frame(Subtype = letters[1:6])
  model <- CCS:::norm_prepare_counts(grid, ranked)
  for (k in 2:6) {
    cuts <- t(utils::combn(1:5, k - 1L))
    auc <- CCS:::norm_auc_batch(cuts, model, "Raw")
    expect_equal(as.vector(auc), rep(0.5, length(auc)))
    fit <- CCS:::norm_solve_dp(model$loss[1, , ], k)
    expect_equal(fit$loss, 0, tolerance = 1e-12)
    result <- CCS:::subtypeNorm(grid, k, "Raw", "RR", 1, FALSE, "dp")
    expect_equal(result$Summary$medianROCAUC, rep(0.5, nrow(result$Summary)))
    expect_false(result$SearchDiagnostics$exact_global_optimum)
  }
})

test_that("input failures are caught before allocating partitions or starting workers", {
  expect_error(CCS:::subtypeNorm(data.frame(), 5), "nonempty")
  df <- data.frame(Cohort = c("A", "A"), Subtype = c("a", "b"),
                   response = c(0, 1), SampleIDs = c("s1", "s2"), tumor_type = "Lung")
  expect_error(CCS:::subtypeNorm(df, 3), "observed")
  expect_error(CCS:::subtypeNorm(df, 2.5), "integer")
  expect_error(CCS:::subtypeNorm(df, 2, numCores = 0), "numCores")
  expect_error(CCS:::subtypeNorm(df, 2, norm_search_method = "unknown"), "arg")
})

test_that("OptimizeRank preserves pROC last-bit ordering and factor response levels", {
  p <- c(12L, 15L, 2L, 0L)
  n <- c(20L, 14L, 14L, 10L)
  df <- do.call(rbind, lapply(1:4, function(i) {
    data.frame(Cohort = "A", Subtype = letters[i],
               response = c(rep(1L, p[i]), rep(0L, n[i])), tumor_type = "Lung")
  }))
  df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
  for (reverse in c(FALSE, TRUE)) {
    if (reverse) df$response <- factor(df$response, levels = c(1, 0))
    model <- CCS:::norm_prepare_counts(df, data.frame(Subtype = letters[1:4]))
    predictor <- 4L - match(df$Subtype, letters[1:4])
    for (method in c("Raw", "OptimizeRank")) {
      fast <- CCS:::norm_auc_batch(matrix(1:3, nrow = 1), model, method)
      exact <- CCS:::get_perform_markers(df$response, predictor, 4L, method)$ROCAUC
      expect_equal(as.vector(fast), rep(as.numeric(exact), 2L), tolerance = 1e-10)
    }
  }
})

test_that("nonmonotone interval costs and DP match complete enumeration", {
  p <- matrix(c(5, 0, 9, 2, 6, 1, 1, 8, 2, 0, 4, 7, 9, 1, 0, 8, 2, 3), 3, byrow = TRUE)
  n <- matrix(c(1, 5, 0, 8, 2, 9, 7, 1, 6, 9, 0, 3, 0, 8, 9, 1, 5, 4), 3, byrow = TRUE)
  df <- do.call(rbind, lapply(1:3, function(c) {
    do.call(rbind, lapply(1:6, function(s) {
      data.frame(Cohort = LETTERS[c], Subtype = letters[s],
                 response = c(rep(1L, p[c, s]), rep(0L, n[c, s])), tumor_type = "Lung")
    }))
  }))
  df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
  model <- CCS:::norm_prepare_counts(df, data.frame(Subtype = letters[1:6]))
  expect_true(any(model$loss < 0))
  cost <- apply(model$loss * c(0, 0.5, 0.3, 0.2), c(2, 3), sum)
  for (k in 2:6) {
    cuts <- t(utils::combn(1:5, k - 1L))
    values <- apply(cuts, 1, function(x) sum(cost[cbind(c(1, x + 1L), c(x, 6))]))
    expect_equal(CCS:::norm_solve_dp(cost, k)$loss, min(values), tolerance = 1e-12)
    fast <- CCS:::norm_auc_batch(cuts, model, "Raw")
    for (i in seq_len(nrow(cuts))) for (c in seq_along(model$rows)) {
      ids <- model$rows[[c]]
      score <- k - 1L - findInterval(match(df$Subtype[ids], letters[1:6]) - 1L, cuts[i, ])
      exact <- pROC::roc(df$response[ids], score, quiet = TRUE)$auc
      expect_equal(unname(fast[i, c]), as.numeric(exact), tolerance = 1e-10)
    }
  }
})
