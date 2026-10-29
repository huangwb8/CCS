test_that("objective definitions match independent discrete probability weights", {
  score <- CCS:::norm_objective_score
  for (nc in c(1L, 2L, 3L, 4L, 5L, 8L)) {
    real <- seq(0.5, 0.95, length.out = nc)
    # Integrate the empirical quantile function over [0, 0.25].
    weight <- pmax(0, pmin(seq_len(nc) / nc, 0.25) - (seq_len(nc) - 1) / nc) / 0.25
    lower <- sum(sort(real) * weight)
    auc <- matrix(c(real, 0.99), nrow = 1L)
    cohorts <- c(paste0("C", seq_len(nc)), "All")
    expected <- c(median = stats::median(c(real, 0.99)),
                  cohort_median = stats::median(real), mean = mean(real),
                  lower_tail = lower, mean_lower_tail = (mean(real) + lower) / 2)
    for (objective in names(expected)) {
      expect_equal(unname(score(auc, cohorts, objective)), unname(expected[objective]), tolerance = 1e-12)
      expect_equal(score(rbind(auc, auc), cohorts, objective), rep(unname(expected[objective]), 2L))
      expect_equal(score(auc[, rev(seq_along(cohorts)), drop = FALSE], rev(cohorts), objective),
                   score(auc, cohorts, objective))
      if (objective != "median") {
        modified <- auc
        modified[, nc + 1L] <- NA_real_
        expect_equal(score(modified, cohorts, objective), score(auc, cohorts, objective))
        for (j in seq_len(nc)) {
          improved <- auc
          improved[, j] <- improved[, j] + 0.01
          expect_gte(score(improved, cohorts, objective), score(auc, cohorts, objective) - 1e-12)
        }
        expect_equal(score(matrix(0.7, 2L, nc + 1L), cohorts, objective), rep(0.7, 2L))
      }
    }
  }
  # Five cohorts require a fractional contribution from the second worst.
  expect_equal(score(matrix(c(0.9, 0.5, 0.6, 0.8, 0.9, 1), 1L),
                     c("All", LETTERS[1:5]), "lower_tail"), 0.52)
})

test_that("invalid objectives and cohort matrices fail explicitly", {
  for (bad in list(NULL, NA_character_, "", "mea", "unknown", 1, c("mean", "median"))) {
    expect_error(CCS:::norm_validate_objective(bad), "norm_objective")
    expect_error(CCS::subtypePerformance(NULL, norm_objective = bad), "norm_objective")
  }
  score <- CCS:::norm_objective_score
  expect_error(score(matrix(0.5, 1, 2), c("A", "A"), "mean"), "unique")
  expect_error(score(matrix(0.5, 1, 1), "All", "mean"), "real cohort")
  expect_error(score(matrix(c(0.8, NA), 1), c("All", "A"), "mean"), "candidate 1, cohort 'A'")
  expect_error(score(matrix(c(0.8, Inf), 1), c("All", "A"), "lower_tail"), "Non-finite")
  expect_equal(score(matrix(c(0.8, NA), 1), c("All", "A"), "median"), 0.8)
})

test_that("cohort objectives change ranking without overwriting the legacy median", {
  raw <- data.frame(Info = rep(c("A", "B"), each = 6L),
                    Cohort = rep(c("All", LETTERS[1:5]), 2L),
                    ROCAUC = c(0.8, 0.8, 0.8, 0.8, 0.51, 0.5, rep(0.75, 6L)))
  result <- CCS:::norm_objective_summary(raw, "mean_lower_tail")
  expect_identical(result$Info, c("B", "A"))
  expect_equal(result$medianROCAUC, c(0.75, 0.8))
  tied <- raw
  tied$ROCAUC <- 0.7
  expect_identical(CCS:::norm_objective_summary(tied, "mean")$Info, c("A", "B"))
  expect_error(CCS:::norm_objective_summary(raw[-1, ], "mean"), "exactly one")
})
