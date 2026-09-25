source(file.path("targets", "statistical_inference.R"), local = TRUE)

values <- c(4, 1, 1, 3, 4, 2)
stopifnot(identical(.asi_weighted_rank(values, rep(1, length(values))),
  as.numeric(rank(values))))

set.seed(40L)
direct <- matrix(rnorm(8L * 4L * 5L), ncol = 5L)
d1 <- cbind(direct[, 1:3, drop = FALSE] + matrix(rnorm(32L * 3L, sd = 0.2),
  ncol = 3L), direct[, 4:5, drop = FALSE])
cohort <- rep(sprintf("cohort_%02d", seq_len(8L)), each = 4L)
first <- sample.int(nrow(direct), 120L, replace = TRUE)
second <- sample.int(nrow(direct), 120L, replace = TRUE)
second[second == first] <- second[second == first] %% nrow(direct) + 1L
distance <- function(x) sqrt(rowSums((x[first, , drop = FALSE] -
  x[second, , drop = FALSE])^2))
run <- function() .asi_geometry_bootstrap(direct, d1, cohort, first, second,
  distance(direct), distance(d1), k = 3L, n_boot = 100L, seed = 40L)
actual <- run()
stopifnot(identical(actual, run()), nrow(actual) == 3L,
  all(actual$valid_resamples == 100L), all(is.finite(actual$ci_low)),
  all(is.finite(actual$ci_high)), all(actual$ci_low <= actual$ci_high),
  all(is.na(actual$p_value)), all(actual$status == "estimable"))
