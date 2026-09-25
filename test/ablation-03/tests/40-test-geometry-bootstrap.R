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

# The kNN interval must resample query cohorts against the complete frozen
# neighbor graph. Rebuilding a smaller atlas inside each draw changes the metric.
knn <- getFromNamespace(".ablation_knn", "CCS")
direct_neighbors <- knn(direct, 3L)
d1_neighbors <- knn(d1, 3L)
agreement <- vapply(seq_len(nrow(direct)), function(i) {
  length(intersect(direct_neighbors[i, ], d1_neighbors[i, ])) /
    length(union(direct_neighbors[i, ], d1_neighbors[i, ]))
}, numeric(1))
expected <- vapply(seq_len(100L), function(i) {
  set.seed(40L + i - 1L)
  selected <- sample(unique(cohort), length(unique(cohort)), replace = TRUE)
  weight <- tabulate(match(selected, unique(cohort)),
    nbins = length(unique(cohort)))
  stats::weighted.mean(agreement, weight[match(cohort, unique(cohort))])
}, numeric(1))
knn_row <- actual[actual$endpoint == "knn_jaccard", ]
stopifnot(identical(knn_row$method, "query_cohort_percentile_bootstrap"),
  isTRUE(all.equal(unname(c(knn_row$ci_low, knn_row$ci_high)),
    unname(stats::quantile(expected, c(0.025, 0.975))), tolerance = 1e-12)))
