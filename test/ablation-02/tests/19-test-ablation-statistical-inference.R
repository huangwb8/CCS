#!/usr/bin/env Rscript

# Purpose: Contract tests for cohort-level bootstrap and paired sign-flip inference.
# Input: Small deterministic paired synthetic data; no genomic data are read.
# Output: stopifnot failures make statistical direction and status regressions clear.

source(file.path("test", "ablation-02", "02-ablation-experiment_functions.R"))

positive <- data.frame(
  cohort = rep(paste0("C", 1:4), each = 2),
  delta = rep(c(0.20, 0.30, 0.10, 0.25), each = 2),
  stringsAsFactors = FALSE
)
fit <- .ae_paired_inference(
  positive, "delta", n_boot = 200L, seed = 11L,
  unit = "cohort", multiplicity_method = "holm"
)
stopifnot(fit$status == "complete", fit$n_cohort == 4L, fit$n_sample == 8L)
stopifnot(fit$estimate > 0, fit$ci_low > 0, fit$ci_high > fit$ci_low)
stopifnot(is.finite(fit$p_value), fit$method == "cluster_bootstrap_sign_flip")

zero <- positive
zero$delta <- 0
zero_fit <- .ae_paired_inference(zero, "delta", n_boot = 200L, seed = 11L)
stopifnot(zero_fit$status == "complete", zero_fit$estimate == 0, zero_fit$p_value == 1)

# Duplicated samples inside one cohort must not inflate the independent cluster count.
clustered <- rbind(
  data.frame(cohort = "A", delta = rep(1, 100)),
  data.frame(cohort = "B", delta = rep(-1, 2))
)
cluster_fit <- .ae_paired_inference(clustered, "delta", n_boot = 200L, seed = 12L)
stopifnot(cluster_fit$n_cohort == 2L, cluster_fit$n_sample == 102L)
stopifnot(abs(cluster_fit$estimate) < 1e-12)

missing_fit <- .ae_paired_inference(data.frame(cohort = "A"), "delta")
stopifnot(missing_fit$status == "not_estimable", missing_fit$reason == "missing_required_columns")

insufficient_fit <- .ae_paired_inference(
  data.frame(cohort = rep("A", 3), delta = 1:3),
  "delta", n_boot = 200L
)
stopifnot(insufficient_fit$status == "not_estimable")
stopifnot(insufficient_fit$reason == "fewer_than_two_independent_clusters")

adjusted <- .ae_adjust_inference(
  data.frame(p_value = c(0.01, 0.04, 0.50)), method = "holm"
)
stopifnot(isTRUE(all.equal(adjusted$p_value_adj, stats::p.adjust(c(0.01, 0.04, 0.50), "holm"))))

cat("19-test-ablation-statistical-inference: PASS\n")
