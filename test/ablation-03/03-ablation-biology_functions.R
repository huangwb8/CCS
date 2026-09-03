# Statistical helpers for biological-anchor readout.
# The inferential unit is the external query cohort, not an individual
# neighbour pair. This prevents the 15 neighbours of one query from being
# treated as independent observations.

.biology_paired_contrast <- function(
    per_query,
    n_boot = 2000L,
    seed = 20260830L,
    min_cohorts = 3L) {
  required <- c(
    "anchor", "representation", "query_sample", "query_cohort", "utility"
  )
  missing_columns <- setdiff(required, names(per_query))
  if (length(missing_columns)) {
    stop(
      "biology inference: missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (n_boot < 100L) {
    stop("biology inference: n_boot must be at least 100.", call. = FALSE)
  }

  per_query <- per_query[is.finite(per_query$utility), required, drop = FALSE]
  if (!nrow(per_query)) {
    return(data.frame())
  }

  # Guard against duplicated rows from upstream joins before cohort means.
  per_query <- stats::aggregate(
    utility ~ anchor + representation + query_sample + query_cohort,
    data = per_query,
    FUN = mean
  )
  cohort_means <- stats::aggregate(
    utility ~ anchor + query_cohort + representation,
    data = per_query,
    FUN = mean
  )
  direct <- cohort_means[
    cohort_means$representation == "Direct-GSClassifier",
    c("anchor", "query_cohort", "utility"),
    drop = FALSE
  ]
  d1 <- cohort_means[
    cohort_means$representation == "Cohort-d1",
    c("anchor", "query_cohort", "utility"),
    drop = FALSE
  ]
  names(direct)[3L] <- "utility_direct"
  names(d1)[3L] <- "utility_d1"
  paired <- merge(direct, d1, by = c("anchor", "query_cohort"))
  if (!nrow(paired)) {
    return(data.frame())
  }
  paired$delta_d1_minus_direct <- paired$utility_d1 - paired$utility_direct

  split_paired <- split(paired, paired$anchor, drop = TRUE)
  out <- lapply(seq_along(split_paired), function(i) {
    part <- split_paired[[i]]
    delta <- part$delta_d1_minus_direct
    n_cohorts <- length(delta)
    estimate <- mean(delta)
    set.seed(seed + i - 1L)
    bootstrap <- replicate(
      n_boot,
      mean(sample(delta, n_cohorts, replace = TRUE))
    )
    if (n_cohorts <= 20L) {
      # Exact paired sign-flip test. With 17 external cohorts this is cheap
      # and avoids a normal approximation with a small number of clusters.
      masks <- 0:(2^n_cohorts - 1L)
      null_distribution <- vapply(masks, function(mask) {
        bits <- as.integer(intToBits(mask))[seq_len(n_cohorts)]
        signs <- ifelse(bits == 1L, 1, -1)
        mean(delta * signs)
      }, numeric(1L))
      p_value <- (
        sum(abs(null_distribution) >= abs(estimate)) + 1
      ) / (length(null_distribution) + 1)
      p_method <- "exact_paired_sign_flip"
    } else {
      signs <- replicate(
        n_boot,
        sample(c(-1, 1), n_cohorts, replace = TRUE)
      )
      null_distribution <- colMeans(signs * delta)
      p_value <- (sum(abs(null_distribution) >= abs(estimate)) + 1) /
        (length(null_distribution) + 1)
      p_method <- "monte_carlo_paired_sign_flip"
    }
    data.frame(
      anchor = part$anchor[1L],
      estimate = estimate,
      ci_low = unname(stats::quantile(bootstrap, 0.025)),
      ci_high = unname(stats::quantile(bootstrap, 0.975)),
      p_value = p_value,
      p_method = p_method,
      cohort_count = n_cohorts,
      query_count = length(intersect(
        per_query$query_sample[
          per_query$anchor == part$anchor[1L] &
            per_query$representation == "Direct-GSClassifier" &
            per_query$query_cohort %in% part$query_cohort
        ],
        per_query$query_sample[
          per_query$anchor == part$anchor[1L] &
            per_query$representation == "Cohort-d1" &
            per_query$query_cohort %in% part$query_cohort
        ]
      )),
      n_boot = n_boot,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, out)
  result$p_value_adj <- stats::p.adjust(result$p_value, method = "BH")
  result$inference_status <- ifelse(
    result$cohort_count >= min_cohorts,
    "estimable",
    "not_estimable"
  )
  result
}
