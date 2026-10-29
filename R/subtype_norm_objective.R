# Shared candidate x cohort scoring for both normCCS search algorithms.
norm_validate_objective <- function(objective) {
  choices <- c("median", "cohort_median", "mean", "lower_tail", "mean_lower_tail")
  if (!is.character(objective) || length(objective) != 1L ||
      is.na(objective) || !objective %in% choices) {
    stop("norm_objective must be one of: ", paste(choices, collapse = ", "), ".")
  }
  objective
}

norm_objective_score <- function(auc, cohorts, objective = "median") {
  # Test
  if (FALSE) {
    auc <- matrix(c(0.9, 0.5, 0.7, 0.8, 0.9, 1), nrow = 1L)
    norm_objective_score(auc, c("All", LETTERS[1:5]), "lower_tail")
  }
  objective <- norm_validate_objective(objective)
  if (!is.matrix(auc) || !is.numeric(auc) || ncol(auc) != length(cohorts) ||
      !is.character(cohorts) || anyNA(cohorts) || any(!nzchar(cohorts)) ||
      anyDuplicated(cohorts)) {
    stop("auc must be a numeric candidate-by-cohort matrix with unique, nonempty cohort identifiers.")
  }
  if (objective == "median") {
    return(apply(auc, 1L, stats::median, na.rm = TRUE))
  }
  real <- auc[, cohorts != "All", drop = FALSE]
  if (!ncol(real)) stop("norm_objective requires at least one real cohort (excluding 'All').")
  if (any(!is.finite(real))) {
    bad <- which(!is.finite(real), arr.ind = TRUE)[1L, ]
    stop("Non-finite AUC at candidate ", bad[1L], ", cohort '",
         cohorts[cohorts != "All"][bad[2L]], "'.")
  }
  if (objective == "mean") return(rowMeans(real))
  if (objective == "cohort_median") return(apply(real, 1L, stats::median))
  mass <- 0.25 * ncol(real)
  whole <- floor(mass)
  fraction <- mass - whole
  tail <- apply(real, 1L, function(x) {
    x <- sort(x)
    (sum(x[seq_len(whole)]) + if (fraction > 0) fraction * x[whole + 1L] else 0) / mass
  })
  if (objective == "lower_tail") return(tail)
  0.5 * rowMeans(real) + 0.5 * tail
}

norm_objective_settings <- function(objective) {
  settings <- list(norm_objective = objective)
  if (objective %in% c("lower_tail", "mean_lower_tail")) settings$tail_fraction <- 0.25
  if (objective == "mean_lower_tail") settings$mean_weight <- 0.5
  settings
}

# Preserve legacy columns; never relabel a new score as medianROCAUC.
norm_objective_summary <- function(raw, objective) {
  infos <- sort(unique(as.character(raw$Info)))
  cohorts <- sort(unique(as.character(raw$Cohort)))
  key <- paste(match(raw$Info, infos), match(raw$Cohort, cohorts), sep = ":")
  if (anyDuplicated(key) || nrow(raw) != length(infos) * length(cohorts)) {
    stop("Each candidate must have exactly one AUC for every cohort.")
  }
  auc <- matrix(NA_real_, length(infos), length(cohorts))
  auc[cbind(match(raw$Info, infos), match(raw$Cohort, cohorts))] <- raw$ROCAUC
  summary <- data.frame(Info = infos,
                        medianROCAUC = norm_objective_score(auc, cohorts, "median"),
                        stringsAsFactors = FALSE)
  summary$objective <- objective
  summary$objectiveScore <- norm_objective_score(auc, cohorts, objective)
  summary$cohortMedianROCAUC <- norm_objective_score(auc, cohorts, "cohort_median")
  summary$cohortMeanROCAUC <- norm_objective_score(auc, cohorts, "mean")
  summary$cohortLowerTailROCAUC <- norm_objective_score(auc, cohorts, "lower_tail")
  summary <- summary[order(-summary$objectiveScore, -summary$cohortMeanROCAUC, summary$Info), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}
