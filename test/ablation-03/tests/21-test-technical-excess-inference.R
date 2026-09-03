# TDD RED test: technical-neighbor excess cohort-level inference.

source(file.path("test", "ablation-03", "02-ablation03-experiment_functions.R"))
if (!exists(".ae_technical_inference", mode = "function")) {
  stop("Expected technical inference helper is missing.", call. = FALSE)
}

sample_ids <- paste0("s", 1:12)
cohorts <- rep(paste0("cohort", 1:4), each = 3L)
make_arm <- function(representation, shift) {
  data.frame(
    sample_id = sample_ids,
    cohort = cohorts,
    label = rep(c("A", "B"), length.out = length(sample_ids)),
    k = 15L,
    assay_type_match_excess = rep(c(0.10, 0.08, 0.12), 4L) + shift,
    platform_id_match_excess = rep(c(0.12, 0.10, 0.14), 4L) + shift,
    source_system_match_excess = rep(c(0.06, 0.04, 0.08), 4L) + shift,
    representation = representation,
    stringsAsFactors = FALSE
  )
}
retrieval <- list(
  per_sample = rbind(
    make_arm("Direct-GSClassifier", 0),
    make_arm("Cohort-d1", -0.02)
  )
)

result <- .ae_technical_inference(
  retrieval,
  k = 15L,
  technical_columns = c("assay_type", "platform_id", "source_system"),
  n_boot = 200L,
  seed = 20260903L
)

stopifnot(nrow(result) == 3L)
stopifnot(all(c(
  "technical_factor", "estimate", "ci_low", "ci_high", "p_value",
  "p_value_adj", "n_cohort", "status"
) %in% names(result)))
stopifnot(all(result$n_cohort == 4L))
stopifnot(all(result$estimate < 0))
stopifnot(all(result$ci_low <= result$estimate))
stopifnot(all(result$estimate <= result$ci_high))
stopifnot(all(is.finite(result$p_value)))
stopifnot(all(is.finite(result$p_value_adj)))
stopifnot(all(result$multiplicity_method == "holm"))
cat("technical excess inference test passed\n")
