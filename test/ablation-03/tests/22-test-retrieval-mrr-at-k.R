# Regression test: MRR@k must not credit hits ranked below k.
source(file.path("R", "ablation.R"))

reference <- matrix(seq_len(6L), ncol = 1L)
query <- matrix(0, nrow = 1L, ncol = 1L)
reference_metadata <- data.frame(
  sample_id = paste0("r", seq_len(6L)),
  cohort = rep("reference", 6L),
  label = c("A", "A", "A", "B", "A", "A"),
  stringsAsFactors = FALSE
)
query_metadata <- data.frame(
  sample_id = "q1",
  cohort = "query",
  label = "B",
  stringsAsFactors = FALSE
)

result <- .ablation_query_reference_retrieval(
  reference = reference,
  query = query,
  reference_metadata = reference_metadata,
  query_metadata = query_metadata,
  label_column = "label",
  k = c(2L, 5L),
  search = "exact"
)

mrr <- result$per_sample$mrr[order(result$per_sample$k)]
stopifnot(identical(as.integer(result$per_sample$k[order(result$per_sample$k)]), c(2L, 5L)))
stopifnot(mrr[1L] == 0)
stopifnot(mrr[2L] == 0.25)
stopifnot(result$summary$mrr[result$summary$k == 2L] == 0)
stopifnot(result$summary$mrr[result$summary$k == 5L] == 0.25)
cat("retrieval MRR@k test passed\n")
