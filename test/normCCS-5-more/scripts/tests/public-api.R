# Rscript test/normCCS-5-more/scripts/tests/public-api.R
source("test/normCCS-5-more/00.Environment.R", encoding = "UTF-8")
devtools::load_all(quiet = TRUE)
grDevices::pdf(NULL)
df <- make_norm_fixture(8L, 4L, 901L, per_state = 24L)
# Existing public Raw ROC expects numeric state scores, as in CCS output.
df$Subtype <- as.integer(sub("S", "", df$Subtype))
fit <- suppressWarnings(CCS::subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 5L, norm_search_method = "dp", norm_roc_method = "Raw",
  numCores = 1L, verbose = FALSE))
stopifnot(identical(fit$Repeat$norm_search_method, "dp"),
          fit$Data$Normalization$SearchDiagnostics$K == 5L,
          is.finite(fit$Data$Normalization$Summary$medianROCAUC[1L]))
none <- suppressWarnings(CCS::subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type", n_norm_subtype = NULL,
  norm_search_method = "dp", norm_roc_method = "Raw", numCores = 1L, verbose = FALSE))
stopifnot(is.null(none$Data$Normalization), identical(none$Repeat$norm_search_method, "dp"))
stopifnot(is.null(fit$Repeat$norm_objective), is.null(none$Repeat$norm_objective))
for (objective in c("cohort_median", "mean", "lower_tail", "mean_lower_tail")) {
  for (search in c("dp", "bruteforce")) {
    result <- suppressWarnings(CCS::subtypePerformance(
      df, col_subtype = "Subtype", col_tissue = "tumor_type", n_norm_subtype = 5L,
      norm_search_method = search, norm_objective = objective, norm_roc_method = "Raw",
      numCores = 1L, verbose = FALSE))
    stopifnot(identical(result$Repeat$norm_objective, objective),
              identical(result$Data$Normalization$Summary$objective[1L], objective))
    if (objective %in% c("lower_tail", "mean_lower_tail")) stopifnot(result$Repeat$tail_fraction == 0.25)
    if (objective == "mean_lower_tail") stopifnot(result$Repeat$mean_weight == 0.5)
  }
}
none_new <- suppressWarnings(CCS::subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type", n_norm_subtype = NULL,
  norm_objective = "mean_lower_tail", norm_roc_method = "Raw", numCores = 1L, verbose = FALSE))
stopifnot(is.null(none_new$Data$Normalization), identical(none_new$Repeat$norm_objective, "mean_lower_tail"))
# Fill old positional formals through verbose, append the new argument last.
positional <- as.list(formals(CCS::subtypePerformance))[-1L]
positional <- lapply(positional, eval, envir = baseenv())
positional <- c(list(df), positional)
names(positional) <- names(formals(CCS::subtypePerformance))
positional$col_subtype <- "Subtype"
positional$col_tissue <- "tumor_type"
positional["n_norm_subtype"] <- list(NULL)
positional$numCores <- 1L
positional$verbose <- FALSE
positional$norm_search_method <- "dp"
positional$norm_objective <- "mean"
position_fit <- suppressWarnings(do.call(CCS::subtypePerformance, unname(positional)))
stopifnot(identical(position_fit$Repeat$norm_objective, "mean"))
grDevices::dev.off()
cat("Public subtypePerformance DP and NULL contracts passed.\n")
