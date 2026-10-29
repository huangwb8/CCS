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
grDevices::dev.off()
cat("Public subtypePerformance DP and NULL contracts passed.\n")
