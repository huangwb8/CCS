# Reproduce the first-iteration quality results, without a separate search runner.
source("test/normCCS-5-more/00.Environment.R", encoding = "UTF-8")
devtools::load_all(quiet = TRUE)
env <- new.env(parent = asNamespace("CCS"))
sys.source("test/normCCS-5-more/reference/subtype-norm-dp-single-cut.R", env)
recorded <- read.csv("test/normCCS-5-more/reference/quality-single-cut.csv")
for (i in seq_len(nrow(recorded))) {
  r <- recorded[i, ]
  df <- make_norm_fixture(r$M, r$C, r$seed, sparse = r$seed %% 2L == 0L)
  ranked <- reference_rank(df, r$rank)
  fit <- env$subtype_norm_dp_single_cut(df, ranked, r$K, r$method, proc.time()[["elapsed"]], FALSE)
  stopifnot(abs(fit$Summary$medianROCAUC[1L] - r$dp_best) <= 1e-10)
}
cat("All 60 archived single-cut quality results reproduced.\n")
