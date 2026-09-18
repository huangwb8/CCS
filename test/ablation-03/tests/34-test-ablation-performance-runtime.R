# Regression checks for the deterministic learning-curve job scheduler.
source(file.path("R", "ablation.R"))

set.seed(1L)
n_train <- 16L
train_metadata <- data.frame(
  sample_id = paste0("r", seq_len(n_train)),
  cohort = rep(paste0("c", 1:8), each = 2L),
  cancer_type = rep(c("A", "B"), 8L),
  stringsAsFactors = FALSE
)
test_metadata <- data.frame(
  sample_id = paste0("q", 1:8),
  cohort = paste0("q", 1:8),
  cancer_type = rep(c("A", "B"), 4L),
  stringsAsFactors = FALSE
)
representations <- list(
  `Direct-GSClassifier` = list(
    train = matrix(rnorm(n_train * 4L), n_train, 4L,
      dimnames = list(train_metadata$sample_id, paste0("x", 1:4))),
    test = matrix(rnorm(8L * 4L), 8L, 4L,
      dimnames = list(test_metadata$sample_id, paste0("x", 1:4))),
    blocks = NULL
  ),
  `Cohort-d1` = list(
    train = matrix(rnorm(n_train * 6L), n_train, 6L,
      dimnames = list(train_metadata$sample_id, paste0("x", 1:6))),
    test = matrix(rnorm(8L * 6L), 8L, 6L,
      dimnames = list(test_metadata$sample_id, paste0("x", 1:6))),
    blocks = list(m1 = 1:3, m2 = 4:6)
  )
)

run_curve <- function(workers) {
  .ablation_learning_curve(
    representations = representations,
    train_metadata = train_metadata,
    test_metadata = test_metadata,
    label_column = "cancer_type",
    fractions = c(0.5, 1),
    repeats = 1L,
    lambda = c(0.1, 1),
    inner_folds = 2L,
    nrounds = 2L,
    numCores = 1L,
    workers = workers,
    seed = 10L
  )
}

serial <- run_curve(1L)
parallel <- run_curve(2L)
stopifnot(identical(serial$metrics, parallel$metrics))
stopifnot(identical(serial$paired, parallel$paired))

checkpoint_dir <- tempfile("ablation-learning-curve-")
dir.create(checkpoint_dir)
on.exit(unlink(checkpoint_dir, recursive = TRUE, force = TRUE), add = TRUE)
checkpointed <- .ablation_learning_curve(
  representations = representations,
  train_metadata = train_metadata,
  test_metadata = test_metadata,
  label_column = "cancer_type",
  fractions = c(0.5, 1),
  repeats = 1L,
  lambda = c(0.1, 1),
  inner_folds = 2L,
  nrounds = 2L,
  numCores = 1L,
  workers = 2L,
  checkpoint_output_dir = checkpoint_dir,
  checkpoint_key = "fixed-parent-key",
  seed = 10L
)
job_files <- list.files(
  file.path(checkpoint_dir, "checkpoints", "learning-curve-job"),
  pattern = "^[a-f0-9]+\\.rds$",
  full.names = TRUE
)
stopifnot(length(job_files) == 4L)
job_md5 <- tools::md5sum(job_files)
resumed <- .ablation_learning_curve(
  representations = representations,
  train_metadata = train_metadata,
  test_metadata = test_metadata,
  label_column = "cancer_type",
  fractions = c(0.5, 1),
  repeats = 1L,
  lambda = c(0.1, 1),
  inner_folds = 2L,
  nrounds = 2L,
  numCores = 1L,
  workers = 1L,
  checkpoint_output_dir = checkpoint_dir,
  checkpoint_key = "fixed-parent-key",
  seed = 10L
)
stopifnot(identical(checkpointed$metrics, resumed$metrics))
stopifnot(identical(job_md5, tools::md5sum(job_files)))

config <- .ablation_resolve_representation_config(
  10L,
  list(validation = list(workers = 2L))
)
stopifnot(config$validation$workers == 2L)
bad <- try(
  .ablation_resolve_representation_config(
    10L,
    list(validation = list(workers = 0L))
  ),
  silent = TRUE
)
stopifnot(inherits(bad, "try-error"))

cat("ablation performance runtime tests passed\n")
