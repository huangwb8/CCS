# Purpose: Load the complete CCS expression atlas, frozen model, and technical metadata.
# Input: PADv20240810 expression RDS, one filtered resCCS, and the BatchInfo workbook.
# Parameters: Environment variables can override machine-specific roots and model hash.
# Output: Full data plus auditable reference/external partitions; no business filtering.

luckyBase::Plus.library(c("CCS", "readxl", "digest"))

# Step 1: Resolve cross-platform roots without embedding user-specific directories.
sysname <- Sys.info()[["sysname"]]
default_roots <- switch(
  sysname,
  Windows = list(data = "E:/", sync = "E:/Sync/"),
  Darwin = list(
    data = "/Volumes/2T01/winE/",
    sync = "/Volumes/share/docker/Sync/data/Sync/"
  ),
  stop(sprintf("Unsupported OS: %s", sysname), call. = FALSE)
)

root_path <- Sys.getenv("CCS_DATA_ROOT", unset = default_roots$data)
root_path_sync <- Sys.getenv("CCS_SYNC_ROOT", unset = default_roots$sync)
n_cores <- as.integer(Sys.getenv(
  "CCS_ABLATION_CORES",
  unset = as.character(max(1L, min(8L, parallel::detectCores(logical = TRUE) - 2L)))
))

project_version <- "PADv20250720"
data_version <- "PADv20240810"
model_version <- "PADv20240911"
paramMD5 <- Sys.getenv(
  "CCS_ABLATION_PARAM_MD5",
  unset = "5ff3a2de76e6cf902e765e8224f9cb66"
)

project.dir <- file.path(
  root_path,
  "iProjects",
  "RCheck",
  "GSClassifier",
  "routine01",
  "ccs",
  project_version
)
model.dir <- Sys.getenv(
  "CCS_ABLATION_MODEL_ROOT",
  unset = file.path(
    root_path,
    "iProjects",
    "RCheck",
    "GSClassifier",
    "routine01",
    "ccs",
    model_version
  )
)
data_path <- file.path(
  root_path,
  "iProjects",
  "CCS_Data",
  "report",
  paste0(
    "DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_",
    data_version,
    ".rds"
  )
)
resccs_path <- file.path(
  root_path_sync,
  "Project",
  project_version,
  "models",
  paste0("resCCS_", paramMD5, ".rds")
)
batch_workbook_path <- file.path(
  "test",
  "pre-train-info",
  paste0(
    "DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_",
    data_version,
    "_BatchInfo.xlsx"
  )
)

required_paths <- c(data_path, resccs_path, batch_workbook_path, model.dir)
missing_paths <- required_paths[!file.exists(required_paths) & !dir.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "ablation-test-data: required path(s) not found: ",
    paste(missing_paths, collapse = ", "),
    call. = FALSE
  )
}

# Step 2: Load the complete expression data and frozen filtered CCS object.
data_all <- readRDS(data_path)
resCCS <- readRDS(resccs_path)
resCCS@Repeat$model.dir <- model.dir
filtered_cohorts <- as.character(resCCS@Data$filtered.cohort)

# Step 3: Read the workbook as the single source for sample identity and technology.
batch_info <- as.data.frame(
  readxl::read_excel(batch_workbook_path, sheet = "BatchInfo"),
  stringsAsFactors = FALSE
)
cohort_summary <- as.data.frame(
  readxl::read_excel(batch_workbook_path, sheet = "CohortSummary"),
  stringsAsFactors = FALSE
)

ablation_metadata <- batch_info
colnames(ablation_metadata)[colnames(ablation_metadata) == "cohort_id"] <- "cohort"
ablation_metadata$cohort_key <- paste(
  ablation_metadata$tissue,
  ablation_metadata$cohort,
  sep = "/"
)
ablation_metadata$analysis_set <- ifelse(
  ablation_metadata$cohort_key %in% filtered_cohorts,
  "external_query",
  "reference_atlas"
)
ablation_metadata$d1_provenance <- ifelse(
  ablation_metadata$analysis_set == "external_query",
  "external_frozen",
  "in_sample"
)

# cancer_type comes from the independent BatchInfo evidence layer and was not
# used to fit the frozen cohort bank. tissue remains the bank-aligned fallback.
ablation_metadata$biology <- ablation_metadata$cancer_type
ablation_metadata$anchor_role <- "independent"
ablation_metadata$duplicate_sample_id_global <- as.logical(
  ablation_metadata$duplicate_sample_id_global
)
ablation_metadata$duplicate_within_cohort <- as.logical(
  ablation_metadata$duplicate_within_cohort
)
ablation_metadata$needs_human_check <- as.logical(
  ablation_metadata$needs_human_check
)

# Step 4: Partition cohorts without altering any sample or expression value.
partition_data <- function(data, cohort_keys, keep) {
  result <- list()
  for (tissue in names(data)) {
    for (cohort in names(data[[tissue]])) {
      key <- paste(tissue, cohort, sep = "/")
      if ((key %in% cohort_keys) == keep) {
        result[[tissue]][[cohort]] <- data[[tissue]][[cohort]]
      }
    }
  }
  result[lengths(result) > 0]
}

external_data <- partition_data(data_all, filtered_cohorts, keep = TRUE)
reference_data <- partition_data(data_all, filtered_cohorts, keep = FALSE)

# Step 5: Fail loudly when the workbook and RDS identity layers drift apart.
cohort_data <- unlist(lapply(data_all, unname), recursive = FALSE)
raw_sample_count <- sum(vapply(
  cohort_data,
  function(leaf) {
    expr <- if (is.list(leaf) && !is.null(leaf$expr)) leaf$expr else leaf
    ncol(expr)
  },
  integer(1)
))
if (raw_sample_count != nrow(ablation_metadata)) {
  stop(
    "ablation-test-data: BatchInfo row count does not match expression samples.",
    call. = FALSE
  )
}
if (!setequal(filtered_cohorts, unique(
  ablation_metadata$cohort_key[ablation_metadata$analysis_set == "external_query"]
))) {
  stop(
    "ablation-test-data: filtered cohort keys do not match BatchInfo.",
    call. = FALSE
  )
}

luckyBase::LuckyVerbose(
  "ablation-test-data: loaded ",
  nrow(ablation_metadata),
  " rows across ",
  nrow(cohort_summary),
  " cohorts; ",
  sum(ablation_metadata$analysis_set == "external_query"),
  " rows belong to ",
  length(filtered_cohorts),
  " external filtered cohorts."
)
