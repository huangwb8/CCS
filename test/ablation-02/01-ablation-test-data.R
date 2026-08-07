# Purpose: Load and audit the complete CCS expression atlas, frozen model, and metadata.
# Input: PADv20240810 expression RDS, one filtered resCCS, and the BatchInfo workbook.
# Parameters: Environment variables can override machine-specific roots and model hash.
# Output: Full data plus auditable reference/external partitions; no business filtering.

luckyBase::Plus.library(c("CCS", "readxl", "digest"))
source(file.path("test", "ablation-02", "01-ablation-test-data_functions.R"))

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
tissue_mapping_path <- file.path(
  "test",
  "pre-train-info",
  "scripts",
  "tissue_mapping.csv"
)

required_paths <- c(
  data_path,
  resccs_path,
  batch_workbook_path,
  tissue_mapping_path,
  model.dir
)
missing_paths <- required_paths[!file.exists(required_paths) & !dir.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "ablation-test-data: required path(s) not found: ",
    paste(missing_paths, collapse = ", "),
    call. = FALSE
  )
}

# Step 2: Load the complete expression data and restore audited tissue groups.
data_all <- readRDS(data_path)
tissue_mapping <- utils::read.csv(
  tissue_mapping_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
raw_cohort_index <- .atd_cohort_index(data_all)
tissue_resolution <- .atd_resolve_tissue_groups(data_all, tissue_mapping)
data_all <- tissue_resolution$data
tissue_resolution_audit <- tissue_resolution$audit
resolved_cohort_index <- .atd_cohort_index(data_all)

resCCS <- readRDS(resccs_path)
resCCS@Repeat$model.dir <- model.dir
filtered_model_cohorts <- as.character(resCCS@Data$filtered.cohort)
filtered_cohorts <- .atd_resolve_cohort_keys(
  filtered_model_cohorts,
  tissue_mapping
)

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
ablation_metadata$bank_tissue <- ifelse(
  ablation_metadata$cohort %in% tissue_mapping$cohort_id,
  "Undefined",
  ablation_metadata$tissue
)
ablation_metadata$bank_cohort_key <- paste(
  ablation_metadata$bank_tissue,
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
# used to fit the frozen cohort bank. bank_tissue preserves the historical model key.
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

# Step 5: Fail loudly when expression, workbook, and frozen-bank identities drift.
raw_sample_count <- sum(resolved_cohort_index$sample_count)
if (raw_sample_count != nrow(ablation_metadata)) {
  stop(
    "ablation-test-data: BatchInfo row count does not match expression samples.",
    call. = FALSE
  )
}
if ("Undefined" %in% names(data_all) || any(ablation_metadata$tissue == "Undefined")) {
  stop(
    "ablation-test-data: audited tissue restoration left Undefined entries.",
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
if (!setequal(filtered_model_cohorts, unique(
  ablation_metadata$bank_cohort_key[
    ablation_metadata$analysis_set == "external_query"
  ]
))) {
  stop(
    "ablation-test-data: frozen-bank cohort keys do not match the model contract.",
    call. = FALSE
  )
}

# Step 6: Build full, unfiltered profiles for the human-readable report.
ablation_data_profile <- .atd_build_data_profile(
  raw_cohort_index = raw_cohort_index,
  resolved_cohort_index = resolved_cohort_index,
  metadata = ablation_metadata,
  tissue_resolution_audit = tissue_resolution_audit,
  filtered_model_cohorts = filtered_model_cohorts,
  filtered_cohorts = filtered_cohorts
)

data_output_dir <- file.path(
  "test", "ablation-02", "tmp", "ablation-experiment"
)
dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  ablation_data_profile,
  file.path(data_output_dir, "data-profile.rds")
)

luckyBase::LuckyVerbose(
  "01-ablation-test-data: loaded ",
  nrow(ablation_metadata),
  " rows across ",
  nrow(cohort_summary),
  " cohorts; ",
  nrow(tissue_resolution_audit),
  " historical Undefined cohorts were restored; ",
  sum(ablation_metadata$analysis_set == "external_query"),
  " rows belong to ",
  length(filtered_cohorts),
  " external filtered cohorts."
)
