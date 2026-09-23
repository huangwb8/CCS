local({
  if (.Platform$OS.type == "windows") {
    previous_locale <- Sys.getlocale("LC_CTYPE")
    on.exit(Sys.setlocale("LC_CTYPE", previous_locale), add = TRUE)
    Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")
  }
  project <- if (file.exists(file.path("test", "ablation-03", "_targets.R"))) {
    file.path(getwd(), "test", "ablation-03")
  } else getwd()
  previous_directory <- getwd()
  setwd(project)
  on.exit(setwd(previous_directory), add = TRUE)
  source(file.path("targets", "functions.R"), local = TRUE)

  cache_root <- tempfile("ablation03-biology-branch-")
  dir.create(file.path(cache_root, "01-representations"), recursive = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE), add = TRUE)
  contract_path <- file.path(cache_root, "01-representations", "sample-contract.rds")
  contract <- list(
    reference = data.frame(sample_id = "reference-1"),
    query = data.frame(sample_id = "query-1")
  )
  saveRDS(contract, contract_path)
  runtime_config <- list(cache_root = cache_root)
  cache <- list(
    schema_version = 2L, status = "complete", sample_key_hash = "stale",
    structural_anchor_cache = list(status = "complete")
  )
  rebuilt <- FALSE
  .ablation03_target_stage <- function(...) {
    rebuilt <<- TRUE
    list(stage = "rebuilt")
  }
  result <- .ablation03_biology_target(
    list(input = list(biology_inputs = cache)),
    representation_target = list(), runtime_config = runtime_config
  )
  stopifnot(rebuilt, identical(result$stage, "rebuilt"))
  stopifnot(!file.exists(file.path(cache_root, "01-biology", "expression-anchor-cache.rds")))

  source_path <- file.path(cache_root, "atlas.rds")
  signature_path <- file.path(cache_root, "signature.rds")
  saveRDS(1, source_path)
  saveRDS(2, signature_path)
  old_source <- Sys.getenv("CCS_FULL_EXPRESSION_RDS", unset = NA_character_)
  old_signature <- Sys.getenv("CCS_GENE_SIGNATURE_RDS", unset = NA_character_)
  on.exit({
    if (is.na(old_source)) Sys.unsetenv("CCS_FULL_EXPRESSION_RDS") else Sys.setenv(CCS_FULL_EXPRESSION_RDS = old_source)
    if (is.na(old_signature)) Sys.unsetenv("CCS_GENE_SIGNATURE_RDS") else Sys.setenv(CCS_GENE_SIGNATURE_RDS = old_signature)
  }, add = TRUE)
  Sys.setenv(CCS_FULL_EXPRESSION_RDS = source_path, CCS_GENE_SIGNATURE_RDS = signature_path)
  cache$sample_key_hash <- digest::digest(
    paste(sort(unique(c(contract$reference$sample_id, contract$query$sample_id))), collapse = "\n"),
    algo = "md5", serialize = FALSE
  )
  cache$sample_contract_md5 <- digest::digest(file = contract_path, algo = "md5")
  cache$builder_md5 <- digest::digest(file = "01.03.00. 生物输入准备.R", algo = "md5")
  cache$source <- list(path = normalizePath(source_path, winslash = "/"),
    md5 = digest::digest(file = source_path, algo = "md5"))
  cache$signature <- list(path = normalizePath(signature_path, winslash = "/"),
    md5 = digest::digest(file = signature_path, algo = "md5"),
    config_md5 = digest::digest(file = file.path("raw", "config", "biological-anchors.yml"), algo = "md5"))
  rebuilt <- FALSE
  .ablation03_biology_target(
    list(input = list(biology_inputs = cache)),
    representation_target = list(), runtime_config = runtime_config
  )
  stopifnot(!rebuilt)
  stopifnot(file.exists(file.path(cache_root, "01-biology", "expression-anchor-cache.rds")))

  cache$builder_md5 <- "stale"
  rebuilt <- FALSE
  .ablation03_biology_target(
    list(input = list(biology_inputs = cache)),
    representation_target = list(), runtime_config = runtime_config
  )
  stopifnot(rebuilt)

  rebuilt <- FALSE
  .ablation03_biology_target(
    list(input = list(biology_inputs = NULL)),
    representation_target = list(), runtime_config = runtime_config
  )
  stopifnot(rebuilt)
})

message("Precomputed biology cache branch isolation passed.")
