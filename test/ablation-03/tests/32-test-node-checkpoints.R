# Verify multi-version checkpoint hits, invalidation and corruption recovery.
source(file.path("R", "ablation.R"), local = TRUE)

cache_root <- tempfile("ablation-node-cache-")
dir.create(cache_root)

prepared <- list(
  input_key = "input-v1",
  cache_key = "representation-v1",
  direct_cache = list(key = "direct-v1")
)
key_a <- .ablation_node_cache_key(
  "readout", prepared, list(lambda = 1), 11L, "test-v1"
)
key_b <- .ablation_node_cache_key(
  "readout", prepared, list(lambda = 2), 11L, "test-v1"
)
stopifnot(!identical(key_a, key_b))
stopifnot(!identical(
  key_a,
  .ablation_node_cache_key("readout", prepared, list(lambda = 1), 12L, "test-v1")
))
stopifnot(!identical(
  key_a,
  .ablation_node_cache_key("readout", prepared, list(lambda = 1), 11L, "test-v2")
))
original_readout <- .ablation_linear_readout
.ablation_linear_readout <- function(...) "changed-for-key-test"
key_code_changed <- .ablation_node_cache_key(
  "readout", prepared, list(lambda = 1), 11L, "test-v1"
)
.ablation_linear_readout <- original_readout
stopifnot(!identical(key_a, key_code_changed))
original_transform <- .ablation_module_balanced_transform
.ablation_module_balanced_transform <- function(...) "changed-transitive-helper"
key_transitive_changed <- .ablation_node_cache_key(
  "readout", prepared, list(lambda = 1), 11L, "test-v1"
)
.ablation_module_balanced_transform <- original_transform
stopifnot(!identical(key_a, key_transitive_changed))

# Unrelated node implementations must not evict a valid readout cache.
original_decoder <- .ablation_decode_direct_features
.ablation_decode_direct_features <- function(...) "unrelated-decoder-change"
key_unrelated_changed <- .ablation_node_cache_key(
  "readout", prepared, list(lambda = 1), 11L, "test-v1"
)
.ablation_decode_direct_features <- original_decoder
stopifnot(identical(key_a, key_unrelated_changed))

# Changing learning-curve-only parameters leaves retrieval/readout identities
# untouched while producing a distinct learning-curve identity.
retrieval_key <- .ablation_node_cache_key(
  "retrieval", prepared, list(k = 5L), 11L, "retrieval-v1"
)
readout_key <- .ablation_node_cache_key(
  "readout", prepared, list(lambda = 1), 11L, "readout-v1"
)
learning_key_a <- .ablation_node_cache_key(
  "learning-curve", prepared, list(fractions = 0.5), 11L, "learning-v1"
)
learning_key_b <- .ablation_node_cache_key(
  "learning-curve", prepared, list(fractions = 0.75), 11L, "learning-v1"
)
stopifnot(
  identical(retrieval_key, .ablation_node_cache_key(
    "retrieval", prepared, list(k = 5L), 11L, "retrieval-v1"
  )),
  identical(readout_key, .ablation_node_cache_key(
    "readout", prepared, list(lambda = 1), 11L, "readout-v1"
  )),
  !identical(learning_key_a, learning_key_b)
)

# Runtime observability changes must not invalidate scientific node keys.  The
# helpers below only update checkpoint metadata (heartbeat, process identity,
# memory, and error summaries); they are intentionally outside the scientific
# implementation graph.
operational_helper_names <- c(
  ".ablation_checkpoint_heartbeat",
  ".ablation_runtime_identity",
  ".ablation_process_memory",
  ".ablation_condition_summary"
)
operational_key_before <- .ablation_node_cache_key(
  "learning-curve", prepared, list(fractions = 0.5), 11L, "learning-v1"
)
operational_originals <- lapply(
  operational_helper_names,
  get,
  envir = environment(.ablation_node_code_identity),
  inherits = FALSE
)
assign(
  ".ablation_checkpoint_heartbeat",
  function(...) invisible(FALSE),
  envir = environment(.ablation_node_code_identity)
)
assign(
  ".ablation_runtime_identity",
  function() list(run_id = "changed-runtime", pid = 1L, hostname = "changed"),
  envir = environment(.ablation_node_code_identity)
)
assign(
  ".ablation_process_memory",
  function() list(working_set_bytes = 0, peak_working_set_bytes = 0),
  envir = environment(.ablation_node_code_identity)
)
assign(
  ".ablation_condition_summary",
  function(condition, limit = 1000L) "changed-summary",
  envir = environment(.ablation_node_code_identity)
)
operational_key_after <- .ablation_node_cache_key(
  "learning-curve", prepared, list(fractions = 0.5), 11L, "learning-v1"
)
stopifnot(identical(operational_key_before, operational_key_after))
invisible(Map(assign, operational_helper_names, operational_originals,
              MoreArgs = list(envir = environment(.ablation_node_code_identity))))

counter <- 0L
first <- .ablation_cached_node(
  "readout", cache_root, key_a,
  compute = function() {
    counter <<- counter + 1L
    list(answer = 42L)
  },
  verbose = FALSE
)
second <- .ablation_cached_node(
  "readout", cache_root, key_a,
  compute = function() stop("matching cache should have been reused"),
  verbose = FALSE
)
stopifnot(counter == 1L, first$status == "written", second$status == "hit")
stopifnot(
  first$lookup_status == "miss", first$reason == "missing-cache",
  second$lookup_status == "hit", second$reason == "valid"
)
stopifnot(identical(first$value, second$value))

cache_path <- file.path(cache_root, second$path)
state_path <- sub("\\.rds$", ".state.rds", cache_path)
expect_rebuild <- function(mutate, expected_reason, answer) {
  mutate()
  rebuilt <- .ablation_cached_node(
    "readout", cache_root, key_a,
    compute = function() {
      counter <<- counter + 1L
      list(answer = answer)
    },
    verbose = FALSE
  )
  stopifnot(
    rebuilt$status == "written", rebuilt$lookup_status == "miss",
    rebuilt$reason == expected_reason, rebuilt$value$answer == answer
  )
  rebuilt
}

third <- expect_rebuild(function() {
  writeLines("corrupted", cache_path)
}, "cache-file-hash-mismatch", 43L)

invisible(expect_rebuild(function() {
  state <- readRDS(state_path)
  state$status <- "running"
  state$pid <- 999999L
  state$hostname <- unname(Sys.info()[["nodename"]])
  .ablation_atomic_save_rds(state, state_path)
}, "stale-running-state", 44L))

invisible(expect_rebuild(function() {
  cached <- readRDS(cache_path)
  cached$schema_version <- 999L
  .ablation_atomic_save_rds(cached, cache_path)
  state <- readRDS(state_path)
  state$cache_md5 <- unname(tools::md5sum(cache_path))
  .ablation_atomic_save_rds(state, state_path)
}, "cache-schema-mismatch", 45L))

invisible(expect_rebuild(function() {
  cached <- readRDS(cache_path)
  cached$key <- "wrong-key"
  .ablation_atomic_save_rds(cached, cache_path)
  state <- readRDS(state_path)
  state$cache_md5 <- unname(tools::md5sum(cache_path))
  .ablation_atomic_save_rds(state, state_path)
}, "cache-key-mismatch", 46L))

invisible(expect_rebuild(function() {
  cached <- readRDS(cache_path)
  cached$value$answer <- 999L
  .ablation_atomic_save_rds(cached, cache_path)
  state <- readRDS(state_path)
  state$cache_md5 <- unname(tools::md5sum(cache_path))
  .ablation_atomic_save_rds(state, state_path)
}, "cache-value-hash-mismatch", 47L))
stopifnot(counter == 6L)

state <- readRDS(state_path)
stopifnot(
  identical(state$schema_version, 2L),
  identical(state$status, "complete"), identical(state$key, key_a),
  is.numeric(state$working_set_start_bytes),
  is.numeric(state$peak_working_set_bytes),
  state$result_bytes > 0
)

failed_key <- digest::digest("failed-node", algo = "md5")
failed <- tryCatch(
  .ablation_cached_node(
    "readout", cache_root, failed_key,
    compute = function() stop("synthetic checkpoint failure"),
    verbose = FALSE
  ),
  error = identity
)
stopifnot(inherits(failed, "error"))
failed_state <- readRDS(file.path(
  cache_root, "checkpoints", "readout", paste0(failed_key, ".state.rds")
))
stopifnot(
  identical(failed_state$status, "failed"),
  identical(failed_state$error_class, "simpleError"),
  grepl("synthetic checkpoint failure", failed_state$error_summary, fixed = TRUE)
)

# Long jobs expose progress without changing their scientific result.  The
# heartbeat must update only the owning running state and remain readable by a
# later resume process.
heartbeat_key <- digest::digest("heartbeat-node", algo = "md5")
heartbeat_path <- file.path(
  cache_root, "checkpoints", "learning-curve-job",
  paste0(heartbeat_key, ".state.rds")
)
.ablation_atomic_save_rds(
  list(
    schema_version = 2L, status = "running",
    node = "learning-curve-job", key = heartbeat_key,
    run_id = "heartbeat-run", pid = Sys.getpid(),
    hostname = unname(Sys.info()[["nodename"]]),
    started_at = "2000-01-01T00:00:00+0000",
    updated_at = "2000-01-01T00:00:00+0000",
    heartbeat_count = 0L
  ),
  heartbeat_path
)
stopifnot(isTRUE(.ablation_checkpoint_heartbeat(
  cache_root, "learning-curve-job", heartbeat_key,
  list(stage = "inner-cv-xgb-complete", index = 2L, total = 13L)
)))
heartbeat_state <- readRDS(heartbeat_path)
stopifnot(
  identical(heartbeat_state$status, "running"),
  identical(heartbeat_state$heartbeat_count, 1L),
  identical(heartbeat_state$progress$stage, "inner-cv-xgb-complete"),
  identical(heartbeat_state$progress$index, 2L),
  identical(heartbeat_state$progress$total, 13L),
  heartbeat_state$updated_at != "2000-01-01T00:00:00+0000"
)

foreign_key <- digest::digest("foreign-heartbeat", algo = "md5")
foreign_path <- file.path(
  cache_root, "checkpoints", "learning-curve-job",
  paste0(foreign_key, ".state.rds")
)
.ablation_atomic_save_rds(
  list(
    schema_version = 2L, status = "running",
    node = "learning-curve-job", key = foreign_key,
    run_id = "foreign-run", pid = Sys.getpid() + 1L,
    hostname = unname(Sys.info()[["nodename"]]),
    heartbeat_count = 0L
  ),
  foreign_path
)
stopifnot(!isTRUE(.ablation_checkpoint_heartbeat(
  cache_root, "learning-curve-job", foreign_key,
  list(stage = "must-not-overwrite")
)))
foreign_state <- readRDS(foreign_path)
stopifnot(identical(foreign_state$heartbeat_count, 0L))

stale_key <- digest::digest("stale-node", algo = "md5")
stale_dir <- file.path(cache_root, "checkpoints", "readout")
.ablation_atomic_save_rds(
  list(
    schema_version = 2L, status = "running", node = "readout", key = stale_key,
    run_id = "dead-run", pid = 999999L, hostname = unname(Sys.info()[["nodename"]]),
    started_at = "2000-01-01T00:00:00+0000", updated_at = "2000-01-01T00:00:00+0000"
  ),
  file.path(stale_dir, paste0(stale_key, ".state.rds"))
)
stale_result <- .ablation_cached_node(
  "readout", cache_root, stale_key,
  compute = function() list(answer = 48L), verbose = FALSE
)
stale_state <- readRDS(file.path(stale_dir, paste0(stale_key, ".state.rds")))
stopifnot(
  stale_result$value$answer == 48L,
  identical(stale_result$reason, "stale-running-state"),
  identical(stale_state$status, "complete"),
  identical(stale_state$recovered_from$status, "stale"),
  identical(stale_state$recovered_from$stale_reason, "owner-process-not-active")
)

legacy_key <- digest::digest("legacy-running-node", algo = "md5")
.ablation_atomic_save_rds(
  list(schema_version = 1L, status = "running", node = "readout", key = legacy_key),
  file.path(stale_dir, paste0(legacy_key, ".state.rds"))
)
legacy_result <- .ablation_cached_node(
  "readout", cache_root, legacy_key,
  compute = function() list(answer = 49L), verbose = FALSE
)
stopifnot(
  legacy_result$value$answer == 49L,
  identical(legacy_result$reason, "stale-running-state")
)

active_key <- digest::digest("active-node", algo = "md5")
.ablation_atomic_save_rds(
  list(
    schema_version = 2L, status = "running", node = "readout", key = active_key,
    run_id = "active-run", pid = Sys.getpid(), hostname = unname(Sys.info()[["nodename"]])
  ),
  file.path(stale_dir, paste0(active_key, ".state.rds"))
)
active_error <- tryCatch(
  .ablation_cached_node(
    "readout", cache_root, active_key,
    compute = function() list(answer = 49L), verbose = FALSE
  ),
  error = identity
)
stopifnot(inherits(active_error, "error"), grepl("active", conditionMessage(active_error), fixed = TRUE))

fit_path <- file.path(cache_root, "fit-cache.rds")
.ablation_atomic_save_rds(
  list(
    schema_version = 1L,
    status = "complete",
    key = "fit-v1",
    direct = list(a = list(
      overall = data.frame(balanced_accuracy = 0.5),
      selected_lambda = 1
    )),
    d1 = list()
  ),
  fit_path
)
stopifnot(!is.null(.ablation_read_fit_cache(fit_path, "fit-v1")))
.ablation_atomic_save_rds(
  list(
    schema_version = 1L,
    status = "complete",
    key = "fit-v1",
    direct = list(a = list(overall = "damaged", selected_lambda = 1)),
    d1 = list()
  ),
  fit_path
)
stopifnot(is.null(.ablation_read_fit_cache(fit_path, "fit-v1")))
writeLines("corrupted", fit_path)
stopifnot(is.null(.ablation_read_fit_cache(fit_path, "fit-v1")))
unlink(cache_root, recursive = TRUE, force = TRUE)
cat("node checkpoint tests passed\n")
