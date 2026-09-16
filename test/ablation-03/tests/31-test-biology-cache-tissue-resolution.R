# Exercise the production extraction loop with two historically Undefined cohorts.
# Restoring the output tissue must not change the outer atlas lookup key.
code <- parse("test/ablation-03/01-ablation03-biology-cache.R", encoding = "UTF-8")
loops <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("for")), as.list(code))
stopifnot(length(loops) == 1L)
fixture <- new.env(parent = baseenv())
fixture$atlas <- list(Undefined = list(
  c1 = matrix(1:6, 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2"))),
  c2 = matrix(7:12, 3, dimnames = list(c("g1", "g2", "g3"), c("s3", "s4")))))
fixture$cohort_lookup <- data.frame(cohort = c("c1", "c2"), cohort_key = c("A/c1", "B/c2"))
fixture$metadata <- data.frame(sample_id = paste0("s", 1:4), cohort_key = rep(c("A/c1", "B/c2"), each = 2))
fixture$required_genes <- c("g1", "g2")
fixture$anchors <- list(a = c("g1", "g2"))
fixture$anchor_names <- "a"
fixture$cohorts <- fixture$coverage <- fixture$missing <- list()
fixture$index <- 0L
eval(loops[[1L]], fixture)
stopifnot(identical(names(fixture$cohorts), c("A/c1", "B/c2")),
  identical(fixture$cohorts[["B/c2"]]$tissue, "B"),
  identical(fixture$cohorts[["B/c2"]]$sample_id, c("s3", "s4")),
  identical(fixture$cohorts[["B/c2"]]$expression, fixture$atlas$Undefined$c2[1:2, , drop = FALSE]))
cat("PASS biology cache preserves raw atlas keys while restoring output tissues\n")
