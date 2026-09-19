source("test/ablation-03/02.02.00. 生物锚点分析_functions.R")
config <- list(anchors = list(a = list(family = "F", names = c("X", "Y"))))
signatures <- list(F = list(Y = c("g2", "g3"), unused = "g4", X = c("g1", "g2")))
stopifnot(identical(.biology_select_anchors(config, signatures)$a, c("g1", "g2", "g3")))
signatures$F$X <- NULL
stopifnot(inherits(tryCatch(.biology_select_anchors(config, signatures), error = identity), "error"))
if (.Platform$OS.type == "windows") {
  old_locale <- Sys.getlocale("LC_CTYPE")
  Sys.setlocale("LC_CTYPE", "C")
  production <- parse("test/ablation-03/01.03.00. 生物输入准备.R", encoding = "UTF-8")
  # Execute the production Windows locale initializer without loading data.
  eval(production[[2L]])
  actual_config <- yaml::read_yaml(
    "test/ablation-03/raw/config/biological-anchors.yml"
  )
  stopifnot(identical(actual_config$anchors$ifn_il6$names[1L], "IFN\u03b3 signaling"))
  Sys.setlocale("LC_CTYPE", old_locale)
}
cat("anchor signature contracts passed\n")
