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
    "test/ablation-03/config/biological-anchors.yml"
  )
  stopifnot(identical(actual_config$anchors$ifn$name, "IFN\u03b3 signaling"))
  Sys.setlocale("LC_CTYPE", old_locale)
}
cat("anchor signature contracts passed\n")

split_config <- yaml::read_yaml("test/ablation-03/config/biological-anchors.yml")
split_signatures <- list(
  `Conserved-PanCan-TME-subtypes` = list(`Tumor proliferation rate` = c("p1", "p2")),
  Zeng2021 = list(TME_A_Immune2 = c("i1", "i2"), TME_B_Stromal2 = c("s1", "s2")),
  `IFN-IL6` = stats::setNames(list(c("shared", "ifn_only"),
    c("shared", "il6_only")), c("IFN\u03b3 signaling", "IL6-JAK-STAT3 signaling")))
selected <- .biology_select_anchors(split_config, split_signatures)
stopifnot(length(selected) == 5L, !"ifn_il6" %in% names(selected),
  setequal(selected$ifn, c("shared", "ifn_only")),
  setequal(selected$il6, c("shared", "il6_only")))
