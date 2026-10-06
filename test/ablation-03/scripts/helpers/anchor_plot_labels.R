# Display names for the five current biological anchors.
.ablation03_anchor_labels <- function(language = "en") {
  if (identical(language, "zh")) {
    return(c(proliferation = "增殖", immune_tme = "免疫微环境",
      stromal_tme = "基质微环境", ifn = "IFNγ", il6 = "IL6–JAK–STAT3"))
  }
  c(proliferation = "Proliferation", immune_tme = "Immune TME",
    stromal_tme = "Stromal TME", ifn = "IFNγ", il6 = "IL6–JAK–STAT3")
}
