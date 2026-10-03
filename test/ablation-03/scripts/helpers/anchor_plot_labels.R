# Display names for biological anchors; scientific keys remain unchanged.
.ablation03_anchor_labels <- function(language = "en") {
  if (identical(language, "zh")) {
    return(c(proliferation = "增殖", immune_tme = "免疫微环境",
      stromal_tme = "基质微环境", ifn_il6 = "IFN / IL6"))
  }
  c(proliferation = "Proliferation", immune_tme = "Immune TME",
    stromal_tme = "Stromal TME", ifn_il6 = "IFN / IL6")
}
