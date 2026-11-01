# CCS: Cohort Congress System

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%203.6.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.8.7-brightgreen.svg)](https://github.com/huangwb8/CCS)

### A computational framework for personalized pan-cancer genomic classification

`CCS` is an R package for personalized pan-cancer genomic classification. It supports cohort-level probability modeling, subtype discovery, and downstream response-oriented analyses, and is designed to work within the `GSClassifier` ecosystem.

## Tutorial

For teaching materials and step-by-step workflow explanations, please visit [Online Tutorial](https://huangwb8.github.io/ccs.principle/). Please note that the tutorial will be made publicly available in its full, fully usable form only after the manuscript is formally published.

## Searching five or more normCCS grades

Use `subtypePerformance(..., n_norm_subtype = 8, norm_search_method = "dp")`
for deterministic, single-core search without enumerating every partition. The
default `"bruteforce"` retains the original complete candidate tables. DP optimizes
fixed-direction mean AUC subproblems and refines the original median AUC objective;
it does not guarantee the global median optimum. Its final candidates are checked
with the original ROC evaluator. `norm_roc_method = "Raw"` preserves the grade
ordering; the legacy `"OptimizeRank"` reorders grades using evaluation labels and
should not be used as a frozen external prediction score.

See the [synthetic benchmark and reproducible report](test/normCCS-5-more/README.md)
for exactness checks, measured speedups, resource usage, and the public API example.

### 选择 normCCS 的目标函数

`norm_objective = "median"` 默认复现包含 pooled `All` 的旧中位数目标。
显式选择 `"cohort_median"`、`"mean"`、`"lower_tail"` 或 `"mean_lower_tail"`
时，仅对真实队列 AUC 评分，队列等权。下尾取最低 25% 概率质量，边界队列按比例计入；
混合目标由均值和下尾均值各占 0.5。穷举和 DP 使用同一目标，DP 整体仍是近似搜索。

```r
fit <- subtypePerformance(df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 8, norm_search_method = "dp", norm_roc_method = "Raw",
  norm_objective = "mean_lower_tail", numCores = 1)
fit$Data$Normalization$Summary[1, ]
```

新增策略以 `objectiveScore` 择优；`medianROCAUC` 保留旧定义。缺少
`Repeat$norm_objective` 表示使用遗留 `median`。参见
[目标函数报告](test/normCCS-5-more/objective-report.html)和
[复现说明](test/normCCS-5-more/README.md)。新目标的训练内得分改善不等同于外部预测收益。

### 连续继承式 normCCS 搜索

`norm_mode = "discrete"` 默认只搜索指定 K，保留既有行为。选择
`"continuous"` 后，先固定 metaCCS 状态顺序，再探索 K=1 到指定最大值。
K1 穷举唯一映射；K≥2 在上一层 top 方案内做受约束 DP，保留父方案的全部
分界，只增加一个分界，并用实际目标核对全部合法单次细分。
`norm_top_n` 控制每层保留数，默认 10；连续模式固定使用单核继承式 DP，
无需另设 `norm_search_method`。

```r
fit <- subtypePerformance(df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 8, norm_mode = "continuous", norm_top_n = 10,
  norm_rank_marker = "ShrunkLogOR", norm_objective = "cohort_median",
  norm_roc_method = "Raw", numCores = 1)
normalization <- fit$Data$Normalization
normalization$Summary  # 各 K 保留方案的合并排名，沿用原列结构
normalization$SearchDiagnostics$layers[[5]]  # K5 的父方案、分界和得分
grades <- newSubtype(df$Subtype, normalization$Summary$Info[1])
```

`Raw` 与 `Summary` 保留各 K 的 top 方案，最终仍由原 ROC 逐队列复核。
映射语法与 `newSubtype()` 不变。K1 为恒定评分，AUC=0.5；有限 top 会剪掉部分
搜索路径，因此不保证全局最优，也不保证增加 K 后得分上升。

`norm_rank_marker` 在原四种规则外增加 `MeanNRR`（队列等权均值）、
`LowerTailNRR`（最低 25% 队列 NRR）、`ShrunkLogOR`（Jeffreys 平滑后的
状态与队列总体 log-odds 差的均值）和 `RankFusion`（中位 RR/NRR 的等权秩融合）。
缺失状态不在该队列补零；新规则同分时按中位 NRR、中位 RR、状态名排序。
这些规则应在训练数据上估计并冻结，效果仍需独立验证。

运行合成测试：`Rscript --vanilla test/normCCS-continuous/run-tests.R`。
源码、DP 内核及目标函数统一位于 `R/phenotype.R`。

## Citation

The `CCS` manuscript is currently under review.

## AI-Assisted Development

Project-wide AI instructions live in [`AGENTS.md`](AGENTS.md); [`CLAUDE.md`](CLAUDE.md) is the Claude Code adapter. Human, AI, and tool contributions are recorded with [bensz-auto-contribution](https://github.com/huangwb8/bensz-auto-contribution) in [`docs/contribution.bac`](docs/contribution.bac). The ledger requires Python 3.10+ and can be verified with:

```bash
python -m bac --root . --bac-file docs/contribution.bac verify --json
```

## Authors

Weibin Huang (黄伟斌); <hwb2012@qq.com>; [Blog](https://blognas.hwb0307.com/); [ORCID: 0000-0003-0159-435X](https://orcid.org/0000-0003-0159-435X)

## License

`CCS` is released under the Apache-2.0 license. See [`license.txt`](license.txt) for details.
