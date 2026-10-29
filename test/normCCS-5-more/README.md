# normCCS 5–12 级搜索测试

本目录展示 CCS 的搜索与目标函数改造。0.8.5 引入显式 DP 搜索，0.8.6 将择优目标扩展为五种可选策略。所有数据由代码生成，不读取真实患者文件。

本次新增内容从 [目标函数报告](objective-report.html) 开始阅读，源码见 [objective-report.Rmd](objective-report.Rmd)，命令与结果见 [本轮验证说明](objective-verification.md)。原有性能报告保留为 0.8.5 的历史证据。

## 0.8.6 目标函数选项

`norm_objective` 追加在公开接口末尾，默认 `"median"` 完整保留旧结果与表结构，包括 pooled `All`。新增 `"cohort_median"`、`"mean"`、`"lower_tail"`、`"mean_lower_tail"` 只评价真实队列，每队列等权。下尾取最低 25% 概率质量，非整数边界按比例计入；混合目标为均值与下尾均值各占 0.5。

```r
fit <- CCS::subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 8, norm_search_method = "dp", norm_roc_method = "Raw",
  norm_objective = "mean_lower_tail", numCores = 1
)
fit$Data$Normalization$Summary[1, ]
```

新策略使用 `objectiveScore` 排序，同分先比较真实队列均值，再比较 `Info`；`medianROCAUC` 仍表示包含 `All` 的旧中位数。非默认目标写入 `Repeat`，下尾和混合设置同时记录。缺少 `Repeat$norm_objective` 表示遗留 `median`；NULL 归一化分支仍验证并记录新目标，但不生成搜索诊断。

新增测试与复现命令如下，均从仓库根目录执行：

```powershell
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/test-objectives.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/scripts/tests/objective-smoke.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/benchmark-objectives.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/diagnose-objective-gaps.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/render-report.R --objectives
```

`test-objectives.R` 可选 `baseline` 或 `objectives` 分段。`benchmark-objectives.R` 可选 `--section=quality`、`timing`、`scaling`；分段结果必须来自同一配置和源码。重新生成已存在结果需显式 `--overwrite`。`objective-benchmark-config.json` 声明正式规模，隔离 smoke 用同一入口降低输入规模并检查正式产品、代码与 raw 未变化。

`reference/baseline-0.8.5.rds` 冻结了 64 组旧输出；两个 `0.8.5.R` 文件保留旧函数供直接计时。正式新增产品为 `objective-quality.csv`、`objective-timing.csv`、`objective-scaling.csv`，另含配置、会话与源码指纹。`diagnose-objective-gaps.R` 对所有差距超过工程比较尺度的案例枚举切点邻域，保留初值覆盖、收敛轮数和所需边界变动数。报告分别解释目标选择带来的分级变化、同目标近似差距、新旧默认路径开销和 DP／穷举加速，不把不同目标的分数直接比较为科学优劣。

建议先阅读 [性能与算法报告](report.html)，源码见 [report.Rmd](report.Rmd)。[验证记录](verification.md) 汇总回归、公开入口和包检查结果。

## 主要改动

`norm_search_method = "dp"` 复用 RR/NRR 排序，把样本压缩为队列与状态的应答计数，用动态规划生成候选，再优化原有的 AUC 中位数。最终候选逐队列通过原 ROC 函数复核。该方法只用单核；没有新增 C++ 或并行依赖。它对加权平均 Raw AUC 子问题精确，对原中位数目标近似，不保证全局最优。

```r
source("test/normCCS-5-more/fixtures.R")
df <- make_norm_fixture(m = 20, cohorts = 4, seed = 42)
# 公开 Raw ROC 沿用原接口的数值状态编码。
df$Subtype <- as.integer(sub("S", "", df$Subtype))
fit <- CCS::subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 8, norm_search_method = "dp",
  norm_roc_method = "Raw", norm_rank_marker = "NRR",
  numCores = 1
)
fit$Data$Normalization$Summary[1, ]
fit$Data$Normalization$SearchDiagnostics
```

默认 `"bruteforce"` 仍返回所有切分的结果。DP 的 `Raw` 和 `Summary` 只包含最终实际复核的少量候选；`SearchDiagnostics` 记录切点、排序、初值、评分模式、精修轮数、计算时间及参考复核次数。过大的穷举会提前报错，避免分配巨型组合矩阵和结果表。

## 运行方法

从 CCS 仓库根目录，用 `C:/R/R-4.3.1/bin/Rscript.exe`：

```powershell
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/test-search.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/scripts/tests/public-api.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/benchmark.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/render-report.R
```

重复生成已存在的基准结果或 HTML，需在相应命令末尾显式加 `--overwrite`。全量组件测试也可分成 `test-search.R legacy 2`、`legacy 3`、`legacy 4`、`kernels`、`state` 执行，便于查看每一部分的退出状态。基准按顺序运行，不使用并行调度。已有包工程沿用 simple 流程，没有为此迁移其它分析的 targets 或补建 renv。

`compare-refinement.R` 可复现保存在 `reference/` 中的单切点初版结果；`benchmark-config.json` 列出正式数据规模和种子，隔离测试通过同一入口换用较小配置。

## 文件与证据

| 文件 | 用途 |
|---|---|
| `fixtures.R` | 固定种子合成数据、独立排序与切点枚举辅助函数 |
| `reference/subtype-norm-0.8.4.R` | 修改前的原始搜索函数，用作回归和性能基线 |
| `test-search.R` | 旧结果等价、DP 精确子问题、pROC 方向和并列、调用者状态测试 |
| `scripts/tests/public-api.R` | 实际运行公开函数及 `n_norm_subtype=NULL` 分支 |
| `scripts/tests/smoke-test.R` | 隔离产品目录的计算与报告轻量测试 |
| `benchmark.R` | 正确性差距、三次直接计时、高维搜索与分级人数检查 |
| `products/quality.csv` | 每个合成问题的完整穷举差距与精修收益 |
| `products/timings.csv` | 原版和 DP 的每次直接计时 |
| `products/scaling.csv` | 107／216 状态、20 队列、5–12 级的实际运行结果 |
| `products/group-support.csv` | 各队列各等级人数与应答人数 |
| `products/session-info.txt` | 实际 R 与依赖版本 |

固定合成数据上的加速和最优值匹配只说明工程行为。真实患者全量验证未运行；更高 K 的临床收益、跨队列稳定性和校准需要另行研究。
