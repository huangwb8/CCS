<div align="center">
  <h1>CCS · Cohort Congress System</h1>
  <p>
    <a href="https://github.com/huangwb8/CCS"><img src="https://img.shields.io/badge/version-0.8.9-brightgreen.svg" alt="版本 0.8.9"></a>
    <a href="https://www.r-project.org/"><img src="https://img.shields.io/badge/R-package-276DC3.svg" alt="R 包"></a>
    <a href="license.txt"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg" alt="Apache 2.0 许可证"></a>
  </p>
  <p><a href="README.md">English</a> · <a href="#快速开始">快速开始</a> · <a href="#api-与文档">文档</a> · <a href="https://github.com/huangwb8/CCS/issues">问题反馈</a></p>
  <p><strong>用于个性化泛癌基因组分类的计算框架。</strong></p>
</div>

CCS 是一个 R 包，通过一组基于 **GSClassifier** 训练的队列特异性模型，将样本表示为队列概率谱。在此基础上，CCS 支持亚型发现、新样本预测，以及跨队列的亚型与治疗应答关联分析。下游 **normCCS** 流程进一步将已有亚型状态归并为有序等级。

下方合成示例无需患者数据或预训练模型库，可用于体验 normCCS。完整的基因组分类训练则需要用户提供表达数据、基因集和基因注释。

## 快速开始

### 安装

项目开发与示例验证使用 **R 4.3.1**。[DESCRIPTION](DESCRIPTION) 未声明最低 R 版本，旧版 R 的兼容性尚未在本次任务中验证。通过 GitHub 安装前，请启用 Bioconductor 仓库，以解析 `ComplexHeatmap` 及 `WGCNA` 所需的依赖：

```r
install.packages(c("remotes", "BiocManager"))
options(repos = BiocManager::repositories())

remotes::install_github("huangwb8/luckyBase", upgrade = "never")
remotes::install_github("huangwb8/GSClassifier", upgrade = "never")
remotes::install_github("huangwb8/CCS", upgrade = "never")

library(CCS)
packageVersion("CCS")
```

依赖包括 `xgboost`、`torch` 和 `luz`。深度学习流程还需要 torch 运行时；使用该后端时，请参阅 [R torch 安装指南](https://torch.mlverse.org/docs/articles/installation.html)。安装失败时，应先确认具体失败的依赖及其仓库或系统要求，再重试。CCS 源码安装需要 C++17 编译器和 Rcpp；Windows 下请安装与 R 版本匹配的 Rtools。精确搜索在支持 OpenMP 时使用多线程，否则使用单线程。

### 运行合成 normCCS 分析

示例构造三个队列、八个亚型状态，共 480 条合成记录，应答率随状态编号递增。随后以单核搜索五级映射。

```r
library(CCS)

df <- expand.grid(
  Cohort = c("A", "B", "C"),
  Subtype = 1:8,
  replicate = 1:20,
  KEEP.OUT.ATTRS = FALSE
)
df$response <- as.integer(df$replicate <= 2 * df$Subtype + 1)
df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
df$tumor_type <- "Lung"

fit <- subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 5, norm_search_method = "dp",
  norm_roc_method = "Raw", norm_objective = "cohort_median",
  numCores = 1, verbose = FALSE
)

normalization <- fit$Data$Normalization
head(normalization$Summary)
grades <- newSubtype(df$Subtype, normalization$Summary$Info[1])
table(grades)
```

`Summary` 按 AUC 目标对候选映射排序，其中的 `Info` 列保存映射记录，可交给 `newSubtype()` 应用于亚型标签。本例将全部 480 条记录分为五级，等级为 0–4。搜索细节位于 `normalization$SearchDiagnostics`，图形位于 `fit$Plot`。合成结果用于演示用法，预测性能仍需在患者数据上验证。

## 从表达数据到亚型

CCS 将队列模型训练、概率表示、亚型发现和下游评价连接起来：

| 任务 | 主要入口 | 输入或输出 |
| --- | --- | --- |
| 准备训练数据 | `createCCSData()` | 按组织、队列嵌套的表达矩阵，以及基因集和基因注释 |
| 训练队列模型 | `ccs()` | 训练数据 → `CCS` S4 对象及保存的模型文件 |
| 探索并定义亚型 | `dr()`、`cluster()`、`optimizeCluster()` | 队列概率表示 → 降维结果与亚型分配 |
| 拟合亚型分类器并预测 | `scaller()`、`predict()` | 已完成聚类的 CCS 对象及新表达数据 |
| 评价应答并构建等级 | `subtypePerformance()`、`newSubtype()` | 亚型标签、队列标签及二分类应答 |
| 比较冻结表示 | `ablation()` | 已拟合的 CCS 对象、表达数据及样本注释 |

`ccs()` 将队列概率表示保存在 `object@Data$Probability$d1` 中，后续方法围绕 CCS 对象展开。完整训练与合成分级示例是两个独立流程，各函数的输入要求见 [API 文档](#api-与文档)。

应答分析需要完整的样本标识、亚型、队列、组织及 **0/1** 应答列，每个队列必须同时包含两类应答。本例使用默认列名 `SampleIDs`、`Subtype`、`Cohort` 和 `response`，并显式设置 `col_tissue = "tumor_type"`。

## 选择 normCCS 搜索方式

### 搜索固定等级数

默认的 `norm_mode = "discrete"` 只搜索指定等级数。默认 `norm_search_method = "bruteforce"` 保留完整候选表，但划分数可能快速增长。快速开始中的 `"dp"` 使用确定性、单核的动态规划与局部细化，仅返回实际评估的最终候选及搜索诊断，不保留全部划分。

DP 精确优化的是固定方向的加权平均 Raw AUC 代理目标。完整搜索还会围绕选定目标细化候选，并由原 ROC 评价器核对，因此 **不保证所选公开目标的全局最优解**。

K = 2–8 时，可以选择 `"stream_exact"`，用编译后的计数评分遍历固定状态顺序下的全部连续分段。它支持两种 ROC 方法和五种目标函数。搜索用有界堆默认保留最佳 `norm_exact_top_n = 2000L` 个映射，K ≥ 5 也适用；组合不足时全部返回。该参数与连续模式的 `norm_top_n` 独立。例如，复用上面的 `df`：

```r
exact_fit <- subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 5, norm_search_method = "stream_exact",
  norm_exact_top_n = 2000L, norm_roc_method = "Raw",
  norm_objective = "cohort_median", numCores = 4, verbose = FALSE
)
exact_fit$Data$Normalization$SearchDiagnostics
```

`verbose = TRUE` 时，`stream_exact` 会在搜索前报告组合数和原生搜索预计耗时；达到一亿个候选或预计一分钟时，还会提醒存在组合爆炸。减小 `norm_exact_top_n` 只减少保留结果，不会减少遍历量。设置 `verbose = FALSE` 可关闭这些消息。

耗时参考已有的 i9-14900K 基准：16 个 OpenMP 线程、20 个队列加 `All`、median 目标和 Top 2000。估计按队列数和请求线程数缩放，线程加速最多按 16 线程计算；OptimizeRank 另以 K6 基准乘以 K/6。这只是参考硬件下的粗估，不是置信区间或对当前机器的实测，预处理和核验还需额外时间。K > 8 即使显示了外推耗时，仍不受支持。

该搜索完整覆盖选定状态顺序下的组合，评分核验容差为 `1e-10`。所有保留候选的 AUC 都与 R 批量评分器核对，并选取最多五个映射调用原始 ROC 评价器复核；K = 2 时，全部保留映射均调用原始评价器，以保留二分类性能指标。浮点得分完全相同时，非默认目标先比较队列平均 AUC，再按 `Info` 的 UTF-8 字节序排序。原始 ROC 计算的末位差异可能改变数学同分映射的顺序。诊断记录完整覆盖、返回数量、线程、核验误差和耗时。后端最多支持 512 个状态，并设置内存保护；运行时间仍随 `choose(M - 1, K - 1)` 增长。`numCores = NULL` 使用可用核心数减一。

### 对精确搜索加入正则化

`stream_exact` 可以降低搜索对小样本等级的偏好。`norm_min_grade_n` 设置每一级的训练记录数底线，`norm_regularization` 控制软惩罚。两者都在完整搜索中生效，再按新得分保留 Top N。默认值分别为 `1L` 和 `0`，复现未正则化的搜索；非默认设置仅支持离散 `stream_exact`。

选择得分为所选 AUC 目标减去 `norm_regularization * mean(1 / sqrt(n_g))`，其中 `n_g` 是第 g 级的合并训练记录数。合并 `All` 只计一次，不要求每个队列或每类结局都达到人数底线。例如，复用上面的合成 `df`：

```r
regularized_fit <- subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 3, norm_search_method = "stream_exact",
  norm_roc_method = "Raw", norm_objective = "cohort_median",
  norm_min_grade_n = 20L, norm_regularization = 0.5,
  numCores = 4, verbose = FALSE
)
regularized_fit$Data$Normalization$Summary
```

这里的数值用于演示接口，不是经过验证的调参建议。原 AUC 列含义保持不变。任一约束启用后，`Summary` 新增 `supportPenalty`、`regularizationPenalty`、`regularizedScore`、`minGradeN` 和 `gradeSampleSizes`（以等级命名的人数向量列表列），按 `regularizedScore` 排序。诊断记录合法搜索空间和被排除的分组数。合法候选少于 Top N 时全部返回；没有合法映射时明确报错。精确最优只针对固定状态顺序、K 和人数约束下的正则化目标。

这是一种偏好充分样本支持、人数较均衡等级的启发式惩罚，可能压制真实稀有亚群，并非置信下界。人数底线和惩罚强度应通过训练侧分组验证选择，每次拟合都在对应训练折内重新学习状态排序与映射。输入前应处理重复患者和重叠来源。该选项不改变排序、Raw 自动方向或 OptimizeRank，因此仍须遵守下文的外部验证要求。上面的耗时估计来自未正则化基准，不包含正则化开销。

### 探索连续继承的等级序列

复用快速开始中的 `df`，探索 K = 1 至 K = 8：

```r
continuous_fit <- subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 8, norm_mode = "continuous", norm_top_n = 10,
  norm_rank_marker = "ShrunkLogOR", norm_objective = "cohort_median",
  norm_roc_method = "Raw", numCores = 1, verbose = FALSE
)

continuous_fit$Data$Normalization$Summary
continuous_fit$Data$Normalization$SearchDiagnostics$layers[[5]]
```

连续模式先固定状态顺序。K = 1 为恒定评分，AUC = 0.5；后续每个子方案保留父方案的全部分界，再新增一个分界。所有合法的单次细分按所选目标评分后，每层最多保留 `norm_top_n` 个候选。两种旧搜索设置在此均使用单核继承式 DP；与 `"stream_exact"` 同用会报错，因为精确后端独立搜索指定 K 的离散分段。

`Summary` 合并各 K 保留的候选，分层诊断记录父映射、分界和得分。有限候选剪枝使整体搜索仍为近似过程，不保证全局最优，也不保证增加 K 后得分上升。

### 选择评分目标

| `norm_objective` | 择优依据 |
| --- | --- |
| `"median"`（默认） | 包含合并队列 `All` 的旧中位 AUC |
| `"cohort_median"` | 真实队列 AUC 的中位数 |
| `"mean"` | 队列等权的平均 AUC |
| `"lower_tail"` | 最低 25% 队列 AUC 概率质量的均值，边界队列按比例计入 |
| `"mean_lower_tail"` | 队列均值与下尾均值的等权平均 |

非默认目标排除合并队列 `All`，要求每个真实队列均有有限 AUC，并按 `objectiveScore` 排序；`medianROCAUC` 保留旧定义。离散和连续搜索均使用选定目标筛选候选。

默认 `norm_rank_marker = "NRR"` 按归一化应答率排序。其他选项包括 `RR`、`NRR+RR`、`RR+NRR`、`MeanNRR`、`LowerTailNRR`、`ShrunkLogOR` 和 `RankFusion`。定义及同分处理规则见 [subtypePerformance 参考文档](man/subtypePerformance.Rd)。

### 正确理解 AUC 与验证

遗留默认值 `norm_roc_method = "OptimizeRank"` 会利用评价标签重排等级。示例显式选择 `"Raw"` 以保留等级顺序，但 ROC 评价器仍会自动选择方向，因此 **仅使用 Raw 评价还不构成冻结的外部预测流程**。

外部验证应在训练数据上估计并冻结状态顺序、分界、映射及评分方向，外部结局不得参与模型或映射选择。训练目标改善不能证明外部预测收益或临床获益。

### 拟合与验证冻结映射

新增 `fitNormCCS()` 在搜索时直接使用固定方向 AUC：等级越高，越倾向预测 response=1。`predict()` 应用冻结映射，不读取结局；`evaluateNormCCS()` 报告独立来源表现和覆盖率。未见状态返回 `NA`。评价始终拒绝与训练相同的患者 ID，默认也拒绝训练中出现过的来源。

```r
# 复用上方合成 df；三个 cohort 在此代表三个独立合成来源。
frozen <- fitNormCCS(df[df$Cohort != "C", ], k = 3)
grades <- predict(frozen, df[df$Cohort == "C", ])
evaluateNormCCS(frozen, df[df$Cohort == "C", ])$summary
validated <- validateNormCCS(df, grid = data.frame(k = 2:3),
                            outer_folds = 3, inner_folds = 2, seed = 489)
validated$evaluation$by_source
```

`validateNormCCS()` 在按来源分组的每个内层训练折重新学习状态顺序和分界，不使用外层测试结局选型；返回最终部署模型、外层预测、调参证据与状态归属稳定性。可选的来源内结局置换会重新执行嵌套选型。排序收缩、等级人数及事件底线、来源支持和惩罚均可放入候选网格调节。示例值不代表经过临床验证的推荐设置。这些接口验证的是给定状态后的映射，上游模型与状态发现仍须独立隔离训练数据。来源信息、覆盖门槛、约束、选型规则和推断边界见[流程说明](docs/normCCS-generalization.md)。

## API 与文档

| 资源 | 内容 |
| --- | --- |
| [训练](man/ccs.Rd)与[数据准备](man/CCS-method.createCCSData.Rd) | 表达输入、基因集、注释及队列模型参数 |
| [亚型分类器](man/CCS-method.scaller.Rd)与[预测](man/CCS-method.predict.Rd) | 分类器后端与预测输入 |
| [应答与分级搜索](man/subtypePerformance.Rd)及[映射应用](man/newSubtype.Rd) | 输入列、目标函数、搜索模式及结果结构 |
| [冻结拟合](man/fitNormCCS.Rd)、[评价](man/evaluateNormCCS.Rd)与[嵌套验证](man/validateNormCCS.Rd) | 固定方向预测、训练内选型、来源支持与诊断 |
| [表示消融](man/ablation.Rd) | 冻结表示比较及参数 |
| [外置缓存说明](docs/external-cache.md) | 本地存储布局与缓存管理 |
| [变更记录](CHANGELOG.md) | 包更新与开发历史 |

安装后可在 R 中查询帮助，例如 `?CCS::ccs` 或 `?CCS::subtypePerformance`。[在线教程](https://huangwb8.github.io/ccs.principle/) 用于逐步讲解；完整、可直接使用的教程计划在论文正式发表后公开。

旧文档引用的部分历史分析目录已不在当前工作副本中。README 示例不附带大型数据集、训练模型库或生成报告；复现具体研究需要相应输入及分析环境。

## 开发与可复现性

在仓库根目录工作，可直接加载当前源码，无需安装：

```r
devtools::load_all(".")
```

项目协作规范见 [AGENTS.md](AGENTS.md)，另有 [Claude Code 适配文件](CLAUDE.md)。Windows 开发检查使用 R 4.3.1。场景测试位于 `test/`，需要在包检查之外单独执行；使用历史文档中的命令前，应先确认对应场景仍在仓库中。

对于已经定义 `targets` 依赖图的分析，执行、恢复和监控均须使用该分析的 `_targets.R` 及实际 store；target 级并行通过 `crew` 配置。大型数据、模型、包库及生成报告放在源码目录之外；本地默认采用 `D:/cache/ccs/<analysis>/`，既有分析继续沿用原缓存根目录。具体安排见外置缓存说明。

项目通过 [bensz-auto-contribution](https://github.com/huangwb8/bensz-auto-contribution) 将人类、AI 和工具贡献记录在 [docs/contribution.bac](docs/contribution.bac)。安装 `bac` 并具备 Python 3.10+ 后，可验证账本：

```bash
python -m bac --root . --bac-file docs/contribution.bac verify --json
```

请勿将患者标识、受保护健康信息或专有基因组数据提交到仓库或问题报告中。反馈问题时使用脱敏或合成示例。

## 引用、支持与许可证

CCS 论文目前正在审稿中，正式引用信息将在可用后补充。

遇到缺陷或使用问题，可提交 [GitHub issue](https://github.com/huangwb8/CCS/issues)，并附上 CCS 版本、R 会话信息和最小可复现示例。

**作者与维护者：**Weibin Huang（黄伟斌）· <hwb2012@qq.com> · [ORCID](https://orcid.org/0000-0003-0159-435X) · [博客](https://blognas.hwb0307.com/)。

CCS 使用 **Apache-2.0** 许可证，详见 [license.txt](license.txt)。依赖包与外部数据集分别遵循各自的许可证及访问条件。
