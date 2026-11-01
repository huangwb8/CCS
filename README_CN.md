<div align="center">
  <h1>CCS · Cohort Congress System</h1>
  <p>
    <a href="https://github.com/huangwb8/CCS"><img src="https://img.shields.io/badge/version-0.8.7-brightgreen.svg" alt="版本 0.8.7"></a>
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

依赖包括 `xgboost`、`torch` 和 `luz`。深度学习流程还需要 torch 运行时；使用该后端时，请参阅 [R torch 安装指南](https://torch.mlverse.org/docs/articles/installation.html)。安装失败时，应先确认具体失败的依赖及其仓库或系统要求，再重试。Windows 下若需从源码编译，还可能需要与 R 版本匹配的 Rtools。

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

连续模式先固定状态顺序。K = 1 为恒定评分，AUC = 0.5；后续每个子方案保留父方案的全部分界，再新增一个分界。所有合法的单次细分按所选目标评分后，每层最多保留 `norm_top_n` 个候选。该模式固定使用单核继承式 DP，不受 `norm_search_method` 设置影响。

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

## API 与文档

| 资源 | 内容 |
| --- | --- |
| [训练](man/ccs.Rd)与[数据准备](man/CCS-method.createCCSData.Rd) | 表达输入、基因集、注释及队列模型参数 |
| [亚型分类器](man/CCS-method.scaller.Rd)与[预测](man/CCS-method.predict.Rd) | 分类器后端与预测输入 |
| [应答与分级搜索](man/subtypePerformance.Rd)及[映射应用](man/newSubtype.Rd) | 输入列、目标函数、搜索模式及结果结构 |
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
