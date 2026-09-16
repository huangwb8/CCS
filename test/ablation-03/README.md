# ablation-03：表示性能、生物学锚点与结构可复现性

这是同一个消融实验的三个分析模块。文件前缀表达依赖顺序：`01*` 准备输入，`02/03/04` 分析已有缓存，同名 `.Rmd` 展示该模块的结果并生成 HTML。

本次重构仅调整计算边界、文件组织与报告归属，保留原有样本选择、Direct/d1 构造、统计方法、参数与随机种子。当前处于**代码审核阶段**：旧实验已停止，尚未执行新流程，也未生成新 HTML。审核通过后再运行下面的命令。

## 文件与职责

| 文件 | 职责 | 主要产物 |
|---|---|---|
| `00.Environment.R` | 环境、包与路径初始化 | 无实验计算 |
| `00-workflow_functions.R` | 缓存读取与递归来源核验 | 缺失或失效时停止，不自动重建 |
| `01a-ablation03-prepare-data.R` | 读取原始输入，恢复组织标签，审计样本与冻结 CCS 对象 | `tmp/01-data/inputs.rds`、`data-profile.rds` |
| `01b-ablation03-prepare-representations.R` | 按原参数准备 Direct、已有 d1 与完整 reference/query 契约 | `tmp/01-representations/` |
| `01c-ablation03-prepare-biology.R` | 按样本契约准备连续 anchor 与结构分析的表达子集 | `tmp/01-biology/` |
| `01-ablation03-data-overview.Rmd` | 汇总三个准备步骤的数据审计 | 同名 HTML |
| `02-ablation03-representation.R` / `.Rmd` | 几何、癌种检索、readout、learning curve、scaling、technical excess、decoder | `tmp/ablation-experiment/`、同名 HTML |
| `03-ablation03-biology.R` / `.Rmd` | 连续表达锚点的邻居效用与配对推断 | `tmp/ablation-biology/`、同名 HTML |
| `04-ablation03-structural-reproducibility.R` / `.Rmd` | 双向生物状态结构复现与 matched-bank 敏感性 | `tmp/ablation-structural-reproducibility/`、同名 HTML |
| `02-ablation03-cohort-scaling.R` | 可选：只重算表示模块的 scaling | 更新表示结果与来源凭据 |

各模块的 `_functions.R` 保留原有计算或报告辅助函数。旧名称入口和旧总报告已移除，不提供兼容入口。图表仍写入 `figures/`，JPG 预览写入 `tmp/plot-previews/`；保留原图号，方便对照既有科研材料。

## 执行顺序

在 RStudio 中打开 `ablation-03.Rproj`，或从仓库根目录执行脚本；不需要手动 `setwd()`。以下命令以仓库根目录为工作目录，Windows 使用 `C:/R/R-4.3.1/bin/Rscript.exe`。

```powershell
# 审核通过后，先准备所有输入。
& C:/R/R-4.3.1/bin/Rscript.exe --vanilla test/ablation-03/01a-ablation03-prepare-data.R
& C:/R/R-4.3.1/bin/Rscript.exe --vanilla test/ablation-03/01b-ablation03-prepare-representations.R
& C:/R/R-4.3.1/bin/Rscript.exe --vanilla test/ablation-03/01c-ablation03-prepare-biology.R

# 再按模块分析；这些入口不执行原始数据准备，也不补建输入缓存。
& C:/R/R-4.3.1/bin/Rscript.exe --vanilla test/ablation-03/02-ablation03-representation.R
& C:/R/R-4.3.1/bin/Rscript.exe --vanilla test/ablation-03/03-ablation03-biology.R
& C:/R/R-4.3.1/bin/Rscript.exe --vanilla test/ablation-03/04-ablation03-structural-reproducibility.R
```

`01a → 01b → 01c` 完成后即可渲染数据概览，无需先跑实验。`03` 使用 `02` 的邻居检索结果；`04` 使用 `02` 的输入一致性凭据及 `01c` 的结构表达缓存，不依赖 `03` 的统计结果。`02` 报告只要求表示分析完成，`03/04` 报告各自要求对应模块完成。

报告渲染也在审核之后执行。例如：

```r
rmarkdown::render("test/ablation-03/01-ablation03-data-overview.Rmd")
rmarkdown::render("test/ablation-03/02-ablation03-representation.Rmd")
rmarkdown::render("test/ablation-03/03-ablation03-biology.Rmd")
rmarkdown::render("test/ablation-03/04-ablation03-structural-reproducibility.Rmd")
```

Rmd 读取完成的 RDS/CSV，沿用原有报告摘要、区间汇总与绘图代码，不运行实验入口。每份报告生成同名 HTML；这里不把旧总报告改名冒充模块报告。

## 输入与分析边界

外部路径继续使用 `CCS_DATA_ROOT`、`CCS_SYNC_ROOT`、`CCS_FULL_RESCCS_RDS`、`CCS_ABLATION_MODEL_ROOT`、`CCS_FULL_EXPRESSION_RDS` 和 `CCS_GENE_SIGNATURE_RDS` 等原有环境变量。`CCS_ABLATION_CORES` 在准备阶段确定并随参数保存。分析阶段修改环境变量不会悄悄重建输入；需要新输入时应显式重跑相应准备步骤。

`01b` 分别保留表示分析的 seed `20260805` 与结构输入准备的 seed `20260912`。`sample-contract.rds` 在准备阶段生成，生物学表达缓存因此不再要求主实验先完成。结构分析的表达提取前移到 `01c`；尺度拟合、状态构造和结构比较仍由 `04` 按原方法执行。

癌种检索、technical excess、readout 与 learning curve 使用预先声明的可估计子集；连续 anchor 与结构分析保留全部候选外部 cohort，实际有效覆盖另行审计。readout 两臂均以原始表示进入函数，由每个训练折拟合尺度；d1 继续直接使用冻结 CCS 对象已有的概率矩阵。

连续 anchor 以完整 reference 拟合共同基因尺度，要求两臂各有完整 top-15 分数，先按 query 配对、再按 cohort 推断。结构比较保留双向 module-bank、分位边界并列处理、cohort-node bootstrap 及 matched-bank 设计；matched-bank 的范围仍是设计敏感性，不是置信区间。

## 缓存与重算

本次编辑没有删除、移动、读取或重写现有实验缓存，也没有启动重算。被中止的旧进程可能留下未完成产物，不能仅凭文件存在认定某个阶段完成。

- `01a/01b/01c` 新增的准备阶段产物与来源凭据，需要审核后显式生成。
- Direct 的内容键缓存仍位于 `tmp/ablation-experiment/direct-feature-cache.rds`。`01b` 使用原缓存校验函数，只有内容键一致才复用；失配时才重建。
- 精确几何缓存与 scaling fit 缓存保持原路径及校验算法。是否复用取决于实际完整性与内容键，不能在静态审核阶段保证全部命中。
- 三个分析阶段的旧来源凭据会因脚本调整失效。审核后需要通过新入口生成完整结果和凭据，不能直接修改哈希“认可”旧结果；计算中已有严格内容键缓存仍可按原规则命中。
- 不需要重新训练 cohort bank 或重新预测原始 d1。本次也没有加入自动全量重跑入口。

只重跑 scaling 会更新表示阶段凭据，下游报告会拒绝混用旧批次结果；此时需按依赖更新生物学与结构结果。普通报告改版不应触发实验计算。

## 静态审核

此次交付只执行 R/Rmd 语法解析、源代码引用检查与迁移前后计算表达式比较；不加载真实基因组数据，不执行实验测试，不渲染图表。审核时重点看：

- `R/ablation.R` 将原表示计算体抽为 `.ablation_run_prepared_representation()`，公共 `ablation()` 的原始输入路径仍保留。
- `01b/01c` 是否完整保存下游所需的输入，`02/03/04` 是否只读取准备好的输入。
- 三份分析 Rmd 是否分别包含原报告中所属模块的全部图表、解释与参数。

历史审核材料反映重构前的文件布局，当前使用方式以本文件为准。
