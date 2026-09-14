# ablation-03 biology 数据审计整合实施计划

## 通俗解释：究竟发生了什么

当前 `01-ablation03-test-data.Rmd` 只展示 32 特征主分析数据画像，biology anchor 则由另一条脚本从完整表达矩阵读取后生成缓存，导致同一消融实验的数据边界分散在两份报告/脚本中。目标是把 biology anchor 的输入覆盖和基因签名审计并入这份 Rmd，同时避免在 Rmd 中重复加载大矩阵。

## 专业判断：问题在哪里

- 现有 Rmd 只读取 `data-profile.rds`，无法说明 biology anchor 实际使用了多少基因、多少 cohort 和多少样本。
- biology 缓存已经保存了完整表达矩阵的来源哈希、所需基因集合、cohort 样本和每个 anchor 的覆盖率；应复用这些可追溯结果，而不是在展示层重算。
- biology 缓存只覆盖 retrieval 邻居涉及的样本/cohort，不等于全量 39,112 个样本；报告必须明确这个范围。

## 要达到什么目标

- 在同一份 Rmd 中增加 biology-anchor 数据审计章节。
- 动态报告缓存状态、完整表达来源、所需基因数、覆盖样本/cohort 数和四类 anchor 的基因覆盖率。
- 展示每个 anchor/cohort 的覆盖审计表，并用文字解释它与 32 特征主数据画像的关系。
- 缓存缺失或不完整时快速失败，避免生成不完整的“合并审计报告”。

## 改进方向

### 复用缓存并集中审计摘要

在 `_functions.R` 增加轻量摘要函数，Rmd 只负责读取缓存、呈现表格和解释结果。摘要不保留完整表达矩阵，只使用缓存中的 metadata、coverage 和 sample key。

### 明确两条数据边界

在报告正文中区分 32 特征主分析画像、完整表达矩阵来源，以及 biology 阶段提取的 316 个唯一 anchor 基因；同时注明 biology 缓存的样本范围来自 top-15 retrieval 邻居。

## 实施范围与顺序

1. 扩展 `test/ablation-03/01-ablation03-test-data_functions.R`，提供缓存校验和 anchor 覆盖摘要。
2. 扩展 `test/ablation-03/01-ablation03-test-data.Rmd` 的 setup、数据审计章节、关键函数表、讨论和数字准确性表。
3. 使用现有 biology cache knit/执行级检查，确认动态数字、覆盖表和缺失缓存错误路径均可追溯。

## 如何确认完成

- Rmd setup 能读取 `tmp/ablation-biology/expression-anchor-cache.rds`，并拒绝缺失或非 complete cache。
- 报告同时出现主数据画像和 biology-anchor 审计，且 biology 样本数不被误报为全量样本数。
- `Rscript` 语法解析、相关脚本检查和 `bensz-rmd-rules` 图表/表格覆盖检查通过。
- BAC verify 无链错误。

## 风险与待确认事项

- biology cache 是外部表达矩阵的派生缓存；若源矩阵或 signature 发生变化，缓存校验应提示重建。
- 本次不修改 biology 计算逻辑、不把完整矩阵写入 Rmd，也不自动重建缓存。
