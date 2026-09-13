# ablation-03 输入说明

外部 query 的 d1 在 Stage 01 的数据入口一次性完成对象适配：从完整的
`PADv20240911/resCCS.rds` 读取所有样本行，按筛选版 resCCS 的冻结训练 bank
列合同合并到 `resCCS_ablation`。`ablation()` 始终只接收一个自洽 CCS 对象，
不写入 `external-d1-cache.rds`；没有预计算 query 行时，`ablation()` 发出警告并取已有 d1 样本交集，不在函数内部重编码。
# d1 输入边界（当前契约）

`ablation()` 不在下游分析中重算 query 的 d1。若 Stage 01 的表达数据包含
`object@Data$Probability$d1` 未覆盖的样本，函数会发出警告并只对已有 d1
样本取交集；被排除的样本 ID 会写入 `excluded-query-d1.csv`。因此，Stage 01
必须显式审计 CCS object 的 d1 覆盖范围，不能依赖函数内部的模型补算。
