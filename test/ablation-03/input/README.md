# ablation-03 输入说明

外部 query 的 d1 在 Stage 01 的数据入口一次性完成对象适配：从完整的
`PADv20240911/resCCS.rds` 读取所有样本行，按筛选版 resCCS 的冻结训练 bank
列合同合并到 `resCCS_ablation`。`ablation()` 始终只接收一个自洽 CCS 对象，
不写入 `external-d1-cache.rds`；没有预计算 query 行的普通调用仍走通用重编码路径。
