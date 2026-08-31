# ablation-03 初始缓存说明

外部 query 的 d1 不再通过 `external-d1-cache.rds` 重算或缓存。03 直接从完整的
`PADv20240911/resCCS.rds` 读取外部样本行，并按筛选版 resCCS 的冻结训练 bank
列对齐；筛选版对象继续负责 reference 边界、filtered cohort 合同和下游分析配置。

若缓存校验失败，`R/ablation.R` 会自动重新编码并覆盖 03 自己的缓存。
