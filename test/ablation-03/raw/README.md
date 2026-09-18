# ablation-03 外部只读输入

本分析不复制基因组原始数据到仓库。输入由 `01.00.00. 数据准备.R` 通过 `CCS_DATA_ROOT`、`CCS_SYNC_ROOT`、`CCS_FULL_RESCCS_RDS`、`CCS_ABLATION_MODEL_ROOT`、`CCS_FULL_EXPRESSION_RDS` 和 `CCS_GENE_SIGNATURE_RDS` 定位，并在读取前做显式路径检查。原始数据与冻结模型只读；所有派生缓存写入 `CCS_ABLATION_CACHE_ROOT`。

