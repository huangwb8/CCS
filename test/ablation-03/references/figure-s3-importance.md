# Figure S3 cohort Gain 的历史来源

此目录的 CSV 是此前已保存结果的只读摘要导入，用于表示报告中的指标对照讨论，不是 ablation-03 新拟合的 importance 结果。

- 图稿：论文 `Figures.pptx` 第 15 页 Figure S3 的 cohort 面板。
- 幻灯片备注对应的原始图文件：`plotImportance_merge_5ff3a2de76e6cf902e765e8224f9cb66.pdf`。
- 已保存对象：`resCCS_5ff3a2de76e6cf902e765e8224f9cb66.rds`。对象内容 MD5、样本数、模块数与 d1 列数见 `figure-s3-provenance.csv`。
- 导入内容：`figure-s3-cohort-gain.csv` 保留整个历史 bank 的模块 ID、cohort ID、两个坐标的 Gain 与 merge Gain，不含患者 ID、表达数据或患者级概率矩阵。
- 历史 d3 列名为 `all|D1`、`all|D2`，图面显示为 allD1、allD2；它们指两个最终坐标，不是 d1/d2 表示层。merge 为两坐标归一化 Gain 的等权平均。
- 历史 importance 表中未出现的 cohort 对应树模型未使用的输入，汇总记为零。导入时已核验 cohort ID 唯一、各坐标 Gain 合计为 1，以及 merge 与两轴均值一致。

计算含义可核对项目 `R/importance.R` 的 `ctImportance()` 与 `mergeImportance()`。本记录定位历史输出并说明汇总定义，不保证历史软件版本、患者范围或概率值与当前分析完全一致。

三个文件均纳入 `representation_report_sources`。报告对比时从 CSV 读取历史 Gain，从正式 `native_geometry` 读取当前方差份额；这些份额不作配对效应量、相减、相除或因果解释。此处未增加模型训练、外部任务验证、bootstrap 区间或 P 值。
