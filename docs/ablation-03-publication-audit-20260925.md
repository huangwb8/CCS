# ablation-03 出版前指标、推断与图表复核

## 审查范围与证据

本次审查沿正式 `test/ablation-03/_targets.R` 的输入、计算、补充推断、四份 R Markdown 报告和 PDF/JPG 预览追踪结果；数值取自 `D:/cache/ccs/_ablation-03` 的正式 targets 产物。没有改动原始输入或 CCS 包版本。图表字体参照 `bensz-rmd-rules` 的静态图建议及 PanPAD `04.06` 系列；检查全部图表参数，并目视核对关键图的标签、颜色、尺度和裁剪。

| 终点 | 实际指标与独立单位 | 本次核对 |
| --- | --- | --- |
| 癌种检索 | Top-k 同癌种邻居比例、MRR@k、top-1；外部 query cohort | MRR 在首次命中排位超过 k 时记零；三项同族 Holm 校正。单 cohort 点无总体显著性含义。 |
| 癌种 readout | 合并 query 的 balanced accuracy 与队列等权 accuracy 是不同估计对象；外部 query cohort | 报告与推断表已分开命名，未将单癌种 cohort accuracy 误称 balanced accuracy。 |
| 技术近邻 | d1−Direct 同技术邻居比例差；外部 query cohort | 两臂扣除相同背景，背景项在差值中严格抵消。图的纵轴已改为实际差值，不能解释为选择性去除技术噪声。 |
| 生物锚点 | 同一 query 的连续表达邻居效用差；外部 query cohort | 四锚点采用 cohort 等权差与配对检验；覆盖率只衡量可计算性，效用不等于临床效益。 |
| 结构复现 | 同一 cohort 对的距离排序 Spearman 差；共享 cohort 节点 | 不能将队列对当独立重复；稀疏同癌种网络的节点检验保留 NA。区间跨零不等于等价。 |
| Decoder | gene-pair balanced accuracy/Brier、连续特征 Spearman/MAE/RMSE；外部 query cohort | 合并 query 后逐特征计算再平均的列已改名；RMSE 基准差先在每个 cohort、特征上开方，再汇总。 |
| 表示几何 | CKA、距离 Spearman、kNN Jaccard；reference cohort | CKA/距离的区间评估队列组成；kNN 固定完整近邻图，只重采样 query cohort。三者均无预设自然零假设，不填造 P 值。 |

## 已纠正的问题

1. Decoder 补充表原列名 `sample_weighted_estimate` 与计算顺序不符。合并 query 的 gene-pair balanced accuracy 为 0.8363，队列等权估计为 0.7864（95% CI 0.7742–0.7970）；两者差异有实际解释价值。现列名为 `pooled_query_estimate`，正文说明其逐特征计算方式，区间只对应 `cohort_equal_estimate`。
2. Decoder 基准比较的原文字把 RMSE 写得像逐 query 损失差；实际是每个 cohort、特征先对平方误差均值开方，再比较并汇总。已改正计算说明，未改变真实计算值。
3. 技术近邻图原纵轴写作 “same-technology excess”，但其两臂差严格等于实际同技术近邻比例差。已改为 “same-technology neighbor rate”，避免将背景扣除误读为混杂控制。
4. kNN Jaccard 原区间重算时删去未抽中的 cohort，改变近邻候选库，得到 0.3356–0.3795 的区间；原完整图点估计 0.3349 不在区间内。现在固定完整图、整组重采样 query cohort，区间为 0.3132–0.3547；代码、方法标签、报告解释和合成回归测试均已同步。
5. 结构热图的完整横轴队列名互相遮挡。两张图共用相同行列顺序和色阶；横轴改用顺序说明，队列名称只显示一次，正文指向明细表以定位具体队列对。

## 出版解释边界

四类表达锚点的队列等权效用差均为负，约 −0.044 至 −0.055，BH 调整后 P 值均约 0.0005。这是当前连续表达近邻任务上的局部损失。正反向结构平均差分别约 0.0019、0.0091，节点区间均跨零；未设可接受损失界值，不能从中宣布结构等价或无损。100% 学习曲线只有一种独立训练设计，设计层 CI/P 仍为 NA；固定训练设计下的外部 query cohort 区间是另一估计对象。

Reference d1 含 in-sample 生成路径，query d1 来自冻结外部预测；表达锚点又来自相关表达数据。当前结果适合报告为冻结表示的诊断性比较，不能据此声称个性化泛癌分类已有独立临床验证、技术噪声已被选择性消除，或生物结构已证明无损。

## 验证记录

- `C:/R/R-4.3.1/bin/Rscript.exe` 运行 `tests/40-test-geometry-bootstrap.R` 与 `tests/42-test-statistical-inference.R`；前者另以独立计算复核固定近邻图的队列区间。
- 正式输入 `D:/cache/ccs/_ablation-03/formal-inputs.rds`、正式 store `D:/cache/ccs/_ablation-03/targets`，由 `scripts/run-targets-renv.ps1 -Action make` 调用 `targets::tar_make()` 增量生成统计产物及四份 HTML。
- 静态图保存于 `test/ablation-03/reports/figures/`；HTML 位于 `test/ablation-03/`。正式运行仍提示既有 renv lock 与环境不同步；跨机器复现前应按项目环境说明补齐包来源。
