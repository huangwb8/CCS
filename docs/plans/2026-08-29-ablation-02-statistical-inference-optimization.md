# ablation-02 统计推断增强实施计划

## 通俗解释：现在缺的是什么

当前结果已经显示 Direct-GSClassifier 与 Cohort-d1 之间存在方向一致的差异：d1 的癌种标签检索、MRR 和线性 readout 都较低。但这些数字目前主要是点估计，不能回答“这种差异在重新抽取 cohort 或 query 后是否仍然稳定”“差异是否可能只是抽样波动”。

本计划不改 CCS 包的基础表示、模型或 S4 接口，而是在 `test/ablation-02` 的实验汇总层增加不确定性估计和配对显著性检验。最终报告同时给出效应量、95% CI、p-value、统计单位和多重比较校正方法；如果某项分析没有足够的独立重复，就明确标记为探索性或 `not_estimable`，不强行生成显著性结论。

## 专业判断：统计单位与主要问题

- 样本属于 cohort，不能把所有 query 当成相互独立的重复。主要推断单位采用 cohort；样本级差异保留为描述性结果。
- 检索的主要终点预先指定为 `k = 30` 的 `top_k_label_rate` 配对差；MRR、`k = 5/15` 和 `top1` 作为次要终点。
- readout 的主要终点预先指定为外部 query 的 balanced accuracy 配对差；macro-AUROC 和各 query cohort 的差异作为次要或探索性终点。
- learning curve 的主要比较点为使用 100% reference cohorts 时的 readout 差；其余训练比例用于判断样本效率，不把每个比例都当作独立主要结论。
- technical-neighbor excess、module scaling 和 cohort-specific 极端值不与主要癌种终点混合解释；它们单独报告并进行相应的多重比较控制。

## 要达到什么目标

完成后，报告应能清楚区分三种说法：

1. “观察到 d1 比 Direct 低多少”（效应量）；
2. “在合理的 cohort/query 重采样下，差异范围是多少”（95% CI）；
3. “在预先指定的零假设下，差异是否有统计学证据”（配对 p-value）。

统计推断只用于判断差异的稳定性，不把 p-value 当作 CCS 新状态有效性的替代终点。

## 改进方向

### 增加统一的 cohort-level 推断 helper

在 `test/ablation-02/02-ablation-experiment_functions.R` 增加分析层 helper，输入现有 `retrieval.rds`、`readout.rds` 或 learning-curve 对象，输出统一的 tidy inference 表。helper 至少包含：

- 点估计、`ci_low`、`ci_high`、`p_value`；
- `n_cohort`、`n_sample`、重采样次数和随机种子；
- `unit`（`cohort`、`query` 或 `design_repeat`）；
- 检验方法、备择方向、零假设和 `multiplicity_method`；
- `status` 与 `reason`，用于记录重复不足、标签不合格或无法估计的情况。

所有重采样都必须按 cohort 整块抽取，并在同一抽样索引下配对 Direct 与 d1。

### 外部癌种标签检索

对每个 query 先保留 d1−Direct 的配对差，再进行两层汇总：

- 描述性层：保留现有样本加权总体均值和 cohort 均值，方便与旧报告对照；
- 推断层：以 cohort 为 cluster 做 stratified/bootstrap 重采样，得到 `top_k_label_rate`、MRR 和 top1 差异的 95% CI。

显著性检验采用 cohort-level paired sign-flip 或 permutation test，零假设为平均配对差等于零。若 query cohort 数量允许，再按癌种或技术来源做预先指定的分层敏感性分析；逐 cohort 的 p-value 只作为探索性结果，并使用 Benjamini–Hochberg 或 Holm 校正。

检索结果表和 Figure 2 应增加总体 CI、p-value、有效 cohort 数，并在图中用零线和误差线展示差异；不再只展示一个没有不确定性的总体点。

### 外部癌种线性 readout

当前 readout 是一次 reference 拟合后在外部 query 评分。统计增强分两层：

- 固定预测的主要推断：对已生成的外部 query 预测按 cohort 做 paired bootstrap，估计 balanced accuracy 差异的 CI，并对 query cohort 做配对 sign-flip/permutation 检验；这回答“在当前冻结模型下，外部 query 的差异是否稳定”。
- 端到端敏感性分析：重复抽取 reference cohort，重新运行已有 `.ablation_linear_readout()`，再在相同 query 上计算 d1−Direct；这部分评估训练 cohort 组成变化带来的不确定性，结果单独标记为 `refit_bootstrap`，不与固定预测 CI 混淆。

readout 的 grouped CV 重复次数从当前 3 次提高到至少 10–20 次；100% reference fraction 的差异作为主要检验，其余 learning fractions 报告 CI 并进行 Holm 校正。若计算预算不足，必须保留固定预测层，并把端到端结果标记为未完成，而不是用 3 次重复宣称显著性。

### 学习曲线与 cohort-bank scaling

- learning curve 以同一 cohort subset、同一 query 和同一随机重复进行配对；对每个 fraction 输出差异 CI、配对 p-value 和调整后 p-value。
- scaling 的现有 bootstrap 只以 5 个 design repeats 为单位，先保留为探索性结果；若要进入确认性统计，增加独立 bank-design repeats（建议至少 15–20 个），再对 breadth/depth slope 与 matched-size contrast 计算 CI 和 sign-flip p-value。
- scaling 的 p-value 不与癌种 readout 的 p-value 合并，不用 scaling 显著性替代独立生物效用。

### 多重比较与报告口径

在 Rmd 的 setup 中集中声明检验清单和校正规则：主要终点使用双侧 alpha = 0.05；同一终点下多个 k、多个 learning fraction 或多个技术因子使用 Holm 校正；逐 cohort 结果使用 BH-FDR 并标记为探索性。报告同时显示未校正 p-value、校正方法和调整后 p-value，避免读者误把多个局部检验当成一个总体检验。

正文使用“估计差异为 X，95% CI 为 [L, U]，p = …”的格式。若 CI 跨零或 p-value 未达到预先阈值，使用“证据不足以拒绝零差异”，不写成“没有差异”。

## 实施范围与顺序

1. 审计现有 `retrieval.rds`、`readout.rds`、`learning_curve.rds` schema，锁定主要终点、cluster 单位、重采样种子和结果字段；不改变 `R/ablation.R` 的基础函数。
2. 在 `02-ablation-experiment_functions.R` 实现 cohort bootstrap、paired sign-flip/permutation、CI 计算、多重校正和失败状态记录。
3. 扩展 `02-ablation-experiment.Rmd` 的 setup、总体结果表、Figure 2/3/4 和讨论文字；所有数字从 inference 对象动态引用。
4. 增加 `test/ablation-02/tests/` 下的统计合同测试，覆盖已知正差异、零差异、cohort 聚类、缺失标签、重复不足和多重校正。
5. 先用小型合成对象验证方向与 p-value，再按计算预算重算真实结果；若端到端 refit 过慢，保留固定预测推断作为主层并在报告中说明边界。
6. 重新渲染 `02-ablation-experiment.html`，检查表格、误差线、p-value、CI、图例和正文动态数字。

## 需要保留的边界

- 不修改 `R/ablation.R`、`DESCRIPTION`、版本号或模型文件；统计代码只服务于本次 ablation-02 分析。
- 不把样本量很大的 query 产生的极小 p-value 误解为生物学效应很大；正文必须同时报告效应量和 CI。
- 不把 reference 内部样本几何当作独立外部验证；readout 与 retrieval 的外部 query 边界继续保持。
- 只有 3 次或 5 次独立重复的现有结果只能作为 pilot/exploratory；增加重复后才能把相应 p-value 作为确认性证据。
- 癌种标签仍是 lineage recoverability 诊断，不能替代 pathway、TME、outcome 或新 pan-cancer state 的外部效用终点。

## 如何确认完成

- 统计结果表对主要终点均包含 estimate、95% CI、p-value、调整后 p-value、统计单位和有效重复数。
- 合成测试能在零差异和已知方向下给出正确的 CI 方向、配对关系和多重校正结果。
- 真实报告中 Direct/d1 使用相同 query、相同 cohort 重采样索引和可追溯随机种子；不存在把样本或同一设计点重复当作独立 cluster 的情况。
- Figure 2/3/4 与正文数字动态一致，CI 和 p-value 不依赖手工硬编码；HTML 成功重新渲染。
- `git diff --check`、相关 `Rscript test/ablation-02/tests/*.R`、报告渲染和 `python -m bac --root . --bac-file docs/contribution.bac verify --json` 均通过。

## 风险与待确认事项

- 当前 query cohort 数和 scaling repeat 数决定 p-value 的分辨率；如果独立 cluster 太少，应报告精确置换 p-value 或 `not_estimable`，而不是套用正态近似。
- 端到端 readout refit 可能增加显著计算量；需要在真实重算前记录预计重复数、缓存策略和失败恢复方式。
- 多个癌种与 cohort/技术来源绑定时，统计显著不等于因果独立；条件化/匹配分析仍是必要敏感性分析。
- bootstrap CI、sign-flip p-value 和 grouped-CV 重复的含义不同，报告必须分别标注，不能把它们写成同一种“95% 置信区间”。
