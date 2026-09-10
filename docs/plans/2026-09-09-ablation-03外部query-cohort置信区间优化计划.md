# ablation-03 外部 query-cohort 置信区间优化计划

## 通俗解释：究竟发生了什么

- **一句话说明：** 当前完整 reference bank 已经在 17 个外部 cohort、5705 个样本上完成预测，数据足以估计“换一个外部 cohort 时，d1 与 Direct 的平均差异有多不确定”，但不能把同一训练集合下的 10 次重复拟合当成 10 个独立训练实验。
- **具体场景：** 这类似于用同一套仪器到 17 家医院做检测。每家医院的患者数不同，但若研究问题是“换一家医院后结果是否仍成立”，统计上应让 17 家医院各有一票，而不是让患者最多的医院决定总体结论。
- **对应到本问题：** 完整 reference bank 对应固定仪器，17 个 query cohorts 对应 17 家医院，每个 cohort 内的 accuracy 差对应该医院中 d1 相对 Direct 的表现变化；cohort bootstrap 则模拟重新抽取一组外部医院。
- **改变前后：** 当前 Figure 4 的 100% 点只能报告点估计，设计层 CI 不可估计；优化后仍保留这一边界，同时在独立的 full-reference 外部泛化结果中明确报告 cohort 等权差、95% CI 和 P 值。当前结果约为差值 `-0.0757`、95% CI `[-0.2209, 0.0427]`、P = `0.321`，应解释为“观察到下降趋势，但外部 cohort 层面的证据仍不确定”。

## 专业判断：问题在哪里

- **当前现象：** `learning_curve.rds` 的 100% fraction 有 10 次模型运行，但只有 1 个唯一 training-cohort subset，因此不能估计训练 cohort 组成变化带来的 CI。与此同时，`readout.rds` 已保存 17 个外部 query cohorts 的配对预测、cohort-level accuracy 和样本数，具备估计固定完整 reference bank 外部泛化 CI 的必要信息。
- **指标边界：** 每个 query cohort 只含一个 cancer type，因此 cohort-level endpoint 必须称为 accuracy，而不能称为 balanced accuracy 或 macro-AUROC。跨 cohort 的主要 estimand 是 17 个 `accuracy_d1 - accuracy_direct` 的等权平均。
- **有效统计单位：** 主要单位是 query cohort，不是 5705 个样本，也不是 10 次算法运行。样本用于计算各 cohort 内的 accuracy；进入总体推断后，每个合格 cohort 等权。
- **当前证据：** 17 个 cohort 已足够生成并报告探索性的 cluster-bootstrap 95% CI，但区间较宽且跨零。该结果能够支持“方向为负但不确定”，不能支持“对所有新 cohort 普遍下降”或“已证明无差异”。
- **主要限制：** cohort 效应异质性较大，样本量从 72 到 1094 不等；多个 cohort 共享 TCGATARGET/CPTAC 等来源，不能假定 17 个 cohort 完全代表 17 个相互独立的数据生成体系。

## 要达到什么目标

- **完成后的变化：** 报告把“训练设计不确定性”和“固定完整 bank 后的外部 cohort 泛化不确定性”严格分开；读者能在主表、Figure 3/Figure 10 与正文中看到同一套 cohort 等权效应量、95% CI、P 值、样本数和解释边界。
- **主要统计对象：** 固定完整 reference bank 与全部既定模型参数后，随机抽取一个符合当前纳入条件的外部 query cohort，其 d1 相对 Direct 的 accuracy 差的平均值。
- **预期论文结论：** “在 17 个外部 query cohorts（5705 个样本）中，d1 相对 Direct 的 cohort 等权 accuracy 差呈负方向，但 95% CI 跨零，因此观察到下降趋势，尚不能确认其在外部 cohort 层面具有普遍性。”
- **不在本次处理范围：** 不重新定义 CCS 的成功标准，不把 cancer-type accuracy 当作新分型或临床效用，不新增外部数据，不用患者级 bootstrap 替代 cohort-level 推断，也不修改项目版本号。

## 改进方向

### 固定主要 estimand 与分析人群

以当前已冻结的完整 reference bank 为条件，纳入现有 17 个通过资格检查的 query cohorts。每个 cohort 分别计算 d1 和 Direct 的 accuracy，并先在 cohort 内配对作差，再对 cohort 等权平均。纳入规则、样本配对、真实标签和两臂预测必须在看到效应量之前固定；缺失预测或两臂样本不一致时早期失败，不做静默删样本。

这意味着论文回答的是“换一个外部 cohort 后平均会怎样”，而不是“随机抽一个患者会怎样”。

### 用 query-cohort bootstrap 估计 95% CI

从 17 个 cohort 中有放回地抽取 17 个 cohort，每次计算 cohort 等权平均差，重复至少 10,000 次，以百分位数法给出双侧 95% CI。固定随机种子并记录 bootstrap 次数，使 HTML、测试和后续复算得到一致结果。现有 2,000 次可用于探索，但出版版本提高到至少 10,000 次，以降低区间端点的 Monte Carlo 抖动。

同时保留精确 paired sign-flip P 值作为“平均差是否偏离零”的辅助证据。CI 是主要不确定性表达；P 值不替代效应量，也不因 CI 跨零而写成“证明没有差异”。若该 endpoint 被预先指定为单独的 full-reference readout，P 值无需与学习曲线各 fraction 混成同一多重比较家族；若与 retrieval endpoints 共同宣称，则继续明确报告相应的 Holm 校正结果。

### 分离 Figure 4 与 full-reference 泛化结果

Figure 4 继续展示训练比例变化：10%、25% 和 50% 的阴影表示不同 training-cohort designs 之间的不确定性；100% 因只有 1 个唯一训练设计而不画设计层 CI，并标注 `not_estimable`。

固定完整 bank 的 query-cohort 95% CI 放在 Figure 3 的 Overall 行、Figure 10 的 `Readout cohort accuracy` 行及配套推断表中。正文增加明确交叉引用，避免读者把这一区间误认为“100% training fraction 的训练设计 CI”。如需在 Figure 4 附近提示，只写“full-bank external cohort CI 见 Figure 3/10”，不在同一 ribbon 中混合两种统计单位。

### 增加异质性与影响诊断

保留 17 个 cohort 的配对差散点，并增加 leave-one-cohort-out 诊断：每次去掉一个 cohort 后重算总体均值，报告估计范围及影响最大的 cohort。该诊断用于说明结论是否由单一 cohort 主导，不作为筛除异常 cohort 的依据。

同时并列保留样本加权 pooled accuracy 差作为描述性结果，但明确标注其回答的是“当前 5705 个样本的总体差异”，不能替代 cohort 等权的外部泛化估计。对于共享数据来源造成的依赖，只做按 `source_system` 的描述性分层和 leave-one-source-out 敏感性检查；来源层级过少时不生成看似精确的 source-level CI。

### 统一报告语言和可追溯字段

所有表格和图注统一使用以下名称：

- `Full-reference external query-cohort accuracy difference`
- 方向统一为 `Cohort-d1 minus Direct-GSClassifier`
- 单位统一为 `query cohort`
- 明确列出 `n_cohort = 17`、`n_query_sample = 5705`、bootstrap 次数、随机种子、CI 方法、P 值方法和推断状态

正文将“显著下降”“稳定下降”等强表述替换为与结果匹配的句子：点估计为负、95% CI 跨零、外部 cohort 层面的方向尚不确定。不得把 CI 跨零改写为“两种表示等价”。

## 实施范围与顺序

1. 在现有 `readout.rds` 上完成预测配对、cohort 资格、17 个 cohort 与 5705 个唯一 query samples 的合同检查，确认无需重新训练模型。
2. 将现有 readout inference helper 固化为 full-reference query-cohort 推断入口，把 bootstrap 次数提升至出版配置，并补充方法、种子和状态字段。
3. 新增 leave-one-cohort-out 与来源分层诊断，确保它们只作为敏感性结果，不改变主要 estimand 或事后删除 cohort。
4. 同步 Figure 3、Figure 10、推断表、Figure 4 交叉引用和讨论段落，确保两类 CI 不混用。
5. 运行回归测试、Rmd 静态门禁、HTML 渲染和十张图预览检查，并通过 BAC 记录文件变化和验证证据。

## 涉及范围

- `test/ablation-03/02-ablation03-experiment_functions.R`：full-reference query-cohort CI、影响诊断和报告字段。
- `test/ablation-03/02-ablation03-experiment.Rmd`：主表、Figure 3、Figure 4 交叉引用、Figure 10 与论文级解释。
- `test/ablation-03/tests/23-test-publication-inference-contracts.R`：补强现有统计合同；如职责过多，可新增相邻编号的独立测试脚本。
- `test/ablation-03/02-ablation03-experiment.html` 与 `test/ablation-03/figures/`：重新渲染的正式交付物。
- `CHANGELOG.md`：实施完成后将行为与报告变化记录在 `[Unreleased]`，不修改 `DESCRIPTION` 版本号。

核心模型训练代码原则上不需要修改，因为 full-reference predictions 和 cohort-level paired results 已存在。只有合同检查发现 `readout.rds` 缺少必要预测字段或两臂样本无法一一配对时，才回到 `R/ablation.R` 补充输出。

## 如何确认完成

- 主要结果严格复现 17 个 query cohorts 和 5705 个唯一 query samples，两臂的 `sample_id`、cohort 与真实标签完全一致。
- cohort 等权估计直接等于 17 个 cohort 内 accuracy 差的算术平均，不能由样本加权结果替代。
- 10,000 次以上 cluster bootstrap 在固定种子下可复现；改变输入行顺序不改变估计和 CI。
- 合成测试证明：当 cohort 样本量极不均衡时，主要结果仍对 cohort 等权；只有 1 个有效 cohort 时返回 `not_estimable`。
- 精确 sign-flip P 值使用完整枚举尾概率，不应用 Monte Carlo `+1` 校正；若未来 cohort 数超过精确枚举门槛，才切换到带校正的 Monte Carlo 方法并明确标记。
- leave-one-cohort-out 输出包含 17 次估计，并能识别影响最大的 cohort；该结果不触发自动剔除。
- Figure 3 与 Figure 10 显示相同的点估计和 95% CI；Figure 4 的 100% 点仍无 training-design CI，并能明确跳转到 full-bank external cohort 结果。
- 最终 HTML 成功渲染，图中无裁切、标签碰撞或过小字号；图表解释门禁、定向回归测试、`git diff --check` 和 BAC 哈希链验证通过。

## 风险与待确认事项

- **外推范围有限：** 17 个 cohort 足以报告当前探索性 CI，但共享来源和癌种构成限制了对所有未来医院、平台或癌种的外推。报告必须保留这一边界。
- **单一高影响 cohort：** CPTAC-PAAD-2021 的 accuracy 差约为 `-0.936`，可能明显影响均值和 CI。应通过 leave-one-out 如实展示影响，不应仅因结果极端而删除。
- **endpoint 的研究地位：** cancer-type accuracy 是 lineage 诊断，不是 CCS 新状态效用的最终指标。即使未来 CI 完全低于零，也只能确认旧癌种标签可恢复性下降，不能单独证明 CCS 成功或失败。
- **分析属性：** 本次优化发生在已有结果之后，若没有事前注册，应在论文中标记为 exploratory 或 post hoc inference，避免使用确认性措辞。

## 计划完成标准

计划实施后，读者应能同时得到三个不混淆的答案：

- 100% training fraction 的点估计仍为负，但训练设计 CI 当前不可估计；
- 固定完整 reference bank 后，外部 query-cohort 的平均 accuracy 差可以报告 95% CI；
- 当前 CI 跨零，因此结论是“观察到下降趋势但证据不确定”，而不是“显著下降”或“证明无差异”。
