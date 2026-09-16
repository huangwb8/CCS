# ablation-03 论文使用前静态审核

审查日期：2026-09-17。对象为当前工作树；审查开始时 Git 工作树干净。仅静态审查，不执行分析、模型训练或场景测试，不读取基因组 RDS，不修改分析源码。源码指纹与语法检查记录见 `.bensz-api/task-20260917-0054-ablation-03-review/`。

## 审查结论

**目前不能认定 ablation-03 已无问题，也不建议将现有全部结果直接用于文章定稿。** 主流程存在确定的 query 重复标准化错误，直接影响 readout 与 learning curve；anchor 的实际样本资格与声明不一致。另有缓存、替代入口和缺失样本配对问题。

以下区分“代码路径确定存在问题”和“是否影响现有数据需要核实”。静态审查能证明错误路径，不能证明磁盘上的旧结果一定由该路径生成，也不能估计修复后的效应大小、方向或显著性。未把条件风险写成已发生的数据污染。

## 范围与验证

重点追踪 `R/ablation.R` 中 representation 路径及其 Direct 重建、预处理、retrieval、readout、learning curve、breadth/depth scaling、decoder、缓存与输出；检查 ablation-03 的数据准备、R/Python biology、双向 structural、统计辅助函数、Rmd 消费逻辑、配置及 20–27 号测试。额外检查 ablation-02 的 readout 测试，解释测试覆盖的边界。旧 layered/metaccs 路径不属于本次结论范围。

实际验证：使用项目指定的 R 4.3.1，只调用 `parse()`，29 个 R 文件／Rmd 代码块集合解析通过；Python 入口通过 `ast.parse()`，未导入执行。R 启动出现四项 locale 设置警告，解析仍全部通过。没有执行 R CMD check、渲染报告、运行分析或运行现有测试；下文 Verify 均为修复验收建议，不是本轮已执行测试。

## 重要发现

### F01 · P1 · readout 和学习曲线的 query 被处理两次

- **位置**：`R/ablation.R:5568–5586`、`5643–5654`、`6076`；变换定义在 `1627–1637`、`4151–4185`、`5869–5876`。
- **What / 证据**：主流程从 `transformed$direct$query`、`transformed$d1$query` 取得 test，却传入原始 `prepared$reference_direct/d1` 作为 train。`.ablation_linear_readout()` 在最终拟合前再次对 train/test 调用 `.ablation_readout_transform()`。学习曲线也接收这两份已变换的 test。
- **Why**：Direct 应使用 `(x−μ)/σ`，实际 query 变成 `(((x−μ)/σ)−μ)/σ`；d1 还重复施加 module 权重。训练和预测坐标不同，且两臂失真不同。子集学习曲线还混用全 reference 与子集尺度。不能据此判断分类性能差距或“更多数据无法弥补差距”。
- **Fix**：向 readout/learning curve 传入同一资格行对应的原始 query；保留函数内部的折内和最终预处理。不要简单删除内部变换，否则破坏 inner CV。
- **Verify**：以非零均值、非单位方差、多个 block 的合成矩阵，验证主入口传到预测器的 test 等于一次正确变换；同时覆盖完整 reference 与训练子集。随后重算 readout、learning curve 及相关图文。现有 `12-test-ablation-readout.R` 直接传原始矩阵给 helper，无法发现主入口错误。

### F02 · P1 · anchor 仍被癌种资格筛选截断

- **位置**：`R/ablation.R:5419–5427`；`01-ablation03-biology-cache.R:45–47`；`03-ablation-biology.R:27–31`。
- **What**：资格表将 anchor 作为无需癌种支持的 endpoint，但唯一保存的 `retrieval.rds` 邻居来自 `cancer_retrieval` 的 estimable view。biology 缓存和评分都只消费这份邻居表，没有对 anchor view 单独检索。
- **Why**：无 reference 癌种支持或支持不足的 candidate 被排除，尽管被宣称可做连续 anchor。实际分析总体比资格表、README 所称范围窄。结构模块单独使用完整 candidates，不会自动补齐局部 anchor 分析。
- **Fix**：为完整 anchor candidates 输出不依赖癌种支持的邻居表，并让缓存和 biology 使用它；若确实只分析癌种可估计子集，则需更改资格表和论文总体定义。
- **Verify**：构造有/无 reference 癌种支持的两个 query cohort，检查癌种 endpoint 排除后者，而 anchor 的邻居、分数与推断仍保留后者。

### F03 · P1 · anchor 的 reference 尺度不是完整、稳定的 reference 合同

- **位置**：`01-ablation03-biology-cache.R:47`、`87–95`；`03-ablation-biology.R:63–83`。
- **What**：缓存只保留 query 与两臂 top-15 命中的 reference 样本；所谓 reference global scale 随后只在这个被命中的 reference 子集拟合。增删 query 或更换检索表示，即使 reference atlas 不变，也可能改变尺度。同时缓存用原始 atlas 的 `tissue/cohort` 生成 key，而排除 external 使用经过 tissue 恢复的 manifest key。
- **Why**：尺度依赖待评价的邻居选择，utility 不再使用独立冻结的 reference 标尺。若 atlas 仍使用 `Undefined/C` 而 manifest 使用 `T/C`，该 external cohort 会错误落入 reference 拟合集合。后一种污染是否在当前外部 atlas 发生，未读取数据确认；代码没有防止它。
- **Fix**：根据审计 metadata 的 sample/cohort 身份确定 reference，缓存固定 reference 拟合样本及完整候选 query，统一 cohort key；断言 fit 与 query 样本不相交。结构模块现有的 `cohort_key_lookup` 可作为映射思路。
- **Verify**：固定 reference、增删 query 后 reference 均值/标准差保持不变；加入原始 `Undefined/C`、规范化 `T/C` 的案例，确认 external 样本不进入拟合。

### F04 · P1 · 缺失 anchor 时，所谓配对差可能比较不同 query

- **位置**：`03-ablation-biology.R:115–123`；`03-ablation-biology_functions.R:31–58`。
- **What**：评分用 inner merge 丢弃无分数的邻居；若某臂一个 query 全部邻居缺分数，这个 query 从该臂消失。统计函数先分别求各臂 cohort 均值，再按 anchor/cohort 合并，而非先按 query 配对。
- **Why**：同一 cohort 内 Direct 均值可能基于 A、B，d1 均值仅基于 A；二者相减包含样本组成差异。`query_count` 后来计算交集，但实际 estimate 未限制到该交集。当前若没有不对称缺失则不触发，静态审查没有确认缺失率。
- **Fix**：先按 anchor、query_cohort、query_sample 合并两臂，计算 query 差，再聚合 cohort。同步冻结有效邻居数规则并报告未配对 query。
- **Verify**：令 B 仅在一臂有分数且 utility 极端，确认加入/移除 B 不改变完整配对估计；检查分母与报告 query_count 一致。

### F05 · P1 · Python 备用入口与 R 方法不同，且会混写结果

- **位置**：`03-ablation-biology.py` 的 `_zscore_score()`、`main()`、输出部分；`README.md:16`；对照 `03-ablation-biology.R:63–101`。
- **What**：Python 每个 cohort 自行求均值和总体标准差，R 则拟合 reference gene-wise 均值和样本标准差。Python 没有生成 cohort-level `anchor_inference.csv`，却覆盖同目录的 coverage、utility、contrasts CSV。
- **Why**：两入口估计的不是同一指标。若先运行 R、再运行 Python，旧 R 推断文件会留在目录内；Rmd 仅检查文件存在，就可能把 Python 效应量和旧 R 的 CI/P 值展示在一起。README 当前将 Python 列为可替代入口。
- **Fix**：统一计算、尺度、缺失处理、配对推断和输出合同，并给结果绑定 run/method hash；或明确停用 Python 的正式产物输出资格。
- **Verify**：同一小型非同尺度数据上比较两入口中间分数、效应和推断；测试先 R 后 Python、失败中断和缺文件场景，禁止混合方法结果。

### F06 · P1 · breadth/depth scaling 使用历史组织名，未接入组织恢复

- **位置**：`R/ablation.R:964–987`、`7200–7208`、`6641–6650`；对照 `04-ablation03-structural-reproducibility_functions.R:5–35`。
- **What**：module manifest 的 tissue 直接从 d1 列名取出。scaling 直接据此分组，没有应用 stage-01 的 tissue_resolution_audit；只有结构模块显式将 `Undefined` 映射回真实 tissue。
- **Why**：当 bank 保留历史 `Undefined` modules 时，来自不同真实组织的 cohort 会被视为同一个 tissue，breadth、depth、matched-size 与 query coverage 随之失真。stage-01 恢复表达矩阵与 metadata 并不会修改 d1 的历史列名。
- **Fix**：保留 module_id 的冻结标识，另提供审计后的组织字段用于分析设计；遇到无法映射的 Undefined 不应继续声称真实组织 breadth/depth。
- **Verify**：两个历史 Undefined modules 分别映射到不同组织，检查 tissue_count、depth 和覆盖率与手算相同；确认映射不改变 d1 数值或列顺序。当前 bank 是否包含此类模块应通过汇总审计核实。

### F07 · P1 · scaling 拟合缓存缺少关键决定因素

- **位置**：`R/ablation.R:4901–4910`、`7257–7274`、`7282–7295`；另见 `5054–5080` 的 legacy geometry promotion。
- **What**：`prepared$cache_key` 绑定 query Direct/d1 和 module IDs；scaling 在此基础上加入样本 ID hash、设计、geometry、scaling，但未绑定 reference 数值、train/test 标签、validation nrounds 等实际拟合参数。Direct contract 的 feature hash 只代表特征名字。
- **Why**：保持样本 ID、query 和特征名不变，仅改 reference 表达/d1、标签或 nrounds，就可能复用旧 fit，与新计算的 geometry/retrieval 拼接。legacy geometry promotion 也仅核对旧 query key、维度与参数，不验证 reference 内容。
- **Fix**：缓存键涵盖所有决定预测与评分的输入内容、metadata 和拟合配置；升级 schema。旧 geometry 没有足够 reference 指纹时重算，不将其提升为可信缓存。
- **Verify**：逐一改变 reference 数值、标签、validation nrounds，要求缓存 miss；完全相同输入才 hit。现有 Direct cache 测试不覆盖 scaling fit cache。

### F08 · P1 · 报告与增量重算没有跨阶段来源一致性门禁

- **位置**：`02-ablation03-experiment.Rmd:79–112`、`130–163`；`02-ablation03-cohort-scaling.R:56–80`；`04-ablation03-structural-reproducibility.R:93–119`、`382–390`。
- **What**：报告独立读取 manifest、各 endpoint RDS、biology CSV、structural CSV/RDS，主要检查存在性。单独 scaling 重算会将新 config 写回旧主 manifest，却不验证旧 retrieval/readout 的配置身份。结构缓存键沿用 biology 缓存中记载的旧 source md5，缓存重建读取 atlas 时也不重新核验该文件；`full_d1_hash` 只哈希维度和行列名。
- **Why**：局部重算、输入原位更新或中断后，报告仍可能正常渲染并显示混合批次结果。hash 字段的存在不等于其证明了内容一致性。这里确认的是校验缺失，不断言当前磁盘结果已混批。
- **Fix**：为每阶段结果记录并验证实际输入内容、算法版本、配置和父产物 hash；报告入口验证依赖链。partial scaling 更新保留自己的配置身份，不能替旧 endpoint 改写 provenance。结构入口验证当前 atlas hash，d1 provenance 绑定内容。
- **Verify**：替换一个结果文件为其他 run、原位改 atlas 或 d1 数值、只重算 scaling，报告应明确拒绝不一致组合。

### F09 · P2 · 常量/大量并列 anchor 会被强行划成高低状态

- **位置**：`04-ablation03-structural-reproducibility_functions.R:420–452`。
- **What**：按 score、sample_id 排序后直接取固定数量两端，没有验证 score 变异，也没有处理切分点的并列值。
- **Why**：一个 cohort 的 anchor 分数全相等时，仍生成 low/high 两个实体，只由样本 ID 决定分组；达到 min_entity_n 后会继续进入几何分析。这些实体没有实际高低生物含义。是否触发需检查分数分布。
- **Fix**：常量 anchor 标记 not_estimable；明确 ties 策略，检查高低区间可分离，记录排除原因。
- **Verify**：全常量与阈值大量 ties 的合成分数不能产生伪高低实体；更换样本 ID 不应改变有意义的状态定义。

### F10 · P2 · signature 配置与缓存有效性未真正冻结

- **位置**：`01-ablation03-biology-cache.R:50–78`；`03-ablation-biology.R:46–55`；`config/biological-anchors.yml`。
- **What**：YAML 只做 token 存在检查，IFN/IL6 实际取 signature 列表的第 1、2 个元素而非声明的名字。缓存记录 signature hash，但 R 消费端只检查 signature 路径，没有重算/比较 signature 或配置内容。
- **Why**：signature 列表重排或同路径内容/配置更新，可能选错基因或继续使用旧锚点，同时产物名称仍显示当前 anchor。这个缺陷不证明当前基因集已选错。
- **Fix**：按配置中的 family/name 精确选择，校验非空基因集；明确 cache 是冻结快照还是要求同步当前源，两种模式都应有可核验内容身份并写入报告。
- **Verify**：列表重排保持结果不变；按名缺失时报错；同路径改 signature/config 后按声明规则失效或明确报告旧快照版本。

## 统计解释与论文表述边界

### Technical excess 不能单独证明去除了技术噪声

`R/ablation.R:4385–4409` 的 expected rate 来自同癌种 reference pool，而 observed rate 来自未限定癌种的近邻。两臂使用同一个 expected，故配对差满足 `(O_d1−E)−(O_direct−E)=O_d1−O_direct`；扣减 baseline 并不消除两臂癌种检索构成差异。报告讨论已部分承认混杂，但摘要和部分解释仍使用“削弱技术来源捷径”等强措辞。

可报告“同技术来源邻居比例/该 excess 指标下降”；若要支持独立技术去偏，应补同癌种条件化、匹配或其他能区分 lineage 与技术构成的分析。三个技术字段本身可能高度相关，不能等同三个独立机制证据。

### 反向投影的 Direct 并非全流程重新在 external bank 冻结

结构入口有意继续使用 150-module Direct feature contract（`04...R:45–47`），仅交换 scale 的 fit/target；Direct 特征支持与 breakVec 则由 reference models 恢复（`R/ablation.R:1007–1090`）。反向 target 是这些 reference cohorts，因此“目标未参与当前 external d1 modules 训练”不等于“目标未参与所有比较表示的特征确定”。

将反向结果描述为固定 Direct 基准下的互惠性/设计敏感性检查较准确；若宣称两臂对称的全流程独立外部验证，应另行冻结不会使用目标队列的 Direct 特征合同。本轮不判定既有 bank 的实际特征选择训练流程是否进一步引入其他泄漏。

### CI 跨零不是等价或非劣效证明

结构模块的节点 bootstrap 比把所有 cohort-pair 当独立样本更合理。但 d1−Direct CI 覆盖零仅说明未确定差异，不能仅凭此证明“结构已保留”或非劣效。可结合绝对相似度做明确标注的描述性陈述；若要确认性地证明保留，应预先定义有科学意义的容许损失界值并评估对应区间。

### 低维 decoder 的误差不证明表示不可逆

`R/ablation.R:7592–7623` 限制 PCA rank 并使用 ridge linear decoder。非零误差既可能来自表示丢失，也可能来自 decoder 容量/正则化。论文宜表述为“本 decoder 下的可恢复性”，不能把它直接等同于 d1 的信息论不可逆性。

## 已有设计中可以保留的部分

- retrieval 的 MRR@k 明确将首个匹配位于 k 之外的 query 记为零；已有对应测试源码。
- 主 readout helper 内部按 cohort 分组 CV，并在训练折拟合预处理；问题在调用方传入的 query 尺度，不需要推翻整个 helper。
- 通用配对推断以 cohort 均值作为单位，小样本使用 exact sign flip，大样本 Monte Carlo 使用加一修正。sign-flip 的“精确”仍以零假设下符号交换/对称性成立为前提。
- learning inference 对完全相同的训练 cohort 设计去重，不将 100% reference 的多次拟合当多个独立设计。
- 结构分析按相同实体交集比较两表示，并用 cohort-node bootstrap 处理共享节点；matched-bank 范围明确不是 CI。
- Direct 特征缓存比早期路径更完整地绑定表达与特征合同；这不能替代其他缓存或最终报告的来源检查。

## 论文使用建议

| 结果模块 | 当前建议 |
|---|---|
| 主 readout、learning curve | 修复 F01 后重算，暂不引用当前性能差、CI/P 值或由此推导的机制 |
| 局部 biological anchor | 修复/澄清 F02–F05、F10，核验来源后重算 |
| breadth/depth、matched-size scaling | 核实并修复 F06/F07，再生成统一来源的结果 |
| retrieval、native geometry、decoder | 未发现 F01 直接传入这些路径；核验 F07/F08 与输入来源后，可按具体指标做描述性报告，仍不等于本轮确认数值正确 |
| 双向结构复现 | 核验 F08/F09 与反向 Direct 合同；保留探索性/描述性解释，不能直接宣称等价、非劣效或普遍生物学增益 |
| 现有 HTML 与出版图 | 应由修复后同批次结果重新生成，避免沿用旧推断与硬编码方向性解释 |

本轮仅交付审查结论，没有实施修复。建议先处理 F01 和 anchor 数据流，再处理 scaling 与跨阶段 provenance，最后使用有针对性的集成回归验证并统一重算。这些工作完成前，不给出“已经可以放心用于文章”的结论。

## 后续入口清理

本审查完成后，已移除旧的 `03-ablation-biology.py` 入口。ablation-03 的 biological-anchor 分析现在仅通过 `03-ablation-biology.R` 执行；本报告中对 Python 实现的描述保留为审查时的历史证据，不代表当前仍存在第二套分析入口。
