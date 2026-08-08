# CCS ablation-02 cohort-bank scaling 优化计划

## 问题是什么

- **当前现象：**`test/ablation-02` 已能比较完整 Direct-GSClassifier 与 Cohort-d1，也已计算癌种标签恢复、技术来源聚集、几何改写、module 方差贡献和输入重建等指标；最近新增的 cohort-bank scaling 却主要把不同 module 数量下的 cancer-type balanced accuracy 与 macro-AUROC 汇总成一条曲线。
- **核心偏差：**这条曲线没有直接回答 CCS 最重要的问题——新增独立 cohort evidence 是否让 d1 的技术稳健性、非冗余信息、跨 cohort 生物结构和状态稳定性出现可重复的 scaling。Cancer type 不是 CCS 的优化目标，因此癌种标签更难恢复只能作为 lineage recoverability 诊断，不能被解释成 CCS 的整体优势或劣势。
- **边界效应：**当前 d1 由完整的 150 个冻结 cohort modules 组成，现有分析只是从这套完整 bank 中截取嵌套子集。子集在 125–150 modules 时接近完整 d1，本来就容易出现收敛或平台；这只能描述固定 bank 内部的边际冗余，不能证明 cohort congress 在 150 附近达到一般性的最佳规模。
- **设计混杂：**现有 tissue-stratified sequence 会交错加入不同 tissue 的 module。早期 module 增加通常同时提高 tissue diversity，后期增加则更多发生在已覆盖 tissue 内，导致 tissue breadth 与 within-tissue cohort depth 被混进同一条 module-count 斜率。
- **比较口径漂移：**主消融使用 529 维完整 Direct-GSClassifier 输入，最近 scaling 默认只使用其中 496 个 ordinary gene-pair。后者仍由 GSClassifier 原生流程构造，但只是完整合同的子集，不能与既往完整 Direct 结果混为同一基线。
- **实现分叉：**`R/ablation.R` 中较早的 scaling 路径已经包含 effective rank、cohort mixing、biology purity 和 embedding stability，当前 representation scaling 则主要输出 cancer-type readout。两套内部逻辑并存，容易继续造成指标、命名和解释漂移。

如果不处理这些问题，报告即使统计计算正确，也可能把固定 150-module bank 的内部收敛误写成 scaling 规律，把 tissue composition 效应误写成 cohort 数量效应，并让非目标终点主导对 CCS 的评价。

## 要达到什么目标

- **完成后的变化：**分析主问题改为“cohort evidence 增加时，d1 的目标相关性质是否获得非冗余、可重复的改善，以及这种改善主要来自新增 tissue，还是来自既有 tissue 内新增 cohort”。
- **证据层级清楚：**技术边界、非冗余信息、跨 cohort 生物保持和状态稳定性进入主分析；几何差异与 cancer-type readout 作为机制或保护性诊断，不再承担 CCS 成败判断。
- **二维 scaling 可辨识：**分别估计 tissue breadth、within-tissue depth 及二者交互，并在相同 module 总数下直接比较 breadth-heavy 与 depth-heavy bank。
- **基线合同统一：**完整 Direct-GSClassifier 作为主基线；496 维 `gene_pair` 子集作为明确标记的敏感性分析。
- **计算与解释一致：**`.R` 保存全量 bank 设计、逐重复指标和审计信息；`.Rmd` 只负责阈值筛选、可视化和基于结果的专家解读。
- **不在本次处理范围：**不训练超出当前 150 modules 的新 cohort model，不据此宣称全局最优 bank 大小，不修改项目版本号，不把尚不存在的 pathway、TME、ICI 或 outcome 标签伪装成已验证终点。

本计划承接 `docs/plans/2026-08-05-ccs-d1-vs-tsp消融源码优化计划.md` 已完成的 Direct-vs-d1 基础工作，重点替换其中尚未正确落实的 cohort-scaling 部分；不会重新设计整个 CCS 下游流程。

## 改进方向

### 重新定义 scaling 的科学对象

每个候选 bank 应被视为一个明确的冻结 module 集合 $B$，由它生成对应的 $d1_B$。分析不再首先问 $d1_B$ 能否预测 cancer type，而是依次回答：

1. 新增 module 是否带来新的独立方向，而不是重复已有 module；
2. 技术来源边界是否减弱，同时跨 cohort 生物结构没有被无选择地抹平；
3. 同等规模、不同 tissue 构成的 bank 是否产生系统性差异；
4. 相同设计在不同随机 bank 组成和外部 query cohort 中是否复现。

指标按作用分层，避免把所有数值都解释成“越高越好”：

| 证据层 | 主要指标 | 允许的解释 |
|---|---|---|
| 非冗余信息 | normalized effective rank、module covariance redundancy、方差贡献集中度、边际独立贡献 | 新增 cohort 是否提供新的可利用方向 |
| 技术稳健性 | assay、platform、source 的 conditional neighbor excess | 技术来源聚集是否在生物或 lineage 背景校正后减弱 |
| 生物保持与外部效用 | 跨 cohort biology/pathway/TME/outcome 邻域一致性或 readout；无合格标签时明确 `not_evaluated` | 重表达是否保留或增强与 CCS 目标一致的外部结构 |
| 表征稳定性 | 同设计重复间的 neighbor Jaccard；固定参数下 d2/d3 状态的 ARI、split/merge 可作为次级扩展 | 结果是否依赖少数 module 或随机组成 |
| 机制诊断 | CKA、distance Spearman、Direct-vs-d1 kNN Jaccard、feature reconstruction | d1 怎样改写输入；方向本身不代表成功 |
| Lineage 诊断 | cancer-type balanced accuracy、macro-AUROC、top-k/MRR | 传统癌种信息的可恢复性，不评价 CCS 整体效用 |

### 构造 tissue breadth × cohort depth 二维设计

在正式计算前先生成可审计的 bank-design manifest。每一行记录设计类型、随机重复、module IDs、tissue 集合、每个 tissue 的 cohort 数、总 module 数、累计 d1 维度、query 覆盖关系和哈希。

二维设计包含三个互补部分：

- **Tissue breadth scaling：**固定每个 tissue 纳入的 cohort 数，逐步增加 tissue 数；回答新增组织背景是否带来跨组织的非冗余证据。
- **Within-tissue depth scaling：**固定 tissue 集合，逐步增加各 tissue 内的 cohort 数；回答同一组织内更多独立队列是否降低 cohort-specific 噪声、提高稳定性。
- **Matched-size contrast：**固定 module 总数，配对比较“较多 tissue、每个 tissue 较少 cohort”和“较少 tissue、每个 tissue 较多 cohort”；隔离 tissue diversity 本身，而不是把总维度变化误当成 diversity 效应。

真实 bank 中各 tissue 的 cohort 数不平衡，因此代码只构造实际可行的网格，不补造 module。每个网格单元使用多组随机但可复现的 tissue/module 组合；统计重复单位是 bank sequence 或 matched design pair，而不是同一 bank 内的样本或多个 module-count 点。

query 结果还要区分：bank 已覆盖 tissue 的 query、新增 tissue 对应 query、以及与新增 tissue 无直接匹配关系的 query。这样可以把“加入与 query 同 tissue 的模型所带来的覆盖收益”和“新增 tissue 对已有组织的跨组织 transfer”分开。

### 统一 Direct 合同与敏感性分析

主分析固定使用完整 529 维 Direct-GSClassifier 合同，即 32 个 `single_bin`、496 个 `gene_pair` 和 1 个 `set_pair`；所有特征继续由冻结模型的 `geneSet`、`breakVec`、feature manifest 和 GSClassifier 原生构造流程确定。

496 维 ordinary gene-pair 保留为 `Direct-GSClassifier-TSP` 敏感性分析，不能再简称为与主分析相同的 `Direct TSP`。两套 Direct 输入共享样本、query、预处理边界和审计信息，报告分别展示，不把两者的差值或斜率直接拼接比较。

对不依赖 Direct 的指标，例如 d1 effective rank、module redundancy 和同设计稳定性，不为了形式强行加入 Direct 差值。Direct 只在需要判断输入空间改写或 lineage recoverability 时出现。

### 收敛 `R/ablation.R` 的 scaling 实现

保留现有 GSClassifier 输入恢复、外部 d1 编码、module block 和 hash 审计逻辑；在此基础上将 bank 设计、指标计算和统计汇总分成三个清晰职责：

- **Bank design：**根据 module manifest 生成 breadth、depth 和 matched-size 设计，并验证嵌套性、组成平衡与可行性。
- **Bank scoring：**对每个 bank 子集构造 $d1_B$，复用已有的 CKA、距离相关、kNN、effective rank、technical excess、module diagnostics 和 stability 组件；对 module block 做一致的等权或标准化处理，避免输出列更宽的 module 自动占优。
- **Scaling summary：**输出 breadth slope、depth slope、交互、matched-size 配对差和边际增益；置信区间以 sequence/design repeat 为 bootstrap 单位。

较早的 `.ablation_experiment_scaling()` 与当前 `.ablation_representation_scaling()` 不再各自维护一套序列和指标规则。实施时应抽取共享内部 helper，由 representation 主入口调用；legacy 路径在不改变公共行为的前提下复用同一 helper，待项目负责人另行决定是否弃用。不要在本次计划中删除公共选项或引入新的 S4 类。

新的 `cohort_scaling.rds` 保留现有文件名以减少下游改动，但提升 schema version，并至少包含：

- `design`：每个 bank 的组成与审计字段；
- `metrics`：逐 bank、逐 query 分层、逐指标的长表；
- `contrasts`：breadth、depth、交互和 matched-size 配对结果；
- `diagnostics`：Direct 合同、lineage recoverability 与缓存资格；
- `status/reasons`：不可估指标和缺失外部锚点的明确原因。

缓存键应包含 bank-design hash、Direct 合同、query hash、指标配置和随机种子。Direct 基线只拟合一次；module-level d1 输出继续复用现有冻结编码缓存，避免为每个网格单元重复预测全部 modules。

### 重写 `test/ablation-02` 的报告主线

`02-ablation-experiment.R` 负责运行和保存所有设计单元的全量结果，`02-ablation-cohort-scaling.R` 继续作为只重算 scaling 的入口，二者必须调用同一个核心实现。`02-ablation-experiment_functions.R` 只保留报告整形和绘图 helper，不承载统计结论或偷偷应用业务阈值。

`02-ablation-experiment.Rmd` 的主线调整为：

- 数据、Direct/d1 合同与 bank 可行网格；
- d1 已观察到的表征改写与技术边界性质；
- tissue breadth scaling；
- within-tissue depth scaling；
- matched-size breadth-vs-depth 对照；
- lineage recoverability 与 pure-TSP 敏感性诊断；
- 讨论、限制与数字准确性验证。

当前 Figure 8 的单轴 `d1–TSP` 癌种差值曲线应被二维性质图和 matched-size 配对图替代。Cancer-type balanced accuracy 与 macro-AUROC 可进入诊断表或次级分面，但图轴和正文统一使用“cancer-type recoverability”，不再使用“d1 优势”“相对劣势”或“超过 Direct”等泛化措辞。

所有图继续输出矢量 PDF，并由 PDF 生成 JPG 预览后逐图检查。二维网格优先使用易读的分面曲线或热图；若真实可行网格稀疏，不用插值曲面制造不存在的数据。所有关键数字从 RDS 对象动态引用，正文不保留会随重跑漂移的硬编码结果。

### 同步方法学文档与审计边界

`docs/关于CCS框架的讨论.md` 中“cohort bank 扩张时 d1 优势如何演化”应改写为 cohort evidence scaling 的可证伪假设，删除把 125–150 modules 平台解释为普遍饱和的表述。新结果只能说明当前冻结 bank 范围内的 breadth/depth 性质，不能外推到新增模型或无限扩张。

`test/ablation-02/README.md` 更新执行顺序、结果文件和测试命令；`R/ablation.R` 的 roxygen 与生成的 `man/ablation.Rd` 同步新参数和输出合同；行为和分析逻辑变更记录到 `CHANGELOG.md` 的 `[Unreleased]`，不修改 `DESCRIPTION` 版本。

## 实施范围与顺序

1. 先冻结科学问题、指标角色、完整 Direct 主合同和 `cohort_scaling.rds` 新 schema，并用合成 module manifest 锁定二维设计不变量。
2. 在 `R/ablation.R` 中实现共享 bank-design 与 scoring helper，复用现有几何、technical excess、effective rank 和 stability 组件，消除两套 scaling 逻辑的规则分叉。
3. 将 representation 主入口和独立 scaling 脚本接到新实现，补齐缓存、manifest、query coverage 和失败原因审计。
4. 重写 Rmd 的 scaling 主线和图表，把 cancer-type 指标降为诊断，并同时呈现完整 Direct 主分析与 gene-pair 敏感性分析。
5. 运行合成合同测试、真实数据重算、报告渲染、视觉检查和包级验证，再同步讨论文档、README、roxygen/Rd 与 `[Unreleased]`。

主要涉及文件：

| 文件 | 计划中的职责 |
|---|---|
| `R/ablation.R` | 二维 bank 设计、共享指标引擎、summary、缓存与公开参数/结果合同 |
| `test/ablation-02/02-ablation-experiment_functions.R` | 实验参数与报告整形/绘图 helper |
| `test/ablation-02/02-ablation-experiment.R` | 完整分析编排与全量产物保存 |
| `test/ablation-02/02-ablation-cohort-scaling.R` | 复用既有输入、单独重算二维 scaling |
| `test/ablation-02/02-ablation-experiment.Rmd` | 图表、阈值、结果解释与数字追溯 |
| `test/ablation-02/tests/17-test-ablation-cohort-scaling.R` | 二维设计、配对、输出 schema 与可复现性合同 |
| `test/ablation-02/README.md` | 执行入口、产物与验证说明 |
| `docs/关于CCS框架的讨论.md` | 方法学结论与证据边界 |
| `man/ablation.Rd`、`CHANGELOG.md` | 生成文档与 `[Unreleased]` 记录 |

## 如何确认完成

### 科学与解释验收

- 报告的主问题明确指向 d1 的技术稳健性、非冗余信息、生物保持和稳定性 scaling，而不是 cancer-type 分类性能。
- 每个结果都能区分 tissue breadth、within-tissue depth 和 matched-size contrast；单一 module-count 斜率不再承担归因。
- 125–150 modules 的接近完整 bank 现象只描述为当前 bank 内部收敛，不宣称全局最佳规模。
- 完整 Direct-GSClassifier 是主合同；`Direct-GSClassifier-TSP` 只作为标记清楚的敏感性分析。
- 没有合格 pathway、TME、outcome 或其它外部锚点时，相应效用指标明确为 `not_evaluated`，不使用 cancer type 替代。
- CKA、距离相关和 kNN Jaccard 只解释为表征改写；cancer-type balanced accuracy、macro-AUROC 和 retrieval 只解释为 lineage recoverability。

### 合成合同测试

扩展 `17-test-ablation-cohort-scaling.R`，至少覆盖：

- breadth sequence 增加 tissue 数时，每个 tissue 的 cohort depth 保持固定；
- depth sequence 增加 cohort 时，tissue 集合保持固定；
- matched-size pair 的 module 总数相同但 tissue diversity 不同；
- cohort 数不平衡时只生成可行网格，并保存被排除设计及原因；
- 相同 seed 复现相同 module IDs、design hash、query hash 和统计结果；
- 同一设计内较大嵌套 bank 严格包含较小 bank；
- module block 宽度不同不会因列数自动改变等权距离贡献；
- full Direct 与 gene-pair sensitivity 的名称、特征数、hash 和结果字段不会混淆；
- bootstrap 以 sequence/design repeat 为单位，不能把同一 sequence 的多个网格点当成独立重复；
- cancer-type 指标在结果 schema 中标记为 diagnostic，缺失外部生物锚点时返回 `not_evaluated`。

### 运行与交付验证

实施完成后至少执行：

- `Rscript test/ablation-02/tests/17-test-ablation-cohort-scaling.R`
- `Rscript test/ablation-02/tests/11-test-ablation-representation.R`
- `Rscript test/ablation-02/tests/14-test-ablation-real-encoder.R`
- `Rscript test/ablation-02/tests/15-test-ablation-public-representation.R`
- `Rscript test/ablation-02/tests/16-test-ablation-public-readout.R`
- `Rscript test/ablation-02/02-ablation-cohort-scaling.R`
- `Rscript -e "devtools::document()"`
- `R CMD build .`
- 使用 `DESCRIPTION` 当前版本对应的 tarball 执行 `R CMD check`
- 使用 `bensz-rmd-rules` 的图表/表格解读、解释质量和 htmlwidget 可见性检查；通过 `knit-rmd-html` 渲染同名 HTML
- 对所有正式 PDF 的 JPG 预览逐图检查标签、字体、图例、颜色和画布比例
- `python -m bac --root . --bac-file docs/contribution.bac verify --json`

当前终端会话未能从 `PATH` 找到 `Rscript`。正式实施前应先恢复项目约定的 R 运行环境；在此之前不能把静态检查当成真实 R 结果验证。

## 风险与待确认事项

- **真实 bank 不平衡：**部分 tissue 可能没有足够 cohort 支撑完整 factorial grid。优先减少网格而不是重复使用 module 或填补虚假组合，并在报告中展示可行性边界。
- **Tissue coverage 混杂：**新增 tissue 可能只是加入了与 query 匹配的模型。必须按 query coverage 分层，才能讨论跨 tissue transfer。
- **维度与 module 宽度：**不同 bank 的 d1 列数必然变化。几何比较使用固定可行 rank或 block 等权距离；effective rank 同时报告绝对值、理论上限和归一化值。
- **计算量：**二维重复会显著增加组合数。先缓存完整 module 输出和 Direct 结果，再计算 bank 子集；通过预先声明的网格和重复数控制预算，不事后删除不利结果。
- **外部锚点不足：**当前 cancer type 和技术变量不能证明新状态有用。若没有独立样本级生物或临床锚点，本轮只能确认机制与稳定性 scaling。
- **150-module 上界：**本设计仍只能研究现有冻结 bank 内部的组成效应。若要判断超过 150 后是否继续改善，必须新增独立训练 cohort models，属于后续研究。
- **科学阈值：**最小重要差异应在查看新结果前由项目负责人确认；源码保存连续效应和不确定性，不把零差异或事后阈值硬编码为成功标准。
- **兼容性：**`cohort_scaling.rds` schema 会变化。通过版本字段、明确字段名和测试迁移降低风险，不静默让旧 Rmd 读取新结构。

