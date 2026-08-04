# 关于 CCS 框架的讨论

本文凝练围绕 CCS（Cohort Congress System）框架的一系列方法学讨论，重点回答三个问题：为什么从 TSP 出发构建 CCS、metaCCS/normCCS 是否属于必要步骤，以及面向免疫检查点抑制剂（immune checkpoint inhibitor，ICI）疗效预测时应如何选择建模路线。

本文是一份研究思路与决策记录，不代表相关假设已经获得实验验证。文中的“稳定”“可迁移”“泛癌”等主张，均须通过预设的跨队列实验加以检验。

## CCS 的原始问题意识

泛癌转录组模型的临床转化首先面临的不是模型复杂度，而是跨平台、跨中心和跨标本类型的可迁移性。RNA 绝对表达量容易受到测量平台、建库流程、批次、标本类型和标准化方法影响。若这一地基不稳，即使模型在单一队列中具有较高精度，也很难成为可复用的单样本临床工具。

CCS 最初的目标可以概括为：

> 利用大量异质转录组队列，构建一种可对单个新样本独立计算、能够跨平台迁移，并可服务不同临床任务的泛癌转录组证据层。

这一目标对于小众癌种尤其重要。单一罕见癌种通常缺少足够样本训练稳定模型，而泛癌、多队列信息有机会提供可迁移的先验与证据。

## 证据层预训练语料的最小数据要求

第一阶段的数据采集目标不是建立字段齐全的临床数据库，而是用现实中能够稳定获得的最少信息，支持 CCS 学习癌种相关证据、识别队列异质性并检验跨平台迁移。临床结局、分期、治疗、肿瘤纯度等信息在公共队列中经常缺失，不应作为预训练样本的纳入条件。“可服务不同临床任务”指冻结的证据层以后可以由带有相应标签的下游数据读取，不意味着预训练语料必须包含这些任务标签。

预训练语料由表达矩阵和一张最小样本信息表组成。样本信息表只要求以下字段：

| 变量 | 定义 | 在 CCS 中的用途 |
|---|---|---|
| `sample_id` | 项目内唯一且与表达矩阵中的样本一一对应的标识；来源数据的样本名发生冲突时可加 `cohort_id` 前缀 | 作为表达数据与样本信息之间的连接键，防止跨队列合并时发生样本错配；它本身不作为生物学或技术特征进入模型 |
| `cancer_type` | 样本所属的项目级主癌种，采用统一名称或编码；组织学类型和分子亚型不混入该字段 | 提供第一阶段必需的生物学身份，用于组织不同癌种的 cohort evidence、描述预训练语料的癌种覆盖，并检查证据层保留的是癌种相关结构还是技术来源；它是否直接作为训练标签由具体预训练方法决定 |
| `cohort_id` | 预先定义的最小独立建模与评估单元；若同一公共数据集包含来源清楚且彼此独立的子队列，应分别编码 | 标记数据生成环境的主要边界，用于按 cohort 划分训练与验证，并检验模型是否只记住某一数据来源；它不能替代患者标识，也不能单独排除跨数据集重复患者造成的泄漏 |
| `assay_type` | 检测技术的大类，第一阶段至少区分 `RNA-seq` 与 `microarray` | 定义跨平台评估中的技术域，用于分层统计、平台留出验证和识别模型是否依赖某一种测量体系；不需要作为单样本评分的输入 |
| `platform_id` | 更具体的测量平台标识，例如 microarray 的 GPL 编号或芯片型号、RNA-seq 的测序仪平台或型号；建库和分析流程不混入该字段 | 描述不同芯片的探针覆盖和其他平台内差异，使“跨平台”结论能够追溯到实际测量系统；若托管者未提供，可记为 `unknown`，但未知平台样本不能用于支持具体 platform-level 的迁移结论 |

表达矩阵是另一项必需资产，而不是样本信息表中的变量。其一维必须能通过 `sample_id` 与样本信息表对齐，另一维的 probe 或 gene 标识必须能够映射到 CCS 使用的统一基因空间，否则无法构造可复现的 TSP，也不应纳入正式预训练。原始队列可以保留 counts、TPM、FPKM 或芯片标准化信号等不同数值形式；第一阶段不要求把它们作为绝对表达量直接拼接，但必须按队列锁定实际使用的数值形式。counts、TPM 与 FPKM 之间包含基因特异变换，不能假定它们产生完全相同的 TSP，后续必须通过同队列敏感性分析检验关键基因对是否受预处理形式影响。数据来源已经提供的数值类型和预处理说明应随队列保留；若托管者未提供，可以标记为未知，不必自动推断，但相应的预处理稳健性主张也应降级。

这套最小结构只保留生物学身份、队列与平台来源以及数据连接所必需的信息：`cancer_type` 回答样本是什么肿瘤，`cohort_id` 回答样本来自哪组研究环境，`assay_type` 与 `platform_id` 回答表达谱由什么测量体系产生，`sample_id` 保证所有信息可以可靠连接。其他元数据有则保留、缺失则接受；只有当队列明确混合肿瘤与正常样本，或同一患者存在多个重复样本时，才应按具体风险补充样本类型或患者标识，而不把它们扩大为所有预训练数据的统一硬门槛。

上述字段定义的是预训练语料的最小数据合同，不同时决定预训练目标。cohort submodel 究竟通过何种无监督、自监督或弱监督目标形成证据，仍须在建模方案中单独定义；不能仅凭这些字段宣称某一种 cohort model 已经可以被训练。这里能够确定的是：临床任务标签不属于公共预训练语料的共同准入条件，后续任务应使用各自可获得的标签读取或检验冻结证据层。

## TSP 为什么适合作为基础表示

对于样本 $i$ 和基因对 $(a,b)$，TSP 特征可写为：

$$
x_{i,(a,b)}=I(g_{ia}>g_{ib}).
$$

TSP 将绝对表达量转换为样本内部的相对排序，主要优势包括：

- 不依赖表达量的绝对尺度；
- 对样本内共同施加的严格递增变换保持排序；基因特异、探针特异或平台特异变换仍可能改变排序；
- 比绝对表达量更容易跨 RNA-seq、microarray 和不同标准化流程迁移；
- 每个特征都能解释为具体的基因排序关系；
- 单个新样本可以独立计算，不要求与训练队列共同标准化。

这些优势不意味着 TSP 对批次效应完全免疫。RNA 质量、细胞组成、探针设计、基因覆盖和系统性测量偏差仍可能改变基因排序。更准确的表述是：TSP 提供了一种具有较强跨平台潜力的低依赖表示，而不是自动消除全部技术差异。

TSP 还具有明确的信息代价。绝对表达差异和变化幅度被压缩为二元关系；共享基因的 TSP、互补 TSP 以及满足排序传递性的 TSP 之间也高度相关。因此，TSP 解决了输入 $x$ 的可比性，但没有自动解决以下问题：

- 什么才是具有生物学或临床意义的目标 $y$；
- 全局泛癌聚类是否主要反映组织来源和细胞组成；
- 高维、相关的 TSP 如何避免距离集中和重复计权；
- 如何证明发现的结构能够在未见队列中复现。

## 纯 TSP 矩阵的聚类选择

当样本表示为 $0/1$ TSP 矩阵时，常见聚类方法包括：

| 方法 | 适用特点 | 主要注意点 |
|---|---|---|
| Hamming 距离 + PAM/层次聚类 | 0 和 1 表示排序的两个对称方向 | 适合作为最直接、可解释的基线 |
| K-modes | 用众数表示典型二元模式 | 需要指定簇数，对初始化敏感 |
| Bernoulli 混合模型 | 输出每个簇中各 TSP 取 1 的概率 | 条件独立假设可能被相关 TSP 破坏 |
| 谱聚类或图社区发现 | 适合非球形和复杂邻域结构 | 依赖相似度、邻居数和分辨率参数 |
| 二元双聚类 | 同时发现患者亚群与 TSP 模块 | 适合亚型只由部分基因对定义的情况 |
| 共识聚类 | 通过样本与特征重采样评价稳定性 | 是稳定性框架，而不是独立的数据模型 |

普通 TSP 中的 0 和 1 通常只是两个排序方向，而不是“事件缺失”和“事件发生”，因此 Hamming 或 Simple Matching 通常比 Jaccard 更自然。只有当 1 被特别定义为稀有、有方向意义的事件时，才应优先考虑 Jaccard。

## CCS 相对于纯 TSP 的理论价值

CCS 与 TSP 可以形成两层异质性控制：

$$
\text{表达矩阵}
\rightarrow \text{样本内 TSP}
\rightarrow \text{队列子模型}
\rightarrow \text{跨队列证据整合}.
$$

第一层由 TSP 处理样本内测量尺度问题；第二层由 CCS 处理队列之间的分布差异。不同队列子模型可被理解为从不同组织、平台和临床背景学习到的局部证据来源。对一个新样本，CCS 不直接比较其绝对表达量，而是汇总多个冻结队列模型对其 TSP 模式的判断。

> **核心判断：在 cohort model bank 固定的前提下，CCS 下的 d1 在信息完备性上不可能优于完整的 raw TSP。它的潜在优势不是创造更多信息，而是重新组织、筛选和压缩已有 TSP 信息，使其中与任务相关且可能跨队列稳定的部分更容易被有限样本、受限模型和正则化方法利用。**

若原始 TSP 矩阵为 $X$，cohort model bank 为 $F$，并且 $F$ 在评价样本及其标签进入分析前已经固定，则：

$$
d1=F(X).
$$

在上述冻结与独立前提下，$d1$ 是 $X$ 的确定性函数，数据处理不等式意味着它不会凭空创造 TSP 中不存在的信息：

$$
I(Y;d1)\leq I(Y;X).
$$

这一结论比“CCS 没有增加特征数量”更强。在无限样本、无限模型容量和理想优化条件下，直接使用完整 $X$ 的分类器可以复现任何基于 $d1$ 的决策，因此 raw TSP 在理论最优预测能力上至少不弱于 d1。若 $F$ 对任务相关信息是可逆或充分的，两者可以等价；若变换丢弃了任务相关信息，d1 的理论上限反而更低。因此，“更容易暴露稳定生物结构”不表示 d1 含有新的生物信息，而只表示原本存在于 TSP 中的结构在 d1 的几何关系或统计形式中更加显著。

有限样本学习讨论的是另一个问题。实际模型受到样本量、假设空间、正则化和优化方法限制，无法保证从高维、相关且含有队列噪声的 raw TSP 中恢复理论最优规则。若 cohort bank 能降低有效维度、减少重复计权、抑制队列特异捷径，并把跨队列一致的证据排列成更简单的决策边界，那么 d1 即使没有增加信息，也可能降低样本复杂度并获得更好的外部测试表现。这种优势应解释为**信息可利用性、归纳偏置或样本效率的改善**，而不是信息量的增加。

因此，评价 CCS 时必须区分三类主张：

- **信息优势：**在固定 $F$ 下不能成立；$d1=F(X)$ 不会比完整 $X$ 包含更多关于 $Y$ 的信息。
- **表征与学习优势：**可能成立；d1 可能让受限学习器更容易提取跨队列稳定信号，但必须通过同一数据边界下的 raw TSP 对照和 cohort-level 外部验证证明。
- **外部先验优势：**若 $F$ 由额外队列学习得到，其参数本身编码了群体经验。此时 d1 的实用增益可能同时来自表征重组和额外训练信息；若 raw TSP 对照未获得等价的数据资源，就不能把全部增益归因于 CCS 表征本身。

若评价样本或其标签参与了 cohort bank 训练，$F$ 还可能直接携带评价信息，此时不能用上述公式证明独立表示价值，必须改用 out-of-fold/out-of-cohort 预测。只有在冻结、无泄漏和公平对照成立后，才可以检验 CCS 是否把同一份 TSP 重表达为更适合有限样本学习、具有更好跨队列归纳偏置和更容易暴露稳定生物结构的空间。

## d1、d2、d3、metaCCS 与 normCCS 的角色

当前 CCS 流程中的各层应严格区分：

- **d1：**所有冻结 cohort submodels 输出块的联合表示。对于固定 TSP 规则、固定模型银行、固定基因匹配和固定软件版本，单个样本的 d1 可以独立、确定地计算。
- **d2：**按 d1 列所属的 tissue model block 分别进行第一阶段降维，再拼接各 tissue 的低维结果。这里分组的是特征块，而不是分别对各 tissue 的患者聚类。当前默认实现为每个 block 降至相同维数，拼接前的清洗、缩放和维数都应进入锁定 manifest。
- **d3：**对合并后的 d2 再做一次全局降维，形成用于可视化和聚类的统一低维坐标。
- **metaCCS：**通常在标准化 d3 上运行 DBSCAN 等聚类方法得到的泛癌状态。
- **metaCCS caller：**以 d1 为输入、以某一次锁定的 metaCCS 聚类标签为目标训练的归纳分类器。
- **normCCS：**根据 ICI response 等任务标签，对 metaCCS 状态进行排序、合并或映射得到的任务定向等级。

相关实现可参见 [`R/dr.R`](../R/dr.R)、[`R/scaller.R`](../R/scaller.R) 和 [`R/ccs.R`](../R/ccs.R)。本文对“当前实现”的判断核查于 2026-08-02；后续代码若保存了可 transform 的降维模型，应重新审计本节结论。

d2 的 tissue-first 设计确实可能减少拥有大量 cohort modules 的 tissue 在原始特征维度上的支配，但它主要平衡的是 tissue model blocks 的低维容量，不等同于直接解决各 tissue 样本量不平衡。每个 tissue 压缩到相同维数也不必然最优：它可能保护小 tissue，也可能丢失复杂 tissue 中的重要方向。

## metaCCS 与 normCCS 是否必要

需要把“计算上是否必要”和“科学叙事上是否有价值”分开。

从计算角度看，metaCCS 与 normCCS 都不是所有下游任务的必经步骤。只要 d1 是稳定、可部署的连续表示，就可以直接建立：

$$
P(Y=1\mid d1)=\sigma(\beta_0+\beta^\top d1),
$$

其中 $Y$ 可以是 ICI response，也可以是其他临床结局。由此可将 CCS 重构为：

$$
\text{TSP}
\rightarrow \underbrace{\text{frozen cohort bank}\rightarrow d1}_{\text{CCS backbone}}
\rightarrow
\begin{cases}
\text{ICI probability head},\\
\text{metaCCS state head},\\
\text{其他任务 head}.
\end{cases}
$$

在这一架构中，metaCCS/normCCS 是可选的状态读出头，而不是 CCS 定义本身不可删除的层。

从科学叙事角度看，metaCCS/normCCS 仍可能具有价值：

- 将 outcome-free 泛癌结构与 outcome-specific 临床解码分开；
- 为跨组学分析提供共同的离散索引；
- 便于报告 ORR 梯度、患者覆盖和临床分层；
- 通过强信息瓶颈降低小样本任务模型的估计方差；
- 为机制解释和临床沟通提供比二维坐标更直观的语言。

因此，发现连续任务头可用，并不自动推出 metaCCS/normCCS 应被删除。更可能的终局是“连续概率作为临床主输出，离散状态作为 companion readout”。

## 当前 d3 为什么不适合直接预测 ICI response

截至 2026-08-02，当前 [`R/dr.R`](../R/dr.R) 中的两级 UWOT/UMAP 流程只保存 d2/d3 坐标，没有保存每个 tissue-level UMAP model 和 global UMAP model。当前 d3 更接近依赖共同样本集合的 transductive atlas coordinate，而不是固定的单样本特征。

为表达其对共同样本集合的依赖，可将其概念性地近似写为：

$$
d3_i=U_{global}\left(\bigoplus_t U_t(d1_{i,t};S);S\right),
$$

其中 $S$ 是本次共同参与拟合的样本集合。改变样本集合、队列组成、随机种子或 UMAP 参数，都可能改变已有样本的坐标。固定 seed 只能提高同一数据集上的重复性，不能保证换一批样本后的可迁移性。

因此，将全部 ICI 样本共同计算 d3，再随机拆分训练集和验证集，会让验证样本的 $x$ 参与表示构造。即使未使用 response 标签，也不属于严格的单样本归纳验证。若 held-out cohort 参与 UMAP 拟合，同样存在 covariate leakage。

若未来仍希望测试 d3，至少应满足：

- 仅在开发数据中拟合所有 tissue-level 和 global UMAP；
- 保存完整的 transformable model、特征顺序、清洗和尺度参数；
- 新样本只能逐级 transform，禁止加入后重新拟合；
- 同一样本单独投影、批量投影和改变输入顺序时，d3 与最终概率应处于分析前预设的数值容差内；
- 评价跨 seed 和重采样后的预测概率稳定性，而不是逐列比较可旋转、翻转和平移的 d3 坐标。

在当前实现满足这些条件之前，不应把 `ICI response ~ d3` 作为正式预测路线。

## 如何理解 metaCCS 的稳定性

“metaCCS 稳定”至少包含两个不同概念：

- **发现稳定性：**d3+DBSCAN 得到的聚类是否对 seed、样本重采样、队列组成和参数变化稳定。
- **部署稳定性：**一旦选定某次聚类解并训练 metaCCS caller，新样本能否通过冻结的 d1→caller 路径重复得到相同标签。

当前 metaCCS caller 主要解决了部署稳定性：它把一次 transductive 聚类解蒸馏成以 d1 为输入的归纳分类器。这不代表产生该标签的原始 d3+DBSCAN 解已经具有生物学真实性或发现稳定性。

同理，normCCS 通过固定规则映射 metaCCS，保证的是映射操作可重复；它可能把上游发现的不稳定隐藏成等级跳变。metaCCS/normCCS 仍需通过 cluster Jaccard、调整兰德指数（adjusted Rand index，ARI）、状态 split/merge、grade flip 和独立生物锚点进行验证。

## 面向 ICI response 的优先模型比较

当前最小、最公平的主比较应只有三条正式路线：

| 模型 | 输入与任务头 | 回答的问题 |
|---|---|---|
| raw TSP model | raw TSP + ridge logistic | cohort bank 是否真的提供净收益 |
| d1 model | frozen d1 + ridge logistic | 连续 CCS evidence backbone 是否适合 ICI 解码 |
| normCCS model | cross-fitted/locked normCCS grade 作为唯一分子预测量，再由训练折内 logistic 映射为概率 | 离散状态路线是否以信息压缩换来稳健性和解释性 |

d3 只有在通过冻结双层 transform 的准入门槛后，才可作为附条件候选；在当前实现中不进入主比较。

三臂必须锁定同一组可用 TSP：统一基因对方向、并列值编码、缺失基因规则和 feature manifest。三臂若需要概率校准，必须采用同一训练折内校准边界。cross-fitted normCCS 遇到训练折中缺失或极少的 meta state 时，也应在查看 held-out response 前规定唯一映射或拒判规则。

主分析优先使用 ridge logistic，而不同时搜索 LASSO、随机森林、XGBoost、神经网络和多种特征筛选器。原因不是复杂模型一定更差，而是现阶段首先要判断表示层是否有价值，过大的模型搜索空间会让表示比较被研究者自由度淹没。

d1 也不是预设胜者。它可能丢失对 ICI response 有用、但对原 cohort 子模型无用的 TSP，也可能将 tissue、平台或 cohort response prevalence 作为捷径。因此，d1 必须先证明优于 raw TSP，之后才有资格讨论是否替代 normCCS。

## 验证设计与主要指标

开发阶段建议在 ICI derivation cohorts 上进行 nested leave-one-cohort-out：

- 外层每次完整留出一个 cohort；
- 内层仍按 cohort 整组选择 ridge 的正则化参数；
- 标准化、缺失处理、阈值和校准只能在外层训练数据中完成；
- 每个 cohort 在主要损失中贡献相同总权重，避免大 cohort 支配；
- normCCS mapping 也必须在每个外层训练折重新推导，不能使用全 derivation 数据优化后的映射预测该折；
- 若 ICI 样本参与过对应 cohort submodel 训练，应使用 out-of-fold/out-of-cohort d1，避免局部模型记忆。

ICI response 还需要跨队列协调：预先锁定 RECIST 或其他影像/病理评价标准、CR/PR 与 SD/PD 的合并规则、评价时间窗、治疗线次、单药或联合方案、baseline 样本定义，以及不可评价病例的处理方式。否则 LOCO 可能主要测量标签体系和治疗组成差异，而不是模型迁移能力。

主要终点建议使用 cohort-macro Brier score，因为临床目标最终是可靠概率，而不仅是病例排序。关键次要指标包括：

- cohort-macro AUROC 和 AUPRC；
- calibration-in-the-large 与 calibration slope；
- 最差 cohort 的 Brier/AUC；
- 预设阈值范围内的 decision-curve net benefit；
- leave-one-tissue-out 和条件允许时的 leave-one-regimen-class-out；
- 无输出率、缺失基因比例和队列间异质性。

不应通过多个指标加权形成事后“综合冠军分”。模型升级所需的最小临床重要差异，应由项目负责人在查看三臂结果前确定。

“d1 优于”应被写成预设的配对比较规则，例如以 cohort-macro Brier 差值为主要估计量、以 cohort 为单位进行 paired bootstrap，并要求置信区间同时越过零和预设的最小临床重要差异；同时不得出现校准、最差 cohort 或运输压力测试中的关键方向性反证。极小或单一结局类别的 cohort 仍保留 Brier 等概率指标，AUROC/AUPRC 标记为不可估，并只在可估 cohort 中汇总且单列 response prevalence。

## 路线决策树

- **d1 不优于 raw TSP：**不支持 cohort bank 对 ICI 任务具有已证净收益。优先选择 raw TSP 或外部表现更好的现有 normCCS，不继续用复杂任务头补救。
- **d1 优于 raw TSP，但不优于 normCCS：**说明 cohort-informed 表示可能有价值，但删除 metaCCS/normCCS 的证据不足，保留现有状态路线。
- **d1 同时优于 raw TSP 和 normCCS，并在锁定外部数据中复现：**可将 ICI 主线改为 `TSP → frozen cohort bank → d1 → response probability`。
- **连续概率胜出，meta/norm 仍提供可复现的跨组学或临床解释：**连续概率作为主输出，meta/norm 降为 companion state readout。
- **连续概率胜出，且 meta/norm 发现不稳定、等级频繁翻转、无独立预测/决策/机制价值：**才删除 metaCCS/normCCS 作为 ICI 推理的必经步骤。
- **三条路线在 tissue、cohort 或 regimen 间方向冲突：**拒绝宣布统一的 pan-ICI winner，报告异质性并考虑收窄适用范围。

只有 d1 胜出后，才需要进一步用 pooled model、mean vote、random/shuffled bank、去除最大队列或去除单一 tissue bank 等对照，证明增益确实来自分散的 cohort congress evidence，而不是普通 stacking、模型容量或单一大队列。

## 论文叙事与主张边界

若最终采用连续任务头，论文不应被叙述成“又一个泛癌 ICI 转录组模型”。更合适的核心定义是：

> CCS 由多个独立队列学习并冻结的局部 TSP 决策函数组成；它将单个新样本转换为可追溯的跨队列证据向量，任务头只负责读取这些预先存在的证据，而不重新训练 cohort bank。

这一叙事成立的前提是：

- d1 对开发样本采用 cross-fitting，对外部样本使用完全冻结的模型银行；
- 新样本不需要与测试队列共同归一化、嵌入或重新校准；
- 预测增益分散于多个 cohort evidence，而非由一个大队列或单一 tissue 支配；
- 相对 raw TSP 和非 CCS 基线存在稳定外部增益；
- 相关贡献可以追溯至 TSP、局部 cohort model 和最终 probability head。

需要谨慎使用若干术语：

- 单样本获得一个概率并不自动等于“personalized congress”；需要样本特异的证据构成或专家贡献支持。
- 只有一个 ICI head 只能称为 CCS 的首个 task-specific instantiation，不能证明 general-purpose multi-task backbone。
- 单臂 ICI-treated cohorts 支持的是 response prediction/association，不能在缺少对照治疗和 treatment-by-model interaction 时声称预测 ICI 相对治疗获益。
- 如果输出只有 response probability，就应使用 prediction 或 risk stratification 语言，而不是继续声称 molecular state inference。
- 若 leave-one-tissue-out 或 leave-one-regimen-out 失败，应收窄 pan-cancer 或 cross-regimen 主张。

## 当前建议

现阶段不应直接删除 metaCCS/normCCS，也不应实施 `ICI response ~ 当前 d3`。最小下一步是：

1. 冻结 TSP manifest、cohort bank、d1 schema、ICI cohort 列表、response 定义和所有验证 folds；
2. 核查 ICI 样本与 cohort bank 训练数据的患者和队列重叠；
3. 在 derivation cohorts 上完成 raw TSP ridge、frozen d1 ridge 和 cross-fitted normCCS 的 nested-LOCO 三臂比较；
4. 将三条完整流程在全部 derivation cohorts 上分别拟合并锁定，再在现有公共 external cohorts 与独立临床队列 QL1706 中同时盲测一次；三臂均不得重训练、重映射、重校准或反向调参；
5. 只有 d1 连续通过“优于 raw TSP、优于 normCCS、锁定外部复现”三道门，才修改现有论文主线；
6. 若准备提出临床部署主张，还需要一个在研究决策层面完全未触碰的新确认队列，并先完成技术重复、平台、福尔马林固定石蜡包埋（formalin-fixed paraffin-embedded，FFPE）质量、缺失基因和拒判机制验证。

因此，本轮讨论的最终结论不是简单地保留或删除某一层，而是重新界定 CCS 的必要部分：

> TSP 是具有跨平台潜力的单样本输入；冻结 cohort bank 及其 d1 是最值得检验的 CCS 核心；d2/d3 是需要归纳化和稳定性证明的表示工具；metaCCS/normCCS 是可能有价值、但应由外部证据决定其核心或辅助地位的状态读出层。

## 相关方法学启示

- UMAP 本质上是非线性邻域嵌入；若要服务新样本预测，必须显式保存和验证 out-of-sample transform：[McInnes et al., 2018](https://doi.org/10.21105/joss.00861)。
- Parametric UMAP 将可复用映射作为明确目标，说明“有坐标”与“有可部署 encoder”不是同一件事：[Sainburg et al., 2021](https://doi.org/10.1162/neco_a_01434)。
- IMPRES 表明基于基因成对关系直接预测 ICI response 具有先例，但其后续可复现性争议也说明排序特征不能替代严格外部验证：[Auslander et al., 2018](https://doi.org/10.1038/s41591-018-0157-9)；[Carter et al., 2019](https://doi.org/10.1038/s41591-019-0671-4)。
- TIDE 等工作说明跨队列转录组可以服务 ICI response prediction，但不意味着必须先建立离散无监督亚型：[Jiang et al., 2018](https://doi.org/10.1038/s41591-018-0136-1)。
