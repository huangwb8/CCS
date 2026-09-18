# Changelog

**重要**：本文件是项目变更的**唯一正式记录**。凡是项目的更新，都要统一在本文件里记录。这是项目管理的强制性要求。

格式基于 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

## [Unreleased]

- 修复包文档检查阻断：移除两处不可执行的示例占位文本，将依赖外部 `resCCS` 对象的示例标记为不自动运行，将 XGBoost 参数范围中的 Unicode 无穷符号改为 LaTeX 可移植的 `Inf`，并补充 `plotImportance()` 的 `nTop` 参数说明。

- 为 `ablation()` 增加可序列化的 `context`、`plan`、`run`、`result` 阶段调度，
  并通过 `cache.root` 将中间节点缓存与正式 `output.dir` 分离；保留默认
  `step = "all"` 的普通调用语义；默认 `.ccs-cache/` 已加入忽略规则，未修改版本号。

- 新增 `ablation-03` 的 targets 串行迁移骨架：`R/ablation.R` 提供不依赖
  workflow/SUCCESS/stage receipt 的阶段化 `ablation()` API；
  `test/ablation-03/_targets.R` 将正式 store 固定到
  `D:\cache\ccs\_ablation-03\targets`，并要求通过已安装的 `CCS::`
  API 运行。本轮只完成代码与接口改造，未运行完整 `ablation-03`，待审阅后
  再进行小规模 smoke 与串行验收；未修改版本号。

- 移除未被正式 targets 流程消费的 `ablation_*` 公共导出和独立 targets API 文档；
  job 规划保留为内部实现，运行身份元数据移回 `ablation-03` wiring，公共入口收敛为 `ablation()`。

- 更新项目指令：明确 ablation 重构通过包、科学等价性和 `ablation-03` 串行验收后，才递增 CCS patch 版本，并使用 `C:\R\R-4.3.1` 对应环境完成构建、检查和安装；本门禁不适用于验证前的开发迭代。
- 明确 ablation-03 的正式运行必须消费已安装的 CCS 包 API；不得直接 source `R/ablation.R`，开发期 `load_all()` 仅用于包测试，不作为正式分析入口。
- 明确新版 targets ablation-03 的正式缓存根目录为 `D:\cache\ccs\_ablation-03`；targets store、分支对象和运行状态写入该外部目录，项目目录仅保留源码、配置和正式报告。

- 新增 `test/ablation-03/scripts/run-fresh-analysis.ps1`，为 Windows R 4.3.1 提供隔离的全新分析缓存根目录、UTF-8 locale、受限 CPU/内存预算、可选 benchmark 与七阶段顺序入口；不覆盖既有 ablation-03 缓存。

- 整理 `test/ablation-03` 目录边界：将重复的 `input/` 说明并入 `raw/`，将生物学锚点配置归入 `raw/config/`，将辅助入口从 `tools/` 统一为 `scripts/`，并把既有正式图表归位到 `reports/figures/`；保留根目录 `analysis-plan.yaml` 作为标准工作流清单，不迁移或删除历史 `tmp/` 缓存。

- 优化 ablation 表示分析执行：readout 内层折叠的标准化/模块平衡矩阵按 fold 只计算一次并复用于全部 lambda；learning curve 改为确定性 job table，按 fraction/repeat/representation 保存可恢复 checkpoint，并通过显式 `CCS_ABLATION_WORKERS` 提供有界 PSOCK 并行。大型只读矩阵在 Windows worker 启动时只传输一次，可用 `CCS_ABLATION_MEMORY_GB` 按估计峰值自动下调 worker；默认仍为单 worker，避免隐式改变既有运行语义。

- 按 `bensz-rmd-rules` 重构 `ablation-03` 为七个 `AA.BB.CC` 分析单元；大体积缓存迁移到 `D:/cache/ccs/_ablation-03`，项目内仅保留可审计的轻量产品清单，正式图表迁入 `reports/figures/`。表示分析新增 retrieval、readout、learning curve 与 decoder 的内容键多版本 checkpoint、运行状态和断点恢复契约；节点键只绑定实际递归代码依赖，state/key/file/value hash 任一异常均安全重算，并在 manifest 审计命中状态及失效原因。四份 Rmd/HTML 同步改用新入口。证据门禁现同时核验 reference/query d1 来源，reference 为 in-sample 时降为描述性；biology 和结构推断对 cohort/dyad 不足的结果只保留点估计，不再输出误导性区间或 P 值。

- 重构 ablation-03 为 `01a/01b/01c` 输入准备与 `02/03/04` 缓存分析；样本契约和结构表达提取前移，原表示计算体抽为可消费已准备输入的内部函数，保留计算方法、参数与 seed。将总报告拆为数据概览、表示性能、生物学锚点和结构可复现性四份 Rmd，移除旧入口与旧 HTML，更新引用及来源核验。本次仅静态修改与检查，实验已按负责人要求停止，新计算和报告渲染留待代码审核后执行。

- ablation-03 biology 阶段收敛为纯 R 工作流，移除无实质计算作用的 Python 入口，避免产生多套分析实现和结果来源歧义。

- 修复 ablation-03 审核发现的表示尺度与评估总体问题：readout/learning curve 统一使用原始 reference/query，继续复用 CCS 已计算的 d1；连续 anchor 覆盖全部候选 query，固定完整 reference 标准化并按 query 配对；按名称解析 signature，拒绝分位边界并列的伪状态。
- 加强 ablation-03 可复现性：恢复后的 tissue 标签用于 bank scaling；缓存增加完整输入与验证参数指纹，取消证据不足的旧几何缓存提升；三个分析阶段增加来源核验，Python biology 统一转交 R 实现；同步回归测试、报告解释与运行文档。

- 更新 ablation-03 实验报告的精确与 Monte Carlo sign-flip 方法说明并重新渲染 HTML；完整枚举保持尾部比例，Monte Carlo 使用加一校正。结构复现设计表新增每次重复的 seed 与 n_repeats，并对重复 cohort key 显式报错。

- 为 ablation representation 的精确原生几何增加严格键控的持久缓存，旧版精确产物需重算后建立完整内容键，不能仅凭 manifest 提升，避免高维 kNN 在重复运行中成为数小时限速步骤。

- 为 `ablation()` representation 流程增加持久化 Direct-GSClassifier 特征缓存：首次运行生成 `direct-feature-cache.rds`，后续按表达矩阵、样本、模型元数据和 feature manifest 哈希安全复用；缓存失配或损坏时才重建，且 d1 继续只复用 CCS 对象已有的 `Data$Probability$d1`。

- 对齐 `docs/ablation.md` 与 `R/ablation.R` 当前实现：明确预计算 d1、矩阵 metadata、默认 `cancer_type` anchor 和 `cohort` 分派契约，移除不存在的 d1 自动编码 helper，并同步源码定位与输出文件清单。
- 将 `ablation-03` 的 biology-anchor 输入审计并入 `01-ablation03-test-data.Rmd`：复用完整表达矩阵派生缓存，报告 anchor 基因覆盖、缓存样本/cohort 范围与 cohort 级明细，新增逐 cohort 的 anchor 覆盖率分布图，并保留 32 特征主数据画像与 biology 输入边界的区分。
- 收紧 `ablation()` 的 d1 输入边界：移除下游函数内部的 query d1 自动重算；缺失 d1 的 query 仅发出警告、按已有 d1 取交集，并将排除样本写入 `excluded-query-d1.csv`，由调用方负责准备一致的 CCS object。
- 重构 `R/ablation.R` 的消融实验编排：公共入口直接区分 representation 与 layered 实验流程，将 layered 运行拆分为上下文准备、实验调度和结果汇总阶段，并补充阶段性审阅注释；保持现有实验顺序、Gate 1 依赖和输出文件契约不变。
- 恢复 `ablation()` 内置的合成数据 smoke fixture：默认关闭，设置 `CCS_ABLATION_RUN_SMOKE=true` 时可直接运行并校验当前 layered 编排、四类 cohort 指标和审计表。

- 优化 `ablation-03` external query 资格边界：共享准备层保留全部候选 cohort，新增 endpoint-specific `candidate`/`estimable`/`not_estimable` 审计表；癌种相关 retrieval、technical excess、readout 与 learning curve 继续使用预先声明的跨 cohort 支持门槛，而 geometry、连续 anchor 与结构复现使用完整候选集合；manifest、audit 和回归测试同步记录候选数与可估计数。
- 修复 `ablation-03` 环境初始化在无活动图形设备的 `Rscript` 会话中调用 `par()`、从而偶发生成空 `Rplots.pdf` 的问题，并新增图形设备副作用回归测试。
- 将 `ablation-03` 独立队列结构验证扩展为双向外部投影：新增 150 个 reference modules → external samples 与 43 个 external modules → reference samples 的分向结果、完整样本/module/tissue 审计，以及共同 tissue 内 20 次 22 对 22 module-bank 匹配敏感性分析；同步新增 Figure 14、回归测试、机器可读产物与 HTML 解读，并明确 externality 由 bank 与 target 是否重叠定义。
- 改写 `ablation-03` 独立队列结构分析的 R Markdown 解读：增加一句话结论、地图类比、癌种协调说明、逐图读法和统计边界的通俗解释，并同步更新报告总结与 HTML。
- 扩展 `ablation-03` 的独立队列结构可复现性分析：以 4 个表达锚点的队列内高低状态构建 8 个共享生物实体，比较 Direct-GSClassifier 与 Cohort-d1 的跨队列质心距离几何；新增 136 个队列对的机器可读结果、cohort-node bootstrap、不依赖分类性能的配对图与双热图，以及同癌种独立队列稳健性分析。
- 完成 `ablation-03` 出版级统计与图形复核：cohort readout 改用单类别场景可解释的 accuracy 差，并将重采样单位明确为 query cohort；澄清 pooled balanced accuracy 是类别等权而非样本加权；scaling 斜率改为按 module 数每翻倍估计；学习曲线按唯一 training-cohort design 推断并将 100% 比例的伪重复 CI/P 值标记为不可估计；精确 sign-flip 检验移除仅适用于 Monte Carlo 的加一校正；同步优化图注、字号、置信区间标记及异质重建指标的分面展示。
- 修复 `ablation-03` HTML 在 Windows R 4.3.1 下的 Unicode 乱码：渲染前显式启用可用的 UTF-8 locale，避免无效 `C.UTF-8` 使 R 输出退化为转义字节。
- 修正 `ablation-03` 出版前一致性问题：移除会混入陈旧 biological-anchor 图表与硬编码数值的 HTML 扩展；同步 README 产物清单；将 readout 的单类别 query-cohort 推断改名为 cohort accuracy，并明确 macro-AUROC 不可估计。
- 将 biological-anchor 分数改为 reference-cohort 全局 gene-wise z-score（同一尺度应用于 query/reference），并记录共同有效邻居与缺失邻居审计，避免跨 cohort 标准化和静默丢失造成不可比效应量。
- 改善 ablation-03 图形可读性：修正 Figure 8 分面尺度、Figure 9 横轴标注、Figure 2 caption 换行、Figure 7 图例布局、Figure 5 技术因子标签、Figure 6 字号及 Figure 1 不可估计注释。

- 修正 `ablation-03` 检索指标的 MRR@k 计算：每个候选 `k` 仅计入前 `k` 名内的首次同标签邻居，并补充回归测试；重新渲染出版前 HTML 与矢量图。

- 严格化 `ablation-03` biological-anchor 比较：保留原始 utility，同时以外部 query cohort 为统计单位新增 d1−Direct 效应量、cohort bootstrap 95% CI、精确配对符号置换 P 值和 Benjamini–Hochberg 多重比较校正，并将推断结果写入 `anchor_inference.csv` 与 HTML 报告。
- 严格化 `ablation-03` technical-neighbor excess 比较：对 assay、platform、source 三个预先指定技术因子复用 cohort-level 配对 bootstrap/sign-flip 推断，并在 HTML 中报告 Holm 校正后的 P 值与区间。

- 补充 Windows 开发环境约定：R 环境路径为 `C:\R\R-4.3.1`，用于 `Rscript`、包构建、检查和场景测试。

- 优化 `ablation-03` 生物学锚点流程：新增一次性表达矩阵/冻结 signature 子集缓存，记录来源、签名与样本键哈希及覆盖审计；R/Python readout 默认只消费缓存，并在缺失或失效时显式报错，不再静默回退读取完整矩阵。

- 重构 `ablation()` 的 d1 输入边界：移除少见场景专用的 `d1_source` 与
  `external-d1-cache.rds`，统一接收单一 CCS 对象；对象已有的 d1 行直接复用，
  缺失的 query 行才按冻结 model bank 按需编码。`ablation-03` 在数据入口将完整
  `resCCS` 的 d1 行按筛选版冻结列合同合并后再调用通用 API。

### Changed

- 修复 `test/ablation-03` 在其本地 `.Rproj` 工作目录下运行时的相对路径解析；统一 R/Rmd/Python 脚本的项目根目录、缓存和外部输入定位，并补充可调试的运行约定。

- 完善 `test/ablation-02` 消融报告的首次指标解读：在各结果模块补充计算口径、通俗含义、动态实际数值与解释边界，并澄清截断 MRR、Jaccard 邻居重合和重建误差分数等易误读概念。
- Figure 2 现以与指标颜色一致的虚线显示样本级总体 top-k 一致率与 MRR 差值。

### Added（新增）

- 新增 `docs/plans/2026-08-08-ablation-02-cohort-bank-scaling优化计划.md`：规划统一 Direct-GSClassifier 主合同、将 cancer-type 指标降为 lineage 诊断，并以 tissue breadth × within-tissue cohort depth 的二维设计重构 d1 cohort-bank scaling。

### Changed（变更）

- 将 `test/ablation-02` 的 `10`–`16` 契约测试脚本集中迁入 `tests/`，并同步更新该目录 README 的运行命令，以明确主分析入口与质量保障脚本的边界。
- 重构 `test/ablation-02` 的分析入口：用 `00`、`01`、`02` 文件名前缀明确环境、数据审计与消融实验顺序，将原单体报告拆分为同名的两组 R/Rmd/HTML，并将契约测试按 `10`–`16` 编号；新增目录 README 说明执行顺序与产物位置。
- 优化 `test/pre-train-info` 的组织标签交付：对源合并 RDS 中 45 个历史 `Undefined` cohort，按各自原始组件 RDS 的顶层标签恢复 xlsx 中的 `tissue` 与 `cancer_type`，同时保留缓存层的源标签与逐行对齐证据。
- 重新校准 `test/ablation-02` 的 d1 消融报告解读：将癌种标签明确为诊断锚点而非 CCS 优化目标，不再以癌种可恢复性下降单独宣判 d1 过度混合，并补充 feature–cohort 融合、新状态效用与 cohort-axis scaling 的证据边界。
- 将预训练元数据流水线入口由 Windows 专用 PowerShell 脚本替换为跨平台 Python 编排器；外部 RDS、cBioPortal、UCSC Xena 与 GEO 路径改为命令行或环境变量配置，并支持 Windows Excel 或三平台 LibreOffice 公式重算。
- 按项目负责人确认的口径，将 assay 与托管来源已经充分确认的 cBioPortal cohort 记为 `metadata_status = confirmed`；数据提供者未公开具体测序仪或芯片型号时保留 `platform_id = unknow`，不再仅因此列入人工审核，并重建预训练元数据工作簿与审核报告。
- 接受有可靠来源支持的平台家族作为预训练元数据的 `platform_id`：Chin 2006 与 Hess 2006 标记为 `Affymetrix U133 family`，Vijver 2002 标记为 `Agilent/Rosetta custom microarray family`；不再因原始研究未公开精确芯片型号而将这些 Xena cohort 保留为 `unknow` 或列入平台人工审核。
- 按项目负责人确认的口径调整预训练元数据状态：已匹配官方 UCSC Xena dataset 的 RNA-seq cohort 即使缺少统一测序仪型号，也保留 `platform_id = unknow` 并将 `metadata_status` 记为 `confirmed`，不再列入平台人工审核；报告同时保留各数据集在官方 metadata 中实际记录的 TPM、FPKM-UQ、FPKM 或 CPM 单位。
- 更新预训练样本信息取证流程：将本地 UCSC Xena 表达矩阵精确匹配到官方 hub/dataset metadata 与数据页，确认 Caldas 2007 队列使用 `Agilent Human 1A (V2) microarray`；证据只支持平台家族时使用家族级 `platform_id`，完全缺少可靠平台证据时保留 `unknow`，并重建相关 CSV、审核报告与工作簿。
- 移除 `test/ablation/ablation.R` 的重复实现；`R/ablation.R` 现为 `ablation()` 的唯一事实来源，消融测试继续直接加载包源码。
- 重构 `ablation()` 参数 schema：将共享参数归入 `general`，并按 `cohort`、`scaling`、`tissue_first`、`metaccs` 组织任务特异参数；保存的 `config.rds` 与 `result$config` 使用规范嵌套结构，同时兼容现有扁平调用，并对未知字段及新旧写法冲突提前报错。
- 更新 `ablation()` 文档与合成场景测试，使参数示例、Gate 1 路径及各实验配置读取与新 schema 一致。
- `ablation()` 的 Direct 表征改为复用 GSClassifier 的原生特征构造合同，按冻结模型的 `bst$feature_names` 重建单基因分箱、普通基因对和 gene-set 对比；Experiment 4 的对应组名更新为 `Direct-GSClassifier`。

### Fixed（修复）

- 修正合成测试中 `table` 属性导致的错误严格比较，并让 metaCCS 特征夹具同步冻结模型实际使用的 `bst$feature_names`。
- 修复消融实验仅手工计算普通 TSP、遗漏 GSClassifier 单基因分箱与 gene-set 对比，以及跨队列表达合并错误使用基因交集的问题。

## [0.8.2] - 2026-08-01

### Added（新增）

- 新增 `docs/plans/2026-08-01-metaccs-end-to-end-ablation.md`：规划在 `CCS::ablation()` 中增加 Direct-TSP 与 Cohort-d1 的端到端、等机会 metaCCS 流程比较。
- 引入 BAC 贡献记录系统：新增 `docs/contribution.bac`，用于记录人类、AI 与工具的协作贡献证据。

### Changed（变更）

- 调整 `test/ablation/ablation.R` 与设计文档中的 Gate 1 策略：将 `params` 扩写为覆盖全部顶层及嵌套字段、默认值、适用实验和约束的完整 API 说明，补充 Gate 判定及 cutoff 用法，并将 `gate1$enforce` 默认值改为 `FALSE`，使 scaling 默认保留探索性结果、确认性分析可按需启用预注册门槛。
- 按 `init-project` 2.3.3 规范优化 `AGENTS.md`、`CLAUDE.md` 与 README：确立 `AGENTS.md` 为通用指令唯一来源，并补充 BAC 使用、安全、验证与维护规则。
- 修改 `AGENTS.md` 的版本管理规则：软件版本号、发布时间、tag 与发布操作只能由项目负责人明确决定，AI 不得自行推断、选择或变更。
- 优化 `test/ablation/ablation.R` 中 `ablation()` 的内嵌测试：统一合成数据、冻结模型和 d1 的 tissue/cohort 语义，并用最小可复现参数覆盖 cohort、scaling、tissue-first 与 metaCCS 四类实验
- 将 `DESCRIPTION` 版本更新为项目负责人指定的 `0.8.2`：记录本次向后兼容的消融工作流扩展与内嵌测试修正
- 调整了 README.md 中教程开放状态的英文措辞：将 “the full tutorial” 细化为“完整且完全可用的形式”，使含义更准确
- 更新了 README.md 中的教程说明：补充论文正式发表前教程不会完全公开的状态说明，避免读者对开放范围产生误解
- 修改了 AGENTS.md 中的"版本号管理规范"章节：将版本号来源从 config.yaml 改为 R 包的 DESCRIPTION 文件，新增 Git Tag 命名规范（`v{Version}` 格式）
- 完善了 README.md：从简单占位符升级为完整的国际化项目文档，包含安装指南、快速开始、API 文档、开发指南、算法概述、引用格式等章节
- 重写了 README.md：首页改为简洁项目入口，突出 `ccs.principle` 教学仓库与 `GSClassifier` 生态链接，并将推荐引用调整为 `Under review` 占位说明
- 调整了 README.md 的作者信息展示：改为与 `GSClassifier` 一致的单行作者样式，并保留邮箱、博客和 ORCID 链接

### Fixed（修复）

- 撤销 AI 未经授权将版本从 `0.8.1` 提升至 `0.9.0` 的变更，并将 `DESCRIPTION` 与 README 版本徽章恢复为项目负责人随后决定的 `0.8.2`。
- 修正 `DESCRIPTION` 与 `AGENTS.md` 中项目英文描述的 `Framwork` 拼写错误，并移除自动合并产生的错误 README 徽章片段与游离版本示例。

## [0.7.4] - 2026-08-01

### Added（新增）

- 为 `test/ablation/ablation.R` 中的 `ablation()` 增加纯代码合成数据示例：在函数开头通过默认关闭的 `if (FALSE)` 区块构造 CCS 对象、表达矩阵、冻结模型与元数据，可直接执行 cohort smoke test

### Changed（变更）

- 将 `DESCRIPTION` 版本从 `0.7.3` 更新为 `0.7.4`：记录本次向后兼容的测试示例增强

## [1.0.0] - 2026-03-01

### Added（新增）

- 初始化 AI 项目指令文件
- 生成 `CLAUDE.md`（Claude Code 项目指令）
- 生成 `AGENTS.md`（OpenAI Codex CLI 项目指令）
- 配置项目工程原则和工作流

### Changed（变更）

### Fixed（修复）

---

## 记录规范（强制性要求）

### 必须记录的变更类型

每次修改以下内容时，**必须**在本文件追加记录：

1. **项目指令文件变更**
   - CLAUDE.md 的任何修改
   - AGENTS.md 的任何修改

2. **项目结构变更**
   - 新增/删除/重命名目录
   - 新增/删除/重命名关键文件（如核心源码文件、配置文件）

3. **工作流变更**
   - 核心工作流程的调整
   - 开发流程的修改

4. **工程原则变更**
   - 新增工程原则
   - 修改或删除现有工程原则

5. **重要配置变更**
   - 影响项目行为的配置文件修改
   - 依赖关系的重大变更

### 记录格式

```markdown
## [版本号] - YYYY-MM-DD

### Added（新增）
- 新增了 XXX 功能/章节：用途是 YYY

### Changed（变更）
- 修改了 XXX 章节：原因是 YYY，具体变更内容是 ZZZ
- 修改了项目目录结构：将 ABC 目录移至 DEF 位置

### Fixed（修复）
- 修复了 XXX 问题：表现是 YYY，修复方式是 ZZZ

### Deprecated（即将弃用）
- XXX 功能将在下一版本移除：原因是 YYY

### Removed（已移除）
- 移除了 XXX 功能：原因是 YYY

### Security（安全）
- 修复了 XXX 安全漏洞：影响是 YYY
```

### 记录时机

- **修改前**：先在 `[Unreleased]` 部分草拟变更内容
- **修改后**：完善变更描述，添加具体细节和影响范围
- **发布时**：将 `[Unreleased]` 内容移至具体版本号下

### 版本号规则

遵循语义化版本（Semantic Versioning）：

- **主版本号（Major）**：重大架构变更、不兼容的 API 修改
- **次版本号（Minor）**：新增功能或章节，向后兼容
- **修订号（Patch）**：修复问题或微调，向后兼容

### 变更类型说明

| 类型 | 说明 | 示例 |
|------|------|------|
| Added | 新增的功能或章节 | "新增了 `## 变更记录规范` 章节" |
| Changed | 对现有功能或内容的变更 | "修改了 `## 工程原则` 章节，增加了早期返回原则" |
| Deprecated | 即将移除的功能（警告） | "旧的目录结构将在下个版本重构" |
| Removed | 已移除的功能 | "移除了已废弃的 `## 代码审查` 章节" |
| Fixed | 修复的问题 | "修复了模板中目录树生成的 bug" |
| Security | 安全相关的修复 | "修复了依赖包的安全漏洞" |

### 质量标准

每条记录应该：
- **清晰具体**：说明改了什么、为什么改
- **可追溯**：包含足够的上下文信息
- **格式统一**：遵循上述模板
- **及时更新**：修改后立即记录，不要拖延
