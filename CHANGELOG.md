# Changelog

**重要**：本文件是项目变更的**唯一正式记录**。凡是项目的更新，都要统一在本文件里记录。这是项目管理的强制性要求。

格式基于 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

## [Unreleased]

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
