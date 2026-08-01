# CCS - Cohort Congress System - 项目指令

本项目是用于个性化泛癌基因组分类的 R 包与计算框架。

## 项目目标

A Computational Framework for Personalized Pan-cancer Genomic Classification.

本工作区的核心目标是开发、验证并维护可复现的个性化泛癌基因组分类方法。

## 核心工作流

当用户提出数据科学或 R 包开发相关需求时，按以下流程执行：

### 任务理解

- 理解用户的真实需求和意图
- 确认任务范围和预期输出
- 识别可能的依赖和约束

### 执行流程

数据获取 → 探索分析 → 模型训练 → 验证评估

### 输出规范

- 代码变更应遵循项目现有风格
- 文档更新应保持一致性
- 测试覆盖应符合项目标准

## 项目目录约定

- `R/`：R 包核心代码与导出 API
- `man/`：由 roxygen2 生成的 Rd 文档
- `test/`：可复现的场景测试、模型与测试夹具
- `docs/`：解释性文档、审计材料与教程
- `docs/plans/`：AI 为特定任务制定的正式计划
- `docs/contribution.bac`：BAC 协作贡献账本，必须通过 `bac` 命令维护
- `.bensz-api/`：AI/skill 中间产物，已由 `.gitignore` 排除

## 工程原则

本项目遵循以下工程原则：

| 原则 | 核心思想 | 在本项目中的体现 |
|------|----------|------------------|
| **KISS** | Keep It Simple, Stupid | 追求极致简洁，避免过度设计；文档标题不使用序号前缀（用 `##` 而非 `## 1)`） |
| **YAGNI** | You Aren't Gonna Need It | 只实现当前需要的功能 |
| **DRY** | Don't Repeat Yourself | 相似逻辑应抽象复用 |
| **SOLID** | 面向对象设计五大原则 | 单一职责、开闭原则等 |
| **关注点分离** | Separation of Concerns | 不同层次逻辑应分离 |
| **奥卡姆剃刀** | 如无必要，勿增实体 | 优先选择最简单的解决方案；Markdown 本身有层级结构，序号是冗余的形式化标记 |
| **最小惊讶原则** | Principle of Least Astonishment | API 行为应符合用户直觉 |
| **早期返回原则** | Early Return | 尽早返回，减少嵌套 |

**原则冲突时的决策优先级**：
1. **正确性 > 一切**
2. **简洁性 > 灵活性**
3. **清晰性 > 性能**
4. **扩展性 > 紧凑性**

## 默认语言

除非用户明确要求其他语言，始终使用 简体中文 与用户对话与撰写文档/说明。

## 联网与搜索

默认优先使用项目内文件与本地上下文；确需联网获取信息时，优先使用本地搜索工具。仅当本地工具不足以满足需求时再使用其它联网手段，并说明原因与保留关键链接。

## 贡献记录

本项目默认且强制基于 [bensz-auto-contribution](https://github.com/huangwb8/bensz-auto-contribution) 使用 `bac` 工具，客观记录人类、AI 与工具的协作贡献边界。默认账本为 `docs/contribution.bac`，运行环境要求 Python 3.10+。

- 开始实质性工作时，用 `bac input record` 记录经过脱敏的用户指令摘要
- 发生文件修改、工具执行或验证时，用 `bac record` 记录对应事件与证据；文件路径必须位于项目内
- 交付前运行 `python -m bac --root . --bac-file docs/contribution.bac verify --json`，确认哈希链完整
- 需要审阅贡献时间线时，运行 `python -m bac --root . --bac-file docs/contribution.bac inspect --limit 20`
- 不直接编辑 `.bac` 文件；初始化或修复必须通过 `bac` 命令完成
- 仅当项目负责人明确要求时，才可在 `init-project` 中使用 `--disable-bac`；关闭状态必须在系统文档与 `CHANGELOG.md` 中说明
- BAC 是过程记录与辅助审计材料，不替代最终署名、责任或合规判断
- 不记录密钥、访问令牌、完整私有提示词、PHI、专有基因组数据或无关个人隐私

## Codex CLI 特定说明

### 文件与输出

- 引用文件时使用可点击路径，并尽量带起始行号：`src/main.py:42`
- 不输出刚写入的大文件内容，只引用路径并说明变更
- 简单确认避免复杂格式；结束时给出简短后续步骤

### 编辑原则

- 只修改与当前任务直接相关的文件
- 保持现有代码风格、结构和类型安全
- 读取足够上下文后再批量处理相关修改
- 无效输入早返回，遵循项目既有日志/通知模式
- 不主动添加用户未要求的功能

## 变更记录与版本

- `DESCRIPTION` 中的 `Version` 是项目版本号的唯一来源，但该事实不授予 AI 修改版本号的权限
- 软件何时发布、采用什么版本号、创建什么 tag 以及是否执行发布，只能由项目负责人本人决定；任何 AI 都不得推断、选择、递增或代替决定
- 除非项目负责人在当前请求中明确给出目标版本并要求修改，否则 AI 严禁修改 `DESCRIPTION` 的 `Version`、README 版本徽章、`CHANGELOG.md` 的发布版本标题、Git tag 或 GitHub Release
- 功能新增、问题修复、文档更新、BAC 引入等变更即使按照 SemVer 可能需要升版，AI 也只能记录到 `[Unreleased]`，不得自行调整版本号；如有必要可说明影响，等待项目负责人决定
- 项目负责人明确要求更新版本时，AI 才可按其指定值同步 `DESCRIPTION`、README、`CHANGELOG.md` 与 tag 说明，不得擅自改成其它版本
- 影响项目行为、结构、工作流、工程原则、指令文件或关键配置的变更，必须更新 `CHANGELOG.md` 的 `[Unreleased]`
- 修改 `AGENTS.md` 后，应同步检查 `CLAUDE.md` 的核心内容是否一致
- `CHANGELOG.md` 遵循 Keep a Changelog；Git tag 使用 `v{Version}`，但只有项目负责人明确授权后才能创建或推送

## Project Structure & Module Organization

Core package code lives in `R/`, with `ccs*.R` owning S4 classes, `hyperTuningGS.R` and `importance.R` for modeling, and utilities such as `normalize.R`, `phenotype.R`, and `simulatePA.R` for preprocessing. Generated Rd files remain in `man/`, while the exploratory scripts, saved models, and fixtures under `test/` mirror real workflows (see `test/ccs/project_01/`). Keep prompts, licensing, and README material at the repo root, and add any new auxiliary data next to the script that consumes it.

## Build, Test, and Development Commands

Refresh documentation before committing with `Rscript -e "devtools::document()"`. Package locally with `R CMD build .`, then run `R CMD check CCS_<version>.tar.gz` against the tarball version read from `DESCRIPTION`. During interactive work `Rscript -e "devtools::load_all()"` mirrors `library(CCS)` without reinstalling. Scenario tests run via plain scripts; e.g., `Rscript test/03.test.Classes.R` validates the class definitions and `Rscript test/test.classifier_performance.R` exercises the classifier benchmarks.

## Coding Style & Naming Conventions

Match the existing tidy R style: two-space indent, `<-` for assignment, and `snake_case` for functions and variables while keeping S4 classes and exported methods in `CamelCase`. Keep roxygen headers complete (`@description`, `@param`, `@return`, `@examples`) so `man/` stays in sync after `devtools::document()`. Reuse verbs from the packages declared in `DESCRIPTION`, and update that file before relying on new imports.

## Testing Guidelines

Each script in `test/` targets a reproducible scenario; follow the numbering scheme (`01.test.*`, `02.test.*`, etc.) to signal execution order. Scripts must be runnable through `Rscript` without interactive prompts and should cleanly write their `.rds` artifacts into dedicated subdirectories (e.g., `test/ccs/project_<id>/`). Favor `stopifnot` or lightweight `testthat` expectations to make failures obvious and add sample data snapshots when a change impacts probability outputs or subtype plots.

## Commit & Pull Request Guidelines

Git history favors concise imperative summaries with optional release tags (`v0.7.0`, `Routine update`). Use the same style, referencing issue IDs when relevant. Pull requests should include: purpose, datasets or parameters touched, confirmation that `R CMD check` and the relevant `test/*.R` scripts succeeded, and screenshots or metrics for any plotting/statistics change. Mention large artifacts deliberately excluded from git and coordinate with reviewers before force-pushing rebases.

## Security & Data Handling

Handle genomic `.rds` files as sensitive. Keep PHI and proprietary data out of the repo, scrub metadata in sample fixtures, and document any secure storage requirements in the PR description. Gate remote downloads behind explicit arguments and never embed credentials or tokens in tracked files or the BAC ledger.

## 有机更新原则

当需要更新本文档时：

- 先理解变更意图，再把规则放到最合适的章节，避免重复堆叠
- 更新工作流、输出规范或术语时，同步检查相关示例、验证清单和引用位置
- 计划文档统一保存在 `./docs/plans/`
- 代码变化导致 `docs/` 中非 `plans/` 文档过时时，必须同步更新
- 修改贡献记录规则时，同步检查 `docs/contribution.bac` 的可验证性与 `CHANGELOG.md` 记录
- 文档标题不使用序号前缀，保持 Markdown 层级清晰
