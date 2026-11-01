# CCS 项目指令

CCS（Cohort Congress System）是用于个性化泛癌基因组分类的 R 包与计算框架，重点是方法开发、验证和可复现性。

## 沟通与协作

- 默认使用简体中文，保留必要的技术名称、命令和路径；README 的双语约定见下文。
- 先明确需求、范围、产物和约束。数据分析通常按“数据获取 → 探索分析 → 模型训练 → 验证评估”推进。
- 只修改任务相关文件，先读上下文再编辑，遵循现有代码、日志和通知风格，不做无关重构或添加未要求的功能。
- 说明简洁自然；引用文件用可点击路径，必要时带起始行号。大文件只给路径和变更说明，交付时简述结果与后续步骤。
- 优先依据项目文件和本地上下文；需要联网时先用本地搜索工具，不足时再用其他联网方式，并说明原因、保留关键链接。

## 工程原则

- 保持简单，只实现当前需要的功能；复用重复逻辑，分离职责，避免过度设计。
- 接口行为符合直觉，无效输入尽早返回，减少嵌套。
- 冲突时优先保证正确性；简洁性优于灵活性，清晰性优于性能，扩展性优于代码紧凑。
- 文档标题不加序号。新增规则放入对应章节，避免重复；更新工作流或术语时，同步检查示例、验证要求和引用。

## 目录与代码规范

| 路径 | 用途 |
| --- | --- |
| `R/` | 包源码与导出 API；`ccs*.R` 定义 S4 类，`hyperTuningGS.R`、`importance.R` 负责建模，`normalize.R`、`phenotype.R`、`simulatePA.R` 等负责预处理 |
| `man/` | roxygen2 生成的 Rd 文档 |
| `test/` | 可复现场景脚本与小型测试夹具；模型和大型数据使用外置缓存 |
| `docs/` | 解释性文档、审计材料与教程；正式任务计划统一放在 `docs/plans/` |
| `docs/contribution.bac` | 通过 `bac` 命令维护的协作贡献账本 |
| `.bensz-api/` | AI 与技能的轻量中间产物，已由 `.gitignore` 排除 |

提示词、许可证和 README 保留在根目录；使用外部数据的脚本旁应说明数据路径。

R 代码使用两空格缩进、`<-` 赋值；函数与变量使用 `snake_case`，S4 类和导出方法沿用 `CamelCase`。保持结构和类型安全，补全 roxygen 的 `@description`、`@param`、`@return`、`@examples`。优先复用 `DESCRIPTION` 已声明的依赖，新增导入前先更新声明。

## targets 分析流程

- 已设计、声明或规划 `targets` 依赖图的分析，必须通过原 `_targets.R`、其加载函数、正式存储目录和 `tar_make()` 执行，通过 `tar_watch()` 监测。不得用自定义运行器、逐脚本调度或临时检查点产生正式结果。
- 正式入口、缓存边界、恢复逻辑或监控面板与 `_targets.R` 不一致时，先停止正式计算并修复，不能另跑“等价”流程。
- `tar_watch()` 必须指向本次实际使用的正式存储目录，不得用旧目录、空目录或自定义检查点冒充进度。删除、替换或降级现有流程须经负责人明确授权，并同步更新计划、README、CHANGELOG 和验证要求。
- 目标级并行必须通过 `crew` 控制器接入 `tar_option_set()`。各工作进程日志、CPU/内存指标及 targets 主进程运行日志统一写入正式缓存根目录的 `observability/`，可由 `autometric::log_read()`、`autometric::log_plot()` 读取或可视化；不得混入科学结果缓存或用自定义并行运行器替代 `crew`。

## 外置缓存与大文件

项目位于坚果云同步目录。编写代码和首次运行前，自动将大型产物规划到 `D:/cache/ccs/<分析名>/`，无需再次询问。按用途划分 `raw/`、`models/`、`targets/`、`tmp/`、`environment/`、`observability/` 和 `exports/`。沿用已启用的正式缓存根目录，如 `D:/cache/ccs/_ablation-03/`、`D:/cache/ccs_ablation-03/ablation-03-max/`，不得为统一命名重建存储目录。

- 单文件预计达到 50 MiB、单目录达到 200 MiB，或会持续增长时，提前外置。训练模型、表达/概率矩阵、RDS 中间对象、targets 存储目录、renv 包库、构建副本和批量试跑产物均适用，不必等到超过阈值。
- `test/`、根目录 `tmp/` 和 `.bensz-api/` 仅保留源码、配置、锁文件、小型夹具、轻量日志和交付说明。正式 HTML/PDF、批量 CSV/CSV.gz、模型结果、展示图及配套资源统一放入对应缓存的 `exports/`；项目内保留 Rmd/Markdown 源稿、位置说明及必要的小型素材。大型 AI 任务副本或归档放入 `D:/cache/ccs/ai-artifacts/<任务名>/`。
- 用统一路径配置和环境变量连接缓存，生产、测试、报告与监控各自使用对应的同一套路径；隔离测试仍走原 `_targets.R`。renv 库路径在激活前设置，`TMPDIR`、`TMP`、`TEMP` 在 R 进程启动前设置。不得以目录联接、符号链接或 Git 忽略规则替代缓存隔离。
- 迁移前确认没有相关写入并检查目标冲突；复制后核对文件数、大小和哈希，再移除原件。不覆盖不同内容，不删除唯一副本。迁移正式 targets 存储目录时，同步全部入口、监控和文档；只读迁移不得重跑分析。
- D 盘不可用或空间不足时停止大文件写入并说明问题，不得回退到同步目录。外置缓存不随坚果云同步，唯一原始输入、不可替代模型和重要结果须另有独立备份。

## 开发与验证

Windows 下统一使用 `C:\R\R-4.3.1` 对应的 R/Rscript。构建副本、包检查和试跑产物遵循外置缓存约定。

| 操作 | 命令 |
| --- | --- |
| 刷新文档（提交前） | `Rscript -e "devtools::document()"` |
| 加载开发源码 | `Rscript -e "devtools::load_all()"` |
| 构建与检查 | `R CMD build .`，再运行 `R CMD check CCS_<version>.tar.gz`；版本从 `DESCRIPTION` 读取 |
| normCCS 单元测试 | `Rscript --vanilla test/normCCS-5-more/run-unit-tests.R` |
| 场景测试示例 | `Rscript test/03.test.Classes.R`、`Rscript test/test.classifier_performance.R` |

- 根目录不保留 `tests/`；normCCS 轻量单元测试放在 `test/normCCS-5-more/tests/`，其他测试归入对应场景目录。`test/` 不随 R 包打包，包检查前必须单独运行相关测试入口。
- 场景脚本按 `01.test.*`、`02.test.*` 等编号表达执行顺序，能通过 `Rscript` 非交互运行。大型 `.rds` 和临时试跑产物写入对应外置缓存，临时产物放在 `tmp/`。
- 使用 `stopifnot` 或轻量 `testthat` 断言；涉及概率输出或亚型图变化时，补充小型样本快照。验证范围应与变更相符。

## 文档与报告

- README 以英文 `README.md` 为主，简体中文 `README_CN.md` 为辅，首屏提供语言切换链接。两份同步维护，章节层级、技术事实、示例、命令、版本和链接对应一致，分别自然表达。未经负责人明确要求，不得反转主辅语言、改名或按技能默认命名覆盖此约定。
- 代码变更使 `docs/` 中非 `plans/` 文档过时时，必须同步更新；修改 `AGENTS.md` 后，检查 `CLAUDE.md` 的核心规则一致性。
- 所有新写或修订的报告，交付前必须使用 `humanization` 技能润色（`C:/Users/Administrator/.codex/skills/humanization/SKILL.md`），减少套话、生硬措辞和重复表达。此规则也适用于非 NSFC 报告，按报告用途组织内容，不套用标书结构。
- 润色只调整表达，保留事实、数字、公式、引用、代码、术语和科学结论，并复核 Markdown/LaTeX/R Markdown 结构。正式分析报告仍通过既有 `targets` 入口更新。

## 版本、变更与提交

- `DESCRIPTION` 的 `Version` 是版本号唯一来源。版本、发布时间、Git 标签及 GitHub Release 由负责人决定；除下述既有授权外，只有当前请求明确指定目标版本并要求修改时，才能同步修改 `DESCRIPTION`、README 版本徽章、CHANGELOG 发布标题及标签说明，不得自行推断或递增版本。
- **ablation 重构既有授权**：新的 `R/ablation.R` 完成改造，通过声明的包检查、科学等价性验证及 `test/ablation-03` 串行验收后，才将当前版本第三位修订号恰好加 1，不改主/次版本，不提前升版。随后用 `C:\R\R-4.3.1` 构建、检查并安装到该 R 环境的可写包库；源码和构建临时文件不得混入 R 安装目录。
- 行为、结构、工作流、工程原则、指令或关键配置的变更，必须记入 `CHANGELOG.md` 的 `[Unreleased]`；功能、修复或文档更新本身不构成升版授权。CHANGELOG 遵循 Keep a Changelog，Git 标签格式为 `v{Version}`，创建或推送均须明确授权。
- 提交说明简短、使用祈使句，必要时关联问题编号。拉取请求说明目的、涉及的数据和参数、包检查及相关场景测试结果；绘图或统计变更附截图或指标，说明未纳入 Git 的大型产物。强制推送变基结果前与审阅者协调。

## 贡献记录

使用 [bensz-auto-contribution](https://github.com/huangwb8/bensz-auto-contribution) 的 `bac` 工具记录人类、AI 与工具的贡献边界，要求 Python 3.10+，账本为 `docs/contribution.bac`。

- 开始实质工作时，用 `bac input record` 记录脱敏后的用户指令摘要；文件修改、工具执行和验证用 `bac record` 记录事件与证据，引用路径必须位于项目内。
- 交付前运行 `python -m bac --root . --bac-file docs/contribution.bac verify --json` 验证哈希链；审阅时间线用 `python -m bac --root . --bac-file docs/contribution.bac inspect --limit 20`。
- 不直接编辑账本，初始化和修复均通过 `bac` 命令。仅负责人明确要求时，才可在 `init-project` 中使用 `--disable-bac`，并在系统文档与 CHANGELOG 说明。修改贡献规则时，同步验证账本并记录变更。
- BAC 仅为过程记录和辅助审计，不替代署名、责任或合规判断；不得记录密钥、令牌、完整私有提示词、个人健康信息、专有基因组数据或无关隐私。

## 数据安全

基因组 `.rds` 按敏感数据处理，个人健康信息和专有数据不得进入仓库；样本夹具清理敏感元数据，拉取请求说明安全存储要求。远程下载须由显式参数启用，凭据或令牌不得写入跟踪文件或贡献账本。
