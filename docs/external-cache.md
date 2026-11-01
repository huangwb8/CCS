# CCS 外置缓存

项目源码位于坚果云同步目录。大体积缓存和生成结果保存在 D 盘，项目内保留代码、配置、锁文件、Rmd/Markdown 源稿、必要的小型输入素材和轻量位置说明。`.gitignore` 只影响 Git，不能阻止坚果云同步，因此缓存隔离直接通过读写路径实现。

## 首轮迁移结果

首轮迁出 **24.945 GiB**，包含 96,460 个普通文件及 146 个历史目录链接，分布在 24 项迁移记录中。普通文件逐一通过大小与 SHA-256 校验，目录链接核对目标后在 D 盘内重建；指向已迁移包库的链接同步更新到新位置。

首轮结束时，项目目录按常规文件统计约 **1.151 GiB**，其中 `.git` 约 0.476 GiB、轻量历史任务及未参与归档的材料约 0.454 GiB、`test/` 约 0.210 GiB。这是首轮快照，统计不跟随目录链接。

三个正式 targets store 的元数据哈希与迁移前一致。170 个 R 文件通过语法解析，相关 Python 脚本通过 AST 检查，两个修改过的 PowerShell 启动器通过语法解析；这些检查没有执行科学分析。

## 继续迁移生成结果

随后又迁出 **170.638 MiB、401 个文件**，包括 normCCS-premium 的穷举 CSV/CSV.gz、HTML 与图表，两个 ablation 分析的生成报告和论文图表，以及旧场景的模型文件。另将一份必要的小型组图输入素材复制到外部结果目录，原件保留在项目中。402 个目标文件均通过 SHA-256 校验。

再次盘点时，项目常规文件约 **0.985 GiB**，其中 `test/` 约 **44.179 MiB**；`.git` 约 487.539 MiB，既有任务记录约 464.769 MiB。盘点排除了 NTFS reparse point，没有跟随目录联接。当前待提交的现存文件总大小约 **5 MiB**：口径为未忽略的新文件及相对 HEAD 已变更文件的完整大小之和，不含已删除文件，也不等于 Git pack 增量。

三个分析的默认输出及其读取路径已改为下表中的 `exports/`，README 和 `RESULTS.md` 提供外部结果入口。源码、报告模板、论文 Markdown 草稿及必要输入仍在项目中；结果说明中的 10 个 HTML 链接均已核对。

现有 targets 元数据和序列化对象未改写。normCCS-premium 的旧 file-target 路径在读取时可解析到新位置；下次正式运行仍由原 `_targets.R` 判断依赖是否需要更新。输出配置变化可能使依赖失效，不能保证仅重新导出报告。本次没有为更新路径而运行分析或渲染报告。

## 路径安排

| 内容 | 默认位置 |
| --- | --- |
| ablation-03 正式科学缓存和 targets store | `D:/cache/ccs/_ablation-03/`，本次保持原位 |
| ablation-03 生成报告、表格和论文图形 | `D:/cache/ccs/_ablation-03/exports/` |
| ablation-03 renv 包库 | `D:/cache/ccs/ablation-03/environment/library/` |
| ablation-03 历史试跑和隔离测试 | `D:/cache/ccs/ablation-03/tmp/` |
| ablation-03-max 正式缓存、包库和测试 | `D:/cache/ccs_ablation-03/ablation-03-max/`；包库位于 `environment/library/` |
| ablation-03-max 生成报告、products 和图形 | `D:/cache/ccs_ablation-03/ablation-03-max/exports/` |
| normCCS-premium 正式缓存和包库 | `D:/cache/ccs/normCCS-premium/`，本次保持原位 |
| normCCS-premium 穷举导出、生成报告和图表 | `D:/cache/ccs/normCCS-premium/exports/` |
| normCCS-premium 历史试跑和隔离测试 | `D:/cache/ccs/normCCS-premium/tmp/` |
| ablation-02 中间结果 | `D:/cache/ccs/ablation-02/tmp/` |
| 预训练信息取证缓存 | `D:/cache/ccs/pre-train-info/cache/` |
| 旧 CCS 场景的 project_01 模型和结果 | `D:/cache/ccs/legacy/ccs/project_01/` |
| 旧 subtypePerformance 场景结果 | `D:/cache/ccs/legacy/PADv20240810_variant_01_Select02/` |
| test 根目录中迁出的较大 RDS 和模型 JSON | `D:/cache/ccs/legacy/test-artifacts/`，文件名不变 |
| 大型历史 AI 任务产物 | `D:/cache/ccs/ai-artifacts/<原任务名>/` |

分析代码不通过目录链接访问这些文件。旧目录中的 `CACHE-MOVED.md` 只提供位置说明；包库仍保留原有 R 版本和平台子目录，不重新安装或更换包版本。

## 后续使用

ablation-03 和 ablation-03-max 的入口在激活 renv 前读取各自的 `renv-paths.R`。包库位置分别可用 `CCS_ABLATION_ENVIRONMENT_ROOT`、`CCS_MAX_ENVIRONMENT_ROOT` 调整，环境根目录下自动分出 `library/`、`cache/`、`source/`、`binary/` 与 `sandbox/`。这些位置独立于科学缓存，隔离测试切换 store 时仍使用固定包库。

ablation-03 的试跑目录可通过 `CCS_ABLATION_SCRATCH_ROOT` 调整，ablation-02 使用 `CCS_ABLATION02_SCRATCH_ROOT`。normCCS-premium 的试跑目录使用其 `CCS_NORMCCS_CACHE_ROOT` 下的 `tmp/tests/`。预训练取证的 R 和 Python 阶段共同读取 `CCS_PRETRAIN_CACHE_ROOT`；在其他操作系统运行前，应将它设为本机同步目录之外的绝对路径。

正式分析继续使用各自的 `_targets.R` 与现有启动脚本；监控连接实际运行的同一 store。为避免意外重算，本次迁移没有调用 `tar_make()`、`tar_watch()`、训练入口、报告渲染或包安装。

## 迁移与校验

2026-11-01 的迁移清单、逐文件 SHA-256 和目录汇总位于 `D:/cache/ccs/migrations/20261101-1157/`。处理顺序为复制、核对大小和 SHA-256、复核源文件未变，再移除已校验的原件。目标发生冲突时停止，不覆盖既有文件。迁移曾因符号链接权限中断；确认本次清单归属后重新校验未完成项，并通过 NTFS 目录联接恢复，原文件在整项核验完成后才移除。

第二轮生成结果的审计位于上述目录的 `generated-results/`：`files.json` 保存源路径、目标路径、字节数、SHA-256 和仅复制标志，`summary.json` 与 `verification.json` 保存迁移及复核摘要。

历史 AI 任务中的私有运行目录、会话记录及凭据类文件不参与归档，原位置保留。项目中没有新增指向外部大文件的 junction 或符号链接。

历史模型及试跑缓存可能保存了旧的绝对路径。本次保留其字节内容，未改写 RDS 或 targets 元数据；后续使用模型时应指定新模型目录，试跑仍由原 targets 入口判断依赖有效性。正式 store 保持原位。

## 恢复与备份

需要恢复时，按迁移清单将文件从 D 盘复制回原路径，核对哈希后再恢复对应路径配置。不要将归档目录直接合并进正式 store。

外置目录不会由坚果云同步。不可替代的原始输入、冻结模型和重要结果应另有备份；D 盘不可用时，停止大文件写入，不回退到项目同步目录。
