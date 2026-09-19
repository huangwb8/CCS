# ablation-03 renv 支持实施计划

## 通俗解释：究竟发生了什么

`ablation-03` 目前能找到全局 R 库中的包，但项目本身没有把这些包版本记录下来。换机器、清理 R 库或更新依赖后，同一 targets 流程可能不再使用同一套软件环境。

本次为 `test/ablation-03` 增加项目级 `renv`：`renv` 负责包环境，`targets` 继续负责分析依赖图和缓存，CCS 仍通过已安装包提供正式 API。

## 专业判断：问题在哪里

- 当前目录缺少 `renv.lock`、`renv/` 和项目激活入口。
- `_targets.R` 未显式恢复环境，正式 launcher 也未检查环境状态。
- 计划要求正式复现前通过 clean-library `renv::restore()`，因此仅补文档不算完成。

## 改进方向

### 项目环境入口

在 `test/ablation-03` 增加一个幂等的 R 引导脚本，支持 `init`、`snapshot`、`restore`、`status` 和 `check`；初始化只发生在分析项目，不改变 CCS 包源码。

### 正式运行入口

增加 PowerShell launcher，统一切换到 `test/ablation-03`、显式激活 renv、检查 CCS 版本和 targets store，再执行 `_targets.R`。保留原有直接调用方式，但把 launcher 作为正式复现入口。

### 运行记录与验证

运行配置记录 renv 是否激活、lockfile 路径和 lockfile 哈希；验证包括语法检查、项目环境状态检查、targets manifest 生成和 BAC 链校验。依赖安装/restore 只由用户明确授权的命令执行。

## 实施范围与顺序

1. 增加 renv 引导脚本、正式 launcher 和项目配置。
2. 更新 `_targets.R` 与 README，要求正式运行通过 renv 入口并记录环境身份。
3. 在当前 R 4.3.1 环境安装 renv、初始化并生成 `renv.lock`，将 CCS 0.8.3 作为已安装包纳入锁定。
4. 执行静态检查、环境状态检查和 targets manifest smoke test；不运行完整大规模分析。
5. 更新 `[Unreleased]` 并验证 BAC。

## 如何确认完成

- `test/ablation-03/renv.lock`、`renv/activate.R` 和 `.Rprofile` 存在并被 Git 跟踪。
- `renv::restore()` 后可加载 `targets`、`rmarkdown`、`yaml` 和 `CCS` 0.8.3。
- launcher 在未恢复环境或 CCS 版本不符时提前失败。
- `targets::tar_manifest(script = "_targets.R")` 可在项目 renv 环境中生成。

## 风险与边界

- 不修改 `DESCRIPTION` 版本号，不安装或升级 CCS 新版本。
- 不把大体积 targets store、输入 RDS 或分析结果写入项目目录。
- 若某个二进制依赖无法在当前平台恢复，保留 lockfile 和失败日志，不绕过环境检查运行正式分析。
