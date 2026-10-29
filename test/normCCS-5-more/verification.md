# 验证记录

本文保留 **0.8.5 的历史验收记录**。当前 0.8.6 的目标函数扩展及最新安装、测试与包检查结果见 [本轮验证说明](objective-verification.md)。

所有输入由 `fixtures.R` 生成，未读取真实患者或基因组文件。测试使用项目指定的 R 4.3.1；包目标版本仅为负责人指定的 0.8.5。

## 已完成的结果核对

冻结 0.8.4 的 `subtypeNorm()` 后，对 K=2、3、4，四种 RR/NRR 排序以及 Raw／OptimizeRank 两种 ROC 模式，分别核对当前默认调用与显式穷举。两者的完整 `Raw`／`Summary` 均与冻结参考一致。

小问题的全部切分用于检查快速评分与原 ROC 函数、加权 DP 与完整穷举的最优损失。用例覆盖稀疏与纯结局状态、缺少部分状态的队列、局部反序、pROC 自动方向、OptimizeRank 的浮点并列，以及反向 factor levels。DP 的 RNG、options、工作目录与并行后端保护另有测试。

公开 `subtypePerformance()` 的五级 DP 与 `n_norm_subtype=NULL` 均已实际调用，检查了参数透传、Repeat 记录和返回结构。相关脚本为 `scripts/tests/public-api.R`。

## 算法迭代与性能

单切点版本在 60 个质量问题中有 53 个达到穷举最优；加入相邻双切点后，60 个问题均在 1e−10 容差内匹配。第一轮记录与冻结实验函数保存在 `reference/`，可运行 `compare-refinement.R` 复现。

独立计时数据采用 M=12、C=4、种子 307，每种配置各运行三次。K=4／5／6、两种 ROC 模式的配对加速比中位数为 21.87–44.00 倍。五级 Raw 为 7.33 → 0.18 秒，五级 OptimizeRank 为 12.83 → 0.42 秒；六级 Raw 为 8.41 → 0.19 秒，六级 OptimizeRank 为 17.08 → 0.44 秒。六级 OptimizeRank 仍有 0.003425495 的 AUC 差距，已在正式报告保留。

这次调试用例中的匹配不构成对任意输入的全局最优保证。`exact_global_optimum=FALSE` 描述最终中位数搜索，`dp_subproblems_exact=TRUE` 只描述生成初值的加权 Raw AUC 子问题。

## 执行状态

107／216 状态、20 队列、K=5/6/8/10/12 的两种 ROC 模式共 20 次高维搜索全部完成，耗时 5.28–31.61 秒，主进程峰值工作集最高 1,238 MiB。所有最终候选均逐队列通过原 ROC 函数核验；高维任务没有运行原版穷举，组合数不作为实测加速比。

| 检查 | 结果与证据 |
|---|---|
| 旧版完整结果回归 | PASS：K=2/3/4，四种排序、两种 ROC 模式，默认与显式穷举 |
| 全切分评分与 DP 子问题 | PASS：包括种子 18 的浮点排序与反向 factor levels |
| 状态与公开入口 | PASS：RNG/options/工作目录/并行后端、DP 与 NULL 分支 |
| 初版结果复现 | PASS：60 个存档单切点目标值全部复现 |
| 合成质量基准 | PASS：最终 60/60 在 1e−10 内达到穷举目标 |
| 直接计时与高维运行 | PASS：完整记录位于 `products/`；保留近似误差 |
| 同入口隔离计算与渲染 | PASS：`tmp/tests/normCCS-20261029-125318/verification.json` |
| 正式报告 | PASS：`report.html` 和 `reports/figures/speedup.pdf` 已生成并检查图形 |
| 真实患者全量验证 | NOT_RUN：本任务只采用代码生成的数据 |

隔离验收比较了配置、环境脚本、冻结参考及代码哈希，也核对了正式产品与报告的文件清单、哈希和修改时间，以及 raw 目录元数据；均未变化。正式报告随后由单独的渲染命令生成。

## 包检查

`devtools::document()` 已刷新 `man/subtypePerformance.Rd`。包源码由 R、man、tests 及明确列出的包元数据构建，暂存文件与项目文件逐项比对哈希，未带入患者数据、历史工作区、环境文件或本地运行历史。源码包为 `CCS_0.8.5.tar.gz`，保存在本任务工作区的 `bensz-r-developer/output/`。

最终 `R CMD check --no-manual --no-build-vignettes` 结果为 **0 errors、7 warnings、1 note**，包测试通过。本次新增的搜索函数没有留下命名空间检查提示。七项 warning 分别涉及既有许可证写法、R 源码非 ASCII 字符、依赖声明与旧调用、S3 plot 签名、Rd 交叉引用、已有导出对象缺文档及旧参数文档；note 汇总原有全局变量和函数绑定提示。它们位于原有代码和文档，本次未扩展到全包清理，因此不能称为无警告的 CRAN 检查。

已将构建包安装到 R 4.3.1 默认可写用户库，确认 `packageVersion("CCS")` 为 **0.8.5**。另起 R 进程从已安装包运行公开五级 DP 示例并通过断言，未使用源码加载代替安装验证。没有创建 tag 或执行发布。

构建、检查、安装及最终组件日志均在本任务工作区 `bensz-r-developer/log/`；包检查完整记录位于 `bensz-r-developer/output/CCS.Rcheck/00check.log`。

BAC 校验的哈希链完整、无错误；账本仍有 10 个历史事件的 actor/source 归属警告，本次没有改写历史记录。

初次直接对整个工作目录构建时，R 的复制阶段遇到历史临时目录的长路径，随后改用上述源码暂存构建。自动审批拒绝了那次失败构建临时副本的清理，返回原因为策略阻止；该副本没有用于交付包或安装。
