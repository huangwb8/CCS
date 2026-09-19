# ablation-03 templates 职责收敛实施计划

## 通俗解释：究竟发生了什么

- **一句话说明：** 一个原本存放网页主题资源的文件夹逐渐混入了分析运行、绘图和测试相关脚本，导致看到目录名时无法判断里面文件的真实用途。
- **具体场景：** 这类似把信纸模板、账本计算规则和验货清单都放进“信纸模板”抽屉；文件仍可能工作，但维护者容易拿错、删错或错误理解依赖关系。
- **对应到本问题：** Liquid Glass CSS/HTML 是信纸模板；workflow/checkpoint helper 是账本规则；测试脚本是验货清单。
- **改变前后：** 当前 Rmd 和分析脚本从 `templates/` 加载多类文件；调整后 `templates/` 只负责 Liquid Glass，运行时 helper 统一从 `scripts/helpers/` 加载，测试仍位于 `tests/`。

## 专业判断：问题在哪里

- **当前现象：** `test/ablation-03/templates/` 同时保存主题资源、Rmd 绘图 helper、checkpoint/恢复 helper 和未引用的 HTML 片段。
- **影响范围：** Rmd 渲染、编号分析脚本、辅助脚本及两个回归测试依赖这些路径；只移动文件而不更新全部引用会导致运行失败。
- **已知原因：** 早期从 `bensz-rmd-rules` 和 `ablation-02` 按需复制的报告模板，后来又承载了 ablation-03 的工作流基础设施，目录职责随功能增长而变宽。

## 要达到什么目标

- `templates/` 只保留 `liquid_glass_theme.css` 和 `liquid_glass_lightbox.html`。
- 仍被正式流程使用的本地 helper 移入 `scripts/helpers/`，保持项目脱离本机 skill 安装后仍可复现。
- 删除全仓库未调用的模板副本和陈旧 HTML 片段。
- 不改变科学计算、缓存语义、报告数值或 CCS 版本。

## 改进方向

### 收窄主题目录

只保留 Rmd YAML 直接引用的 Liquid Glass CSS 与 after-body HTML。这样看到 `templates/` 即可确定它只服务网页主题。

### 保留必要的项目级 helper

`nature_colors`、`theme_nature`、`render_dt` 和图表交付函数当前被 Rmd 实际调用；workflow、checkpoint 和 stage receipt helper 被正式分析入口调用。这些文件不能改为读取本机 skill 目录，否则项目在另一台机器上无法独立运行，因此迁入 `scripts/helpers/` 而不是删除。

### 删除未使用文件

`complexheatmap_template.R`、`plotly_template.R` 和 `ablation03-extension.html` 没有运行时调用点，删除它们不会改变现有流程。

## 实施范围与顺序

1. 建立 `scripts/helpers/`，迁移仍有运行时调用的 helper。
2. 同步编号脚本、Rmd、辅助脚本、测试和 helper 内部的路径与代码身份清单。
3. 删除未使用文件，并更新 README 与 `[Unreleased]` 变更记录。
4. 执行静态引用检查、helper 加载检查和相关回归测试；不启动完整昂贵分析。

## 如何确认完成

- `templates/` 只剩两个 Liquid Glass 文件。
- 除 Rmd YAML 中的 Liquid Glass 引用外，`test/ablation-03` 不再从 `templates/` 加载 R 脚本。
- 所有移动后的 R helper 可以解析并按预期加载。
- `29-test-stage-provenance.R` 与 `33-test-workflow-resume.R` 通过；其它直接受路径变化影响的轻量检查通过。
- BAC 哈希链验证无错误。

## 风险与待确认事项

- 实施前已确认本轮涉及的 ablation-03 路径不存在其它未提交冲突；后续修改继续保持在目录职责和直接引用范围内。
- 完整 Rmd 渲染依赖既有产品和 R 环境，本次先做路径与 helper 加载验证；若现有产品可用，再选择一个报告做最小渲染验证。
- `bensz-rmd-rules` v0.23.0 的 workflow checker 将 `templates/checkpoint_helpers.R` 硬编码为唯一合法位置，尚不支持项目级替代路径；本项目保持更清晰的 `scripts/helpers/` 边界，并将该兼容性问题作为 skill 设计缺陷本地留痕。
