# bensz-router 接入说明

2026-10-07 已按[官方安装指南](https://router.benszresearch.com/install.zh-CN.md)配置 Codex 用户级 MCP 与 hooks。沿用已有专用 Python 环境中的客户端 1.1.4 和现有授权，未重新登录，也未改动已安装技能源码。

## 已配置内容

MCP 服务以固定名称 `bensz_router` 注册在用户级 `~/.codex/config.toml`，通过绝对客户端路径连接 `https://router.benszresearch.com`。用户级 `~/.codex/hooks.json` 包含 SessionStart、PostToolUse、SessionEnd、SubagentStart 和 SubagentStop；PostToolUse 只匹配 router 的受管 MCP 工具，相关上下文预算保持 `additionalContextLimit: 0`。

客户端位于 `%LOCALAPPDATA%\bensz-router\venv\Scripts`，该目录已加入用户 PATH。为确保 Windows MCP 和 hooks 的中文消息使用 UTF-8，用户环境已设置 `PYTHONUTF8=1`；这也会影响此用户后续启动的其他 Python 进程的默认文本编码。已有进程需要重启才能继承环境变更。

## 验证结果

`hooks doctor` 返回 `registered: true`、无配置冲突、`setup_complete: true`；MCP 状态为 `ready`，进程启动和全部七个工具定义验证通过。另用字节流独立校验了 MCP 初始化与工具列表响应，严格 UTF-8 解码和 JSON 解析均通过。远程目录和本人默认方案读取成功。

七个工具为 `catalog`、`profiles`、`preflight`、`run`、`current`、`advance` 和 `release`。指南中的 `activate` 属于未发布源码构建，已发布的 1.1.4 不提供该工具。本次没有发起预检或启动工作流。

当前会话尚未加载新注册的工具和 hooks，诊断仍为 `host_verified: false`、`host_loaded: false`、`reload_required: true`。这些状态不能用安装或进程测试替代。

## 首次使用

完整退出并重新启动 Codex 或 IDE，使用户环境生效，然后新建会话。在 `/hooks` 审核并信任本次安装的 hooks，在 `/mcp` 核对 `bensz_router` 和上述七个工具。信任审核是宿主要求，不能由安装脚本代为标记完成；参见 [OpenAI hooks 文档](https://developers.openai.com/codex/hooks/)与 [MCP 文档](https://developers.openai.com/codex/mcp/)。

如需在尚未继承新环境的 PowerShell 中复查，可运行：

```powershell
$env:PYTHONUTF8 = '1'
$routerCli = Join-Path $env:LOCALAPPDATA 'bensz-router\venv\Scripts\bensz-router.exe'
& $routerCli --server 'https://router.benszresearch.com' hooks doctor --harness codex --scope user --transport mcp --json
```

## Windows 配置兼容处理

客户端 1.1.4 在无空格的 Windows 绝对路径上生成未加引号的 hook 命令，而所有权校验使用 POSIX `shlex.split`，导致反斜杠被移除，`hooks doctor` 报 `invalid managed hook definitions`。本次仅给五个 hook 命令的可执行文件路径加双引号，并同步其所有权记录，未改动事件、匹配规则或客户端源码。复查通过。

该兼容处理会使 `upgrade_available: true`，因为命令格式与安装器原始输出不同；它不表示存在新版软件。对 1.1.4 再次执行 `hooks install` 可能恢复有问题的命令，升级或重装后需检查路径引用并重新运行诊断。

此外，默认 Windows Python 编码下，MCP 工具列表可能不是 UTF-8，而同样使用本地编码的进程自检仍会通过。本次通过用户级 UTF-8 模式解决，并增加独立的协议字节校验。两项缺陷均已通过 bensz-collect-bugs 留作本地记录，未公开上传。

脱敏诊断证据保存在 `.bensz-api/task-20261007-1802-router-check/bensz-router/log/doctor-sanitized.json`。
