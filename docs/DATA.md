# 获取数据与核验结果

仓库包含本项目产生的 CSV、图表及统计脚本，不包含 Watermark 原始信道、官方 MATLAB 源码或论文 PDF。来源和归属见 [THIRD_PARTY.md](../THIRD_PARTY.md)。

## 不下载数据，核验已有结果

在仓库根目录运行，只需 Python 标准库：

```sh
python tools/verify_public_release.py
```

该命令检查公开文件清单、冻结 MATLAB 源码、NOF1 误码计数、连续实验的完整32秒分母及共同初相。它核验已发布快照，不替代重新进行信道仿真。

## 获取 Watermark

请先阅读 [FFI 官方页面](https://www.ffi.no/en/research/watermark)及其手册中的说明。可以下载官方完整压缩包并解压到 `vendor/Watermark/`。本仓库也保留了按 ZIP 条目下载的工具，例如：

```sh
python tools/fetch_watermark.py --fetch Watermark/doc/ Watermark/matlab/ NOF1_001.mat NCS1_001.mat
```

各阶段所需文件：

|阶段|信道与记录|
|---|---|
|原论文对齐、NOF1主诊断|NOF1/NCS1 001–003；NOF1诊断仅用NOF1|
|NOF1追加确认|NOF1 004–005|
|开发与冻结策略|NOF1/NCS1 001–010|
|原留出评测、初相错误修正复跑|NOF1/NCS1 031–040|

PowerShell 批量获取已用到的记录：

```powershell
$channelNames = foreach ($channel in @('NOF1','NCS1')) {
    foreach ($record in ((1..10) + (31..40))) {
        '{0}_{1:000}.mat' -f $channel,$record
    }
}
python tools/fetch_watermark.py --fetch @channelNames
```

例如第001条应位于 `vendor/Watermark/input/channels/NOF1/mat/NOF1_001.mat`。官方下载源、条目清单和下载文件 SHA256 保存在 `references/watermark_*.json`。下载脚本检查 ZIP CRC，并记录 SHA256。

## 重跑与公开链接

依赖见 `requirements.txt`；MATLAB实验已在R2020b与Signal Processing Toolbox上运行。完整命令见 [运行指南](v2_run_guide.md)。报告脚本原本面向本机阅读，会生成本机文件链接；重新生成后运行以下命令转为GitHub可浏览的相对链接：

```sh
python tools/portable_reports.py
```

重新生成图表、报告或元数据会改变快照字节，因此原发布清单的哈希检查可能不再通过。请在干净克隆上核验已发布快照；将自行复跑的数据与原结果分开保存。冻结策略已查看过031–040，重新调参后不能再称这些记录为未见测试数据。

