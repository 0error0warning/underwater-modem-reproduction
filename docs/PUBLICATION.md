# 公开快照说明

此仓库是研究工作目录的精选公开快照，仅包含本项目代码、文档和自产实验结果。简历、个人材料、原始Watermark文件、论文PDF、整页论文截图、运行日志、缓存、smoke输出和作废诊断不在发布范围。

2026-11-03 将首页改写为第一人称的复现记录：说明作者最初的问题、排查中改变的判断、尚未验证的硬件问题，以及 AI 工具与作者的分工；同时添加 MIT 许可证（仅覆盖本项目自己的代码与文字）。本次只改文档、许可证和公开清单，不产生新的算法或实验结果，也不修改任何实验 CSV 或 MATLAB 源码。

2026-11-02 的叙述整理将首页改为围绕问题、诊断取舍与消息收发展开；详细证据和运行入口移至 `docs/experiment_index.md`，新增 `docs/research_notes.md` 区分已有观察与尚未实施的研究设想。本次只改文档及公开清单，不产生新的算法或实验结果。

V4 在此公开仓库内继续开发：新增冻结主峰定时验证和 CRC16 双软件节点，不修改 V1–V3 的实验 CSV 或冻结接收器。新增文件的清单条目直接记录本仓库产物的 SHA-256；历史文件仍保留原始来源摘要。README、数据指南与公开核验入口随新功能更新，清单说明这些变更。

公开整理不改实验CSV数值、接收算法、冻结策略的决策参数或MATLAB源码字节。Markdown中的本机链接已改为相对链接；JSON中的本机项目路径改为仓库相对路径。冻结策略JSON仅将`development_source`改为相对路径。

`publication_manifest.json`分别记录原始实验文件的`source_sha256`与公开文件的`sha256`，并说明公开转换。历史核验JSON里的原始摘要用于描述当时的实验文件，不应假定它们就是路径清理后的文件摘要。公开快照的完整性请使用`tools/verify_public_release.py`核验。

公开版`report_continuous_phasecheck_v3.py`增加相对路径解析和公开摘要映射支持，统计函数与策略不变。`publication.py`仅接受两种情况：文件仍等于原摘要，或其当前摘要及原摘要均与明确记录的公开映射一致；它不忽略其他差异。

`results/v2/records_continuous.csv`是包含初相生成错误的历史运行，仅用于更正对照；当前连续结果是`results/v3/continuous_phase_corrected/records_continuous.csv`。NOF1主峰定时是V3诊断候选，尚未替换冻结接收器或重新拟合自适应策略。

V4 的 [预注册提交](https://github.com/0error0warning/underwater-modem-reproduction/commit/21ffd69) 先于 011–020 的下载及运行，固定 22 个代码、协议和回归依据文件；`results/v4/timing_holdout/freeze.json` 记录其摘要和当时新数据文件尚不存在的检查。正式运行前后检查相同摘要，并保留每个收到波形的配对摘要。021–030 继续保留。本证据说明本次工作流程中的先后顺序，不把同一测量序列的片段表述为独立海域。

`results/v4/packet_link/provenance.json` 记录软件协议的源码依赖与最终检查结果摘要。协议测试初次出现过空数组形状比较失败，已在 `development_checks/` 保留失败收据和说明；修正的是测试比较方式，最终验证与演示结果单独保存。历史无错载荷率仍由原离线真值模型统计，没有被新的 CRC/ARQ 演示替换。

