# 公开快照说明

此仓库是研究工作目录的精选公开快照，仅包含本项目代码、文档和自产实验结果。简历、个人材料、原始Watermark文件、论文PDF、整页论文截图、运行日志、缓存、smoke输出和作废诊断不在发布范围。

公开整理不改实验CSV数值、接收算法、冻结策略的决策参数或MATLAB源码字节。Markdown中的本机链接已改为相对链接；JSON中的本机项目路径改为仓库相对路径。冻结策略JSON仅将`development_source`改为相对路径。

`publication_manifest.json`分别记录原始实验文件的`source_sha256`与公开文件的`sha256`，并说明公开转换。历史核验JSON里的原始摘要用于描述当时的实验文件，不应假定它们就是路径清理后的文件摘要。公开快照的完整性请使用`tools/verify_public_release.py`核验。

公开版`report_continuous_phasecheck_v3.py`增加相对路径解析和公开摘要映射支持，统计函数与策略不变。`publication.py`仅接受两种情况：文件仍等于原摘要，或其当前摘要及原摘要均与明确记录的公开映射一致；它不忽略其他差异。

`results/v2/records_continuous.csv`是包含初相生成错误的历史运行，仅用于更正对照；当前连续结果是`results/v3/continuous_phase_corrected/records_continuous.csv`。NOF1主峰定时是V3诊断候选，尚未替换冻结接收器或重新拟合自适应策略。

