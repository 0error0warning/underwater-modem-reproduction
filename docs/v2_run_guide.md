# V2：接收同步、多普勒与自适应保护时间

> 2026-11-02 更新：NOF1 4 ms 无保护的异常主要来自最早前导峰定时与符号相位复位的交互，见 [V3 诊断](../results/v3/nof1_diagnostic/report.md)。另修正了旧连续驱动的 `pi` 变量遮蔽；下面改用修正驱动，原 `run_continuous_v2` 仅为冻结历史代码。修正连续结果仍使用 V2 接收器，没有引入 V3 定时候选。

V2 在 V1 的已知时序物理层基线上加入前导捕获和多普勒估计，再研究因果自适应保护时间。旧版数据保留，方便区分每一步改变的条件。

## 接收链路

发送帧包含 HFM 上扫频/下扫频、两段已知单音导频和 BFSK 载荷。HFM 每段 64 ms，间隔 64 ms；后续导频及静默使载荷在帧开始后 528 ms 发出。2024 原论文只说明 HFM 脉冲对，没有给全套参数；64 ms 参考其引用的 [2023 同步研究](https://vibsys.put.poznan.pl/_journal/2023-34-1/articles/vps_2023102.pdf)，双音导频与具体保护安排是本实现扩展。

接收函数 `uac2.receive_packet` 只接收波形、公开帧格式、载荷长度与补偿开关：

1. 将收到的实数信号变成解析信号，下变频并降采样。
2. 用 HFM 上/下扫频的归一化匹配相关，检测未知起点的帧。
3. 只从两段已知导频估计时间尺度，重新细化 HFM 捕获位置；从前导相关峰估计多径扩展。
4. 按估计尺度对载荷作分数采样校正，再比较两个 BFSK 频率的非相干能量。

真正信道 `h`、实际速度 `V0`、发送起点和未知载荷均不会传入该函数。驱动另设有已知时序/尺度的诊断接收器，单独标记为 `oracle_timing_scale`；它不参与实用接收器或策略的决策，也不是保证最优的性能上界。

正速度的尺度约定为 `alpha = 1/(1-v/1500)`：信号时长增加，频率降低。Watermark 原始 TVIR 已去掉整体速度，因此先按官方 `sfetch` 的方式重新加回逐包尺度。连续实验在真实 TVIR 时间上推进，但没有模拟全程累计的传播延迟和完整往返声学链路。

## 自适应的因果顺序

候选保护时间为 0/4/8/16/32 ms，有效发射时间固定 8 ms。当前包的保护时间由**上一包实际收到的前导特征**决定；首包或上包未捕获时回退。模式通知假设能装入统一的 16 bit、100 bit/s 控制预算，即 0.16 s。它是明确的理想反馈模型，并未实现回传物理层或控制包解码。

训练使用每种信道 001–005，验证使用 006–010；最终 031–040 在模型冻结后才计算，中间 011–030 不使用。规则只允许一个阈值、两个动作，且不使用信道名称或真实信道参数。完整预定协议见 [adaptive_protocol_v2.md](adaptive_protocol_v2.md)。

两类实验回答不同问题：

- **固定槽位对照：**所有保护时间在相同信道时刻、相同比特与噪声下测量。报告包含 6.4 s 空槽的实验时间轴速率，同时另列按空口占用成本计算的反事实效率。
- **连续时序仿真：**下一包发送时刻实际由上一包帧长、128 ms 包尾和 0.16 s 控制预算决定。重新通过相应时刻的信道，重新捕获、估计和解调。每条记录以完整 32 s 为分母，包括末尾未用时间。这个对照能检验固定槽位上的选择是否在动态时序下仍有效。

策略是否带来提升以独立结果为准；验证选择固定策略时会保留这个决定，自适应候选仍作为预先声明的比较对象。

## 运行顺序

在项目根目录运行 MATLAB：

```matlab
addpath('matlab');
validate_receiver_v2;
validate_replay_at;
run_doppler_stress_v2;
run_link_v2('development');
```

然后在项目根目录 PowerShell 冻结模型：

```powershell
python -X utf8 tools/analyze_adaptive_v2.py --stage fit --records results/v2/records_development.csv
```

冻结后运行最终实验：

```matlab
run_link_v2('test');
run_continuous_phasecheck_v3;
```

最后生成分析：

```powershell
python -X utf8 tools/analyze_adaptive_v2.py --stage evaluate --records results/v2/records_test.csv
python -X utf8 tools/report_link_v2.py
```

`report_link_v2.py` 保留旧连续运行并显示更正提示。当前修正结果请读取 [correction_report.md](../results/v3/continuous_phase_corrected/correction_report.md)。该复跑使用已记录的原冻结策略和 `preflight.json`，不是一轮新留出评测。

只重跑本轮 NOF1 原因诊断，不必重新拟合策略：

```matlab
run_nof1_diagnostic_v3('full');          % 001–003；内置原基线逐包回归
run_nof1_diagnostic_v3('confirmation');  % 004–005；固定候选规则的小范围确认
```

然后执行 `python -X utf8 tools/report_nof1_diagnostic_v3.py` 生成诊断报告；需要约数分钟。诊断脚本不更改冻结接收器、自适应模型或原始 V2 数据。

原论文对齐诊断由 `run_paper_alignment_v2` 与 `tools/report_paper_alignment_v2.py` 单独运行，只使用开发记录 001–003，避免最终数据影响实现选择。它使用 64 bit 包，不能与主实验的 128 bit 包直接比较 PER。

`run_link_v2` 默认从完整记录检查点继续；要强制重新运行可调用第二参数 `false`。应保留冻结模型和原始结果副本，重新拟合后的模型不能再把已经查看的测试记录称为未见数据。

## 看结果时优先检查什么

- 同时看捕获丢失率、条件 BER 与整包错误率。未捕获包标为擦除，不把它虚构成零比特错误。
- 固定参考样本功率 SNR、每比特能量 Eb/N0、整帧采样 SNR 分开显示。相同数字标签不是相同噪声。
- 比较补偿与未补偿时，它们共享同一次捕获；真实时变信道中，单一尺度估计未必对每包都有益。
- 无错载荷率依据仿真真值判定整包正确。CRC、ACK 错误和重传尚未建模，因此是这个模型下的结果。
- 导频速度和恢复的 `V0` 之间的差异还包含信道残余多普勒，不应都解释为估计器错误。
- 前导相对平均功率时延谱峰值的偏移只是诊断；最早到达路径与最强路径本来可能不同。
- 最终记录来自同一次测量的较晚时间段。记录重采样范围能说明本批结果的波动，不能证明跨海域普遍有效。

## 结果入口

- [接收链路与连续时序报告](../results/v2/link_report.md)
- [自适应策略独立评测](../results/v2/adaptive/report.md)
- [原论文对齐与 FFT 诊断](../results/v2/paper_alignment/report.md)
- `results/receiver_v2/validation.json`、`long_packet_doppler.json`
- `results/v2/adaptive/frozen_policy.json`、`results/v2/receiver_code_freeze.json`

只有对应元数据标记 `COMPLETE` 或验证标记 `PASS` 的结果才代表该项已跑完。源代码、随机种子、冻结规则和逐包记录共同提供复跑依据。
