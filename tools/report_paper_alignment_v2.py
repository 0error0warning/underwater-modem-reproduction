"""Report all disclosed paper-alignment conditions without selecting a fit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEYS = ['channel', 'on_ms', 'guard_ms', 'snr_label_db', 'noise_mode']


def aggregate(data, keys):
    grouped = data.groupby(keys, sort=True, dropna=False)
    result = grouped.agg(packets=('payload_bits', 'size'), detected_packets=('detected', 'sum'),
                         payload_bits=('payload_bits', 'sum'), bit_errors=('bit_errors', 'sum'),
                         bit_erasures=('bit_erasures', 'sum'), packet_errors=('packet_error', 'sum')).reset_index()
    result['decoded_bits'] = result.payload_bits-result.bit_erasures
    result['conditional_ber'] = result.bit_errors/result.decoded_bits.replace(0, np.nan)
    result['capture_loss'] = 1-result.detected_packets/result.packets
    result['packet_loss'] = result.packet_errors/result.packets
    result['success_payload_bits'] = (result.packets-result.packet_errors)*64
    return result


def md_table(data):
    def fmt(value):
        if isinstance(value, (float, np.floating)):
            return 'NA' if pd.isna(value) else f'{value:.5g}'
        return str(value)
    lines = ['|'+'|'.join(data.columns)+'|', '|'+'|'.join(['---']*len(data.columns))+'|']
    lines.extend('|'+ '|'.join(fmt(x) for x in row)+'|' for row in data.itertuples(index=False, name=None))
    return '\n'.join(lines)


def save_figure(fig, out, stem):
    fig.savefig(out/f'{stem}.png', dpi=170, bbox_inches='tight')
    fig.savefig(out/f'{stem}.svg', bbox_inches='tight')
    plt.close(fig)


def plot_alignment(summary, reference, out):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5), sharex=True, sharey=True)
    modes = ['sample_snr_packet', 'ebn0_received']
    for ci, channel in enumerate(['NOF1', 'NCS1']):
        for ni, mode in enumerate(modes):
            ax = axes[ci, ni]
            for equal_guard, color in [(False, '#c75f33'), (True, '#217b8b')]:
                data = summary[(summary.channel == channel) & (summary.noise_mode == mode) & (summary.snr_label_db == 20)]
                data = data[(data.guard_ms == data.on_ms) if equal_guard else (data.guard_ms == 0)]
                label = 'Tg=Ton' if equal_guard else 'Tg=0'
                ax.plot(data.on_ms, data.conditional_ber, 'o-', color=color, label=f'V2 {label}')
                paper = reference[(reference.channel == channel) & (reference.snr_label_db == 20)]
                paper = paper[(paper.guard_ms == paper.on_ms) if equal_guard else (paper.guard_ms == 0)]
                values = pd.to_numeric(paper.ber.str.replace('<', '', regex=False))
                ax.plot(paper.on_ms, values, 'x--', color=color, alpha=.65, label=f'Paper {label}')
            ax.set_title(f'{channel} | {mode} | label 20 dB', fontsize=10)
            ax.set_yscale('symlog', linthresh=0.001)
            ax.grid(alpha=.25); ax.set_xticks([4, 8, 16, 32])
            ax.set_xlabel('Active symbol duration (ms)'); ax.set_ylabel('Conditional payload BER')
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Descriptive comparison; noise definitions are not inferred from closeness\nPaper <0.0001 is plotted at its upper boundary; V2 zero remains zero', fontsize=11)
    fig.tight_layout()
    save_figure(fig, out, '01_paper_comparison')

    compact = aggregate(pd.read_csv(out/'per_packet.csv'), ['channel', 'noise_mode', 'snr_label_db'])
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), sharey=True)
    colors = {'NOF1': '#217b8b', 'NCS1': '#c75f33'}
    for ni, mode in enumerate(modes):
        ax = axes[ni]
        for channel in ['NOF1', 'NCS1']:
            data = compact[(compact.channel == channel) & (compact.noise_mode == mode)]
            ax.plot(data.snr_label_db, data.capture_loss, 'o-', color=colors[channel], label=f'{channel} capture loss')
            ax.plot(data.snr_label_db, data.packet_loss, 's--', color=colors[channel], label=f'{channel} packet loss')
        ax.set_title(mode); ax.set_xticks([0, 10, 20]); ax.set_ylim(-.03, 1.03)
        ax.set_xlabel('Declared dB label'); ax.set_ylabel('Fraction of transmitted packets'); ax.grid(alpha=.25)
    axes[0].legend(fontsize=8)
    fig.tight_layout(); save_figure(fig, out, '02_capture_and_packet_loss')


def plot_fft(data, out):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, on_ms in zip(axes, [4, 8]):
        subset = data[(data.on_ms == on_ms) & (data.ebn0_db == 10)]
        labels = [f'{int(x.f0_hz)}/{int(x.f1_hz)}' for x in subset.itertuples()]
        x = np.arange(len(subset))
        ax.bar(x-.18, subset.exact_ber, .36, color='#217b8b', label='Exact-tone energy')
        ax.bar(x+.18, subset.fft_ber, .36, color='#c75f33', label='Nearest N-point FFT bin')
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15); ax.set_title(f'Ton={on_ms} ms; Eb/N0=10 dB')
        ax.set_xlabel('Transmit tone pair (Hz)'); ax.set_ylabel('AWGN BER'); ax.grid(axis='y', alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Identical received samples; rectangular windows; no zero padding', fontsize=11)
    fig.tight_layout(); save_figure(fig, out, '03_fft_alignment_awgn')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=ROOT/'results/v2/paper_alignment')
    args = parser.parse_args(); out = args.results.resolve()
    metadata = json.loads((out/'metadata.json').read_text(encoding='utf-8-sig'))
    if metadata['status'] != 'COMPLETE':
        raise ValueError('Simulation has not completed')
    data = pd.read_csv(out/'per_packet.csv')
    assert set(data.record) <= {1, 2, 3}, 'Alignment must not read held-out records'
    assert set(data.noise_mode) == {'sample_snr_packet', 'ebn0_received'}
    assert (data.payload_bits == 64).all()
    assert ((data.bit_erasures == 0) | (data.bit_erasures == 64)).all()
    assert (data.bit_errors <= data.payload_bits-data.bit_erasures).all()
    assert not data.duplicated(KEYS+['record', 'slot']).any()
    is_full = metadata['stage'] == 'full'
    if is_full:
        assert len(data) == 1440, 'Expected 2 channels x 3 records x 8 configurations x 5 slots x 3 labels x 2 definitions'
    # Two definitions differ only by their declared power/energy denominator.
    paired = data.pivot(index=KEYS[:-1]+['record', 'slot'], columns='noise_mode', values='noise_variance')
    lengths = data[data.noise_mode == 'sample_snr_packet'].set_index(KEYS[:-1]+['record', 'slot']).power_reference_samples
    ratio = paired.ebn0_received/paired.sample_snr_packet
    assert np.allclose(ratio, lengths.reindex(ratio.index)/(2*64))
    summary = aggregate(data, KEYS)
    per_record = aggregate(data, KEYS+['record'])
    summary.to_csv(out/'summary.csv', index=False)
    per_record.to_csv(out/'per_record_summary.csv', index=False)
    reference = pd.read_csv(ROOT/'references/paper_bfsk_reference.csv', dtype={'ber': str})
    reference['paper_is_upper_bound'] = reference.ber.str.startswith('<')
    reference['paper_value'] = pd.to_numeric(reference.ber.str.replace('<', '', regex=False))
    compared = summary.merge(reference.rename(columns={'ber': 'paper_ber_reported'}), on=KEYS[:-1], validate='many_to_one')
    if is_full:
        assert len(summary) == 96 and len(compared) == 96
        assert (summary.packets == 15).all()
    compared['ber_difference_for_exact_paper_values'] = np.where(compared.paper_is_upper_bound, np.nan, compared.conditional_ber-compared.paper_value)
    compared.to_csv(out/'paper_comparison_all_conditions.csv', index=False)
    alignment_metrics = []
    worst_conditions = []
    for mode, values in compared.groupby('noise_mode', sort=True):
        eligible = values[~values.paper_is_upper_bound & values.conditional_ber.notna()].copy()
        eligible['absolute_ber_difference'] = eligible.ber_difference_for_exact_paper_values.abs()
        alignment_metrics.append({
            'noise_mode': mode,
            'all_conditions': len(values),
            'exact_reference_conditions': int((~values.paper_is_upper_bound).sum()),
            'compared_exact_conditions': len(eligible),
            'no_decoded_ber_conditions': int(values.conditional_ber.isna().sum()),
            'paper_upper_bound_conditions': int(values.paper_is_upper_bound.sum()),
            'mean_absolute_ber_difference': eligible.absolute_ber_difference.mean(),
            'max_absolute_ber_difference': eligible.absolute_ber_difference.max(),
        })
        if not eligible.empty:
            worst_conditions.append(eligible.loc[eligible.absolute_ber_difference.idxmax()])
    alignment_metrics = pd.DataFrame(alignment_metrics)
    worst_conditions = pd.DataFrame(worst_conditions)
    alignment_metrics.to_csv(out/'alignment_error_summary.csv', index=False)
    worst_conditions.to_csv(out/'largest_ber_deviations.csv', index=False)
    compact = aggregate(data, ['noise_mode', 'snr_label_db'])
    compact.to_csv(out/'noise_definition_summary.csv', index=False)
    fft = pd.read_csv(out/'fft_awgn.csv')
    plot_fft(fft, out)
    if is_full:
        plot_alignment(summary, reference, out)
    shown = compact[['noise_mode', 'snr_label_db', 'packets', 'bit_errors', 'decoded_bits', 'conditional_ber', 'capture_loss', 'packet_loss']]
    fft_shown = fft[fft.ebn0_db == 10][['on_ms', 'f0_hz', 'f1_hz', 'bits', 'exact_errors', 'fft_errors', 'exact_ber', 'fft_ber']]
    zeros = summary[(summary.bit_errors == 0) & (summary.decoded_bits > 0)]
    unavailable = int((summary.decoded_bits == 0).sum())
    all_rows = compared[['channel', 'noise_mode', 'on_ms', 'guard_ms', 'snr_label_db', 'paper_ber_reported', 'bit_errors', 'decoded_bits', 'conditional_ber', 'capture_loss', 'packet_loss']]
    channels_label = '/'.join(sorted(data.channel.unique()))
    record_label = '、'.join(f'{int(x):03d}' for x in sorted(data.record.unique()))
    alignment_figures = (f'![20 dB 对齐图]({(out/"01_paper_comparison.png").as_posix()})\n\n'
                         f'![捕获与整包损失]({(out/"02_capture_and_packet_loss.png").as_posix()})') if is_full else '本目录是 smoke 自检输出，不包含完整对齐图。'
    recovery_note = ''
    if (out/'runtime_recovery.json').exists():
        recovery = json.loads((out/'runtime_recovery.json').read_text(encoding='utf-8'))
        recovery_note = (f'运行中 MATLAB R2020b 进程在 {recovery["completed_rows_before_restart"]} 条完整观测落盘后出现 '
                         '`0xc0000374 Heap corruption`。按完整配置检查点续跑剩余记录；没有修改实验条件、种子、噪声公式或接收器。'
                         '恢复范围和原检查点 SHA256 保存在 `runtime_recovery.json`。这次进程恢复不涉及最终留出数据。')
    standout = compared[(compared.channel == 'NOF1') & (compared.noise_mode == 'sample_snr_packet') &
                        (compared.on_ms == 4) & (compared.guard_ms == 0) & (compared.snr_label_db == 20)]
    deviation_note = '本 smoke 子集不含完整的最大偏离诊断。'
    if not standout.empty:
        row = standout.iloc[0]
        deviation_note = (f'**仍有明显未对齐：NOF1、4 ms、无保护、`sample_snr_packet`、20 dB，'
                          f'本次条件 BER 为 {int(row.bit_errors)}/{int(row.decoded_bits)}={row.conditional_ber:.5f}，'
                          f'原文为 {row.paper_value:.4f}，绝对差 {abs(row.ber_difference_for_exact_paper_values):.5f}。** '
                          '不能仅凭 8 ms 等长保护等条件接近原表，就认为整体数值已经对齐。')
    worst_shown = worst_conditions[['noise_mode', 'channel', 'on_ms', 'guard_ms', 'snr_label_db',
                                   'bit_errors', 'decoded_bits', 'conditional_ber', 'paper_ber_reported',
                                   'absolute_ber_difference', 'capture_loss']]
    ebn0_worst = worst_conditions[worst_conditions.noise_mode == 'ebn0_received'].iloc[0]
    worst_capture_note = (f'`ebn0_received` 的最大偏差条件只解出 {int(ebn0_worst.detected_packets)} 包，'
                          f'失锁 {int(ebn0_worst.packets-ebn0_worst.detected_packets)}/{int(ebn0_worst.packets)} 包，'
                          f'条件 BER 使用 {int(ebn0_worst.decoded_bits)} bit；它不能代表全部发送载荷的误码率。')
    report = f'''# V2 原论文参数对齐诊断

本次为开发数据上的描述性诊断（阶段 `{metadata['stage']}`），使用 {channels_label} 的 {record_label} 记录，**未读取最终留出 031–040**。完整诊断与原文相同的是信道名称、250 Hz 双音间隔、4/8/16/32 ms 有效发射时间解释、无保护/等长保护，以及 0/10/20 dB 标签；包长、频率位置、前导与导频、相位、同步、噪声含义仍包含自行指定的条件，不能称为精确数值复现。

{deviation_note}

## 尚未对齐的幅度

下表按两种噪声定义分别汇总。MAE 是各可比条件的条件 BER 绝对差等权平均，只纳入原文给出精确 BER 且本次至少有一个载荷获得判决的条件；不把 `<0.0001` 当成精确值，也不替无判决条件补零。无 BER 条件数与原文界限条件数可能重叠。

{md_table(alignment_metrics)}

各噪声定义下绝对差最大的可比条件为：

{md_table(worst_shown)}

{worst_capture_note} MAE 也不计入失锁代价，且两口径可比较条件集合不同，不能用 MAE 排名来认定原论文采用了哪种噪声定义。失锁及整包损失仍须与下文一起看。

一个**待验证解释**是定时参考不同：本接收器按收到的 HFM 相关图中首个显著路径确定载荷窗，而平均 PDP 峰、相关最强路径及 oracle 的参考并不必然指向同一时刻。相同时移占 4 ms 符号的比例更大，可能使短符号结果更敏感；当前数据尚未隔离这一因素，不能断言它已解释上述偏差。还需要在相同接收波形上比较这些参考的定时偏移及误码关联，才能判断其贡献。本次没有为此修改接收器或选择新的频率。

## 本次样本和接收链路

- {len(data)} 个包级观测；完整运行应为 96 个条件，每条件 3 条记录 × 5 包 × 64 bit = 960 个发送载荷 bit。
- 每包在 6.4 s 完整监听槽内具有 0.12–0.34 s 随机未知起点；所有候选共用载荷、起点、随机载波相位及标准正态噪声向量。符号相位重置为该包相位。
- `replayfilter_fast` 后按官方 `alpha=1/(1-V0/1500)` **逐槽恢复整体时间尺度**。每个槽的局部时间原点独立，因此本对齐诊断不声称复现整段 32 s 连续传播的累计漂移。接收器只获取收到的波形、公开帧配置和载荷长度；真实 TVIR、V0、起点、载荷不输入接收器。
- HFM 包含 64 ms UP、64 ms 间隔及 64 ms DOWN，扫频 10–18 kHz；另加两个 64 ms 已知双音导频和保护段，载荷标称偏移 528 ms。后者和频带/加窗是本项目选择。
- 接收器必需自检 {metadata['receiver_required_passed']}/71 PASS。这里没有按诊断 BER 调捕获阈值或调发射频率。

{recovery_note}

## 两种噪声定义及结果

令 `E=sum(clean_slot.^2)`，包括前导、导频、载荷和信道拖尾；`L=ceil(alpha*(完整发送帧时长+128 ms尾段)*fs)`，排除随机监听前后背景。

1. `sample_snr_packet`：`sigma²=E/L/10^(label/10)`，含前导、保护和尾段的整帧平均功率。保护时间改变 L，因此也会改变噪声；它不是固定噪声的纯保护时间对照。
2. `ebn0_received`：`sigma²=E/(2*64*10^(label/10))`，整帧接收能量按 64 个信息 bit 分摊，明确计入前导和导频能量。这与按载荷有效能量定义的 Eb/N0 也不同。

同一条件的两个方差满足 `sigma²_EbN0/sigma²_sample=L/(2*64)`，本报告脚本已逐包验证。**原论文仅给 SNR 标签；下面的差异不能用于反推其未披露的噪声实现。**

{md_table(shown)}

条件 BER 仅以实际得到判决的载荷 bit 为分母；捕获失败仍计入包损失并记录 64 bit 擦除，不把擦除虚构成误码或零误码。共有 {unavailable} 个条件完全没有已判决 bit，其条件 BER 记 NA。`capture_loss` 指接收器未输出完整载荷，并未另行通过真值验证多径中的定时关联正确率；已返回但错误的载荷进入误码和包损失。

观察到 {len(zeros)} 个有判决且零比特错误的条件；这不证明真实 BER=0。每条件最多只有 960 bit，单个错误的分辨率为 1/960≈0.00104，远不足以确认原表 `<0.0001` 的水平。位、包、相邻记录相关，因此没有把 960 bit 当成独立信道样本计算置信区间；`per_record_summary.csv` 保留 3 条记录各自结果。

{alignment_figures}

## FFT bin 对齐诊断

单独的 AWGN 实验在同一接收采样波形上比较精确双音能量与无补零 N 点矩形窗 FFT 的最近一个频点；最近 bin 半整数按 MATLAB `round` 向上取整，判决能量相等时选 0。不是原文未披露的窄带积分实现。三组频率事先固定，所有频率组共用比特、逐符号随机相位和噪声实现。每个配置 {int(fft.bits.min())} bit。

{md_table(fft_shown)}

13,875/14,125 Hz 在 4 ms 时为 55.5/56.5 个周期，在 8 ms 为 111/113 整数周期；较低、较高两组在这两种时长均落在整 bin。4 ms 下最近-bin 对称频率可能把一个音的能量分散到两个候选 bin，结果还取决于半-bin 舍入与有限窗相位；精确频率相关不受同一种频点舍入限制。补零可以细化采样网格，但不增加观察时长决定的频率分辨率。

![AWGN FFT 对齐诊断]({(out/"03_fft_alignment_awgn.png").as_posix()})

这个实验隔离了检测器与 bin 位置关系；若在 Watermark 实测信道中改发射频率，也会改变频率选择性衰落，不能把变化全部归给 FFT。未在留出数据上挑频率，也不根据与原表接近的程度选择噪声定义。

## 原文证据边界

[2024 主论文](https://vibsys.put.poznan.pl/_journal/2024-35-2/articles/vps_2024207.pdf)介绍 HFM-UP/HFM-DOWN 前导和 FFT/窄带判决；第 6 页将 BFSK 的 FFT 输入点数 N 与 Ts 内的采样点数关联，FFH 则与 Th 内采样关联，不能说其完全没有交代点数关系。它仍没有给出实际采样率、具体窗型、补零设置、窄频带边界及完整频带统计实现，也没有披露 HFM 具体时长、双音绝对频率、逐符号相位策略、完整捕获设置或噪声归一化。[2023 引用同步论文](https://vibsys.put.poznan.pl/_journal/2023-34-1/articles/vps_2023102.pdf)提供 64 ms 脉冲与 64 ms 间隔示例，不代表主论文采用了相同参数；其印刷的多普勒表达式存在量纲/零速解释疑点，本项目没有直接照抄，而以公开导频估计时间尺度并通过正负速度合成测试核验。详见 `docs/protocol_v2_review.md` 和 `docs/paper_audit.md`。

完整 48 行原表在两种噪声定义下分别对照，写入 `paper_comparison_all_conditions.csv`，原文三条 `<0.0001` 原样保留且不计算精确差值。下面全量列出；所有结果均是本次已披露假设下的实验观察。

{md_table(all_rows)}
'''
    (out/'report.md').write_text(report, encoding='utf-8')
    files = [ROOT/'matlab/run_paper_alignment_v2.m', Path(__file__), ROOT/'matlab/+uac2/receive_packet.m',
             ROOT/'matlab/+uac2/default_config.m', ROOT/'references/paper_bfsk_reference.csv', out/'per_packet.csv']
    status = {'status': 'PASS', 'full_run': is_full, 'packet_rows': len(data), 'condition_rows': len(summary),
              'paper_comparison_rows': len(compared), 'heldout_records_used': [],
              'noise_variance_ratio_verified': True,
              'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    (out/'report_validation.json').write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in status.items() if k != 'sha256'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
