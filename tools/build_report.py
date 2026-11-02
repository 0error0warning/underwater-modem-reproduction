"""Build figures and a Chinese report only from a completed finite-scope run."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEYS = ['channel', 'on_ms', 'guard_ms', 'noise_mode', 'level_db']
MODES = ['ebn0_received', 'sample_snr_packet', 'fixed_noise_active_reference']
MODE_TITLES = {
    MODES[0]: 'Received Eb/N0',
    MODES[1]: 'Whole-packet sample SNR',
    MODES[2]: 'Fixed reference noise',
}
MODE_ZH = {
    MODES[0]: '接收端 Eb/N0',
    MODES[1]: '整包采样 SNR',
    MODES[2]: '同记录固定噪声',
}
COLORS = ['#1768a6', '#d75a27']


def load_completed(run_dir: Path):
    metadata = json.loads((run_dir/'run_metadata.json').read_text(encoding='utf-8'))
    if metadata.get('status') != 'COMPLETE_LIMITED_SCOPE':
        raise RuntimeError('Run is not COMPLETE_LIMITED_SCOPE; refusing to create a completion report.')
    data = pd.read_csv(run_dir/'summary.csv')
    per = pd.read_csv(run_dir/'per_sounding.csv')
    reference = pd.read_csv(ROOT/'references'/'paper_bfsk_reference.csv', dtype={'ber': str})
    n = int(metadata['n_soundings_per_channel'])
    if len(data) != 144 or len(per) != 144*n or len(reference) != 48:
        raise ValueError('Expected 144 summaries, 144*n sounding rows, and 48 paper references.')
    if data.duplicated(KEYS).any() or per.duplicated(KEYS+['sounding']).any():
        raise ValueError('Duplicate experiment keys.')
    for channel in ['NOF1', 'NCS1']:
        if set(per.loc[per.channel == channel, 'sounding']) != set(range(1, n+1)):
            raise ValueError(f'Unexpected sounding selection for {channel}.')
    totals = per.groupby(KEYS, as_index=False).agg(
        checked_bits=('bits', 'sum'), checked_errors=('errors', 'sum'),
        checked_soundings=('sounding', 'nunique'),
        checked_packets=('packets', 'sum'), checked_packet_errors=('packet_errors', 'sum'))
    checked = data.merge(totals, on=KEYS, validate='one_to_one')
    if len(checked) != 144:
        raise ValueError('Summary/per-sounding configurations differ.')
    for field in ['bits', 'errors', 'soundings']:
        if not np.array_equal(checked[field], checked['checked_'+field]):
            raise ValueError(f'Inconsistent summary {field}.')
    if not np.allclose(checked.ber, checked.checked_errors/checked.checked_bits, atol=1e-12):
        raise ValueError('Inconsistent BER aggregation.')
    if not np.allclose(checked.per, checked.checked_packet_errors/checked.checked_packets, atol=1e-12):
        raise ValueError('Inconsistent PER aggregation.')
    return metadata, data.sort_values(KEYS), per, reference


def style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.titlesize': 12, 'axes.labelsize': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.alpha': .2, 'grid.linewidth': .7,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
        'svg.fonttype': 'none',
    })


def save_figure(fig, destination: Path, stem: str):
    fig.savefig(destination/f'{stem}.png', dpi=190, bbox_inches='tight')
    fig.savefig(destination/f'{stem}.svg', bbox_inches='tight')
    plt.close(fig)


def mark_zero_observations(ax, rows: pd.DataFrame, x_column: str):
    for _, row in rows.loc[rows.errors == 0].iterrows():
        ax.annotate(f'0/{int(row.bits):,} obs.', (row[x_column], 0),
                    xytext=(4, 7), textcoords='offset points', fontsize=7,
                    color='#555555', rotation=15)


def figures(data: pd.DataFrame, destination: Path, n: int):
    destination.mkdir(exist_ok=True)
    style()
    fixed = data.loc[(data.noise_mode == MODES[2]) & (data.level_db == 20)]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5), layout='constrained')
    for ax, channel in zip(axes, ['NOF1', 'NCS1']):
        for guard, color, marker in zip([False, True], COLORS, ['o', 's']):
            rows = fixed.loc[(fixed.channel == channel) & ((fixed.guard_ms > 0) == guard)].sort_values('on_ms')
            ax.fill_between(rows.on_ms, rows.min_sounding_ber, rows.max_sounding_ber, color=color, alpha=.10)
            ax.plot(rows.on_ms, rows.ber, color=color, marker=marker, linewidth=1.7,
                    label='Guard = active time' if guard else 'No guard')
            mark_zero_observations(ax, rows, 'on_ms')
        ax.set(title=channel, xlabel='Active symbol duration (ms)', ylabel='Observed BER')
        ax.set_xticks([4, 8, 16, 32])
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.legend(frameon=False)
    fig.suptitle('Guard-time comparison | fixed reference noise, 20 dB', fontsize=14)
    fig.supxlabel(f'First {n}/60 recordings per channel. Shading: observed recording range, not a confidence interval.', fontsize=9)
    save_figure(fig, destination, '01_guard_time_ber')

    eight = data.loc[data.on_ms == 8]
    ymax = min(1, max(.02, float(eight.ber.max())*1.13))
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.7), sharex=True, sharey=True, layout='constrained')
    for ri, channel in enumerate(['NOF1', 'NCS1']):
        for mi, mode in enumerate(MODES):
            ax = axes[ri, mi]
            for guard, color, marker in zip([False, True], COLORS, ['o', 's']):
                rows = eight.loc[(eight.channel == channel) & (eight.noise_mode == mode) &
                                 ((eight.guard_ms > 0) == guard)].sort_values('level_db')
                ax.plot(rows.level_db, rows.ber, color=color, marker=marker, linewidth=1.7,
                        label='Guard = 8 ms' if guard else 'No guard')
                mark_zero_observations(ax, rows, 'level_db')
            ax.set_title(f'{channel} | {MODE_TITLES[mode]}', fontsize=10)
            ax.set_xticks([0, 10, 20])
            ax.set_ylim(0, ymax)
            ax.yaxis.set_major_formatter(PercentFormatter(1))
            if mi == 0:
                ax.set_ylabel('Observed BER')
            if ri == 1:
                ax.set_xlabel('Noise-level label (dB)')
            if ri == 0 and mi == 0:
                ax.legend(frameon=False, fontsize=8)
    fig.suptitle('Noise definitions produce different experiments | active time = 8 ms', fontsize=14)
    fig.supxlabel('The same dB label does not imply the same noise variance. All panels use the same BER scale.', fontsize=9)
    save_figure(fig, destination, '02_noise_conventions')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout='constrained')
    for ax, channel in zip(axes, ['NOF1', 'NCS1']):
        subset = fixed.loc[fixed.channel == channel]
        for active, rows in subset.groupby('on_ms'):
            rows = rows.sort_values('guard_ms')
            ax.plot(rows.error_free_packet_goodput, rows.ber, color='#b5bec5', linewidth=1, zorder=1)
            for _, row in rows.iterrows():
                is_guard = row.guard_ms > 0
                ax.scatter(row.error_free_packet_goodput, row.ber, s=54, marker='s' if is_guard else 'o',
                           color=COLORS[int(is_guard)], edgecolor='white', linewidth=.6, zorder=3)
                offset = (5, -11 if is_guard else 6)
                alignment = 'left'
                if channel == 'NCS1' and int(active) == 8 and is_guard:
                    offset = (9, 8)
                if channel == 'NCS1' and int(active) == 32 and not is_guard:
                    offset = (9, -12)
                if row.error_free_packet_goodput >= .95*subset.error_free_packet_goodput.max():
                    offset = (-8, 8)
                    alignment = 'right'
                ax.annotate(f'{int(active)} ms', (row.error_free_packet_goodput, row.ber),
                            xytext=offset, textcoords='offset points', fontsize=8, ha=alignment)
            mark_zero_observations(ax, rows, 'error_free_packet_goodput')
        for color, marker, label in zip(COLORS, ['o', 's'], ['No guard', 'Guard = active time']):
            ax.scatter([], [], color=color, marker=marker, label=label)
        ax.set(title=channel, xlabel='Ideal error-free payload rate (bit/s)', ylabel='Observed BER')
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=-max(1, subset.error_free_packet_goodput.max()*.04))
        ax.margins(x=.15, y=.15)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle('Reliability and useful payload rate | fixed reference noise, 20 dB', fontsize=14)
    fig.supxlabel('A failed 128-bit packet contributes no useful bits. No CRC, retransmission or acquisition overhead is modeled.', fontsize=9)
    save_figure(fig, destination, '03_ber_payload_tradeoff')


def comparison_csv(data: pd.DataFrame, reference: pd.DataFrame, run_dir: Path):
    ref = reference.rename(columns={'ber': 'paper_ber_text', 'table': 'paper_table',
                                    'page': 'paper_page', 'snr_label_db': 'paper_snr_label_db'})
    ref['paper_ber_exact'] = pd.to_numeric(ref.paper_ber_text, errors='coerce')
    ref['paper_ber_reported_upper_bound'] = ref.paper_ber_text.where(
        ref.paper_ber_text.str.startswith('<')).str.slice(1).astype(float)
    joined = data.merge(ref, left_on=['channel', 'on_ms', 'guard_ms', 'level_db'],
                        right_on=['channel', 'on_ms', 'guard_ms', 'paper_snr_label_db'],
                        validate='many_to_one', how='left')
    if len(joined) != 144 or joined.paper_ber_text.isna().any():
        raise ValueError('Paper reference join did not cover all 144 simulated conditions.')
    joined['raw_ber_difference_from_exact_paper_value'] = joined.ber-joined.paper_ber_exact
    joined['comparison_scope'] = 'Same numerical label; noise/receiver conditions are not established identical'
    joined = joined.rename(columns={
        'noise_mode': 'simulation_noise_definition', 'level_db': 'simulation_noise_level_db',
        'ber': 'simulation_observed_ber', 'bits': 'simulation_tested_bits',
        'errors': 'simulation_error_bits', 'per': 'simulation_observed_per',
        'error_free_packet_goodput': 'simulation_ideal_error_free_payload_bps',
    })
    ordered = ['channel', 'on_ms', 'guard_ms', 'simulation_noise_definition',
               'simulation_noise_level_db', 'paper_snr_label_db', 'paper_ber_text',
               'paper_ber_exact', 'paper_ber_reported_upper_bound', 'paper_table', 'paper_page',
               'simulation_error_bits', 'simulation_tested_bits', 'simulation_observed_ber',
               'raw_ber_difference_from_exact_paper_value', 'simulation_observed_per',
               'simulation_ideal_error_free_payload_bps', 'soundings', 'min_sounding_ber',
               'max_sounding_ber', 'noiseless_ber', 'nominal_bit_rate',
               'payload_rate_with_padding', 'comparison_scope']
    joined[ordered].to_csv(run_dir/'paper_comparison_all_conditions.csv', index=False, encoding='utf-8-sig')


def observed(row):
    errors, bits = int(row.errors), int(row.bits)
    if errors == 0:
        return f'0/{bits:,}（本次未观察到错误）'
    return f'{row.ber:.3%}（{errors:,}/{bits:,}）'


def make_report(metadata, data, per, reference, run_dir):
    n = int(metadata['n_soundings_per_channel'])
    fixed = data.loc[(data.noise_mode == MODES[2]) & (data.level_db == 20)]
    base = fixed.loc[fixed.guard_ms == 0].set_index(['channel', 'on_ms'])
    guarded = fixed.loc[fixed.guard_ms > 0].set_index(['channel', 'on_ms'])
    improved = int((guarded.ber < base.ber).sum())
    equal = int((guarded.ber == base.ber).sum())
    rate_up = int((guarded.error_free_packet_goodput > base.error_free_packet_goodput).sum())
    case_table = ['|信道|噪声定义，标签均为 20 dB|无保护 BER（错误数/比特数）|8 ms 保护 BER（错误数/比特数）|',
                  '|---|---|---|---|']
    for channel in ['NOF1', 'NCS1']:
        for mode in MODES:
            selected = data.loc[(data.channel == channel) & (data.on_ms == 8) &
                                (data.noise_mode == mode) & (data.level_db == 20)].set_index('guard_ms')
            case_table.append(f'|{channel}|{MODE_ZH[mode]}|{observed(selected.loc[0])}|{observed(selected.loc[8])}|')
    paper_table = ['|信道|原表 2：8 ms，无保护，SNR 标签 20 dB|原表 3：8 ms，8 ms 保护，SNR 标签 20 dB|',
                   '|---|---|---|']
    for channel in ['NOF1', 'NCS1']:
        selected = reference.loc[(reference.channel == channel) & (reference.on_ms == 8) &
                                 (reference.snr_label_db == 20)].set_index('guard_ms')
        paper_table.append(f'|{channel}|{selected.loc[0, "ber"]}|{selected.loc[8, "ber"]}|')
    observations = []
    for channel in ['NOF1', 'NCS1']:
        left, right = base.loc[(channel, 8)], guarded.loc[(channel, 8)]
        observations.append(f'- {channel} 的 8 ms 配置：无保护时 {observed(left)}，加入 8 ms 保护后 {observed(right)}；'
                            f'理想无错载荷率从 {left.error_free_packet_goodput:.2f} 变为 {right.error_free_packet_goodput:.2f} bit/s。')
    delta_sentences = []
    for channel in ['NOF1', 'NCS1']:
        selected = data.loc[(data.channel == channel) & (data.on_ms == 8) &
                            (data.guard_ms == 8) & (data.level_db == 0)].set_index('noise_mode')
        delta_sentences.append(f'- {channel}，8 ms 发射与 8 ms 保护，标签 0 dB：接收端 Eb/N0 口径为 '
                               f'{observed(selected.loc[MODES[0]])}；整包采样 SNR 口径为 '
                               f'{observed(selected.loc[MODES[1]])}；固定参考噪声口径为 {observed(selected.loc[MODES[2]])}。')
    zero_count = int((data.errors == 0).sum())
    zero_sentence = (f'本轮有 {zero_count} 个汇总配置未观察到错误。报告将此情况写成 `0/N`，不据此声称真实 BER 等于零。'
                     if zero_count else '本轮 144 个汇总配置均观察到错误。若某配置未观察到错误，应写成 `0/N`，不能据此声称真实 BER 等于零。')
    elapsed = float(metadata.get('elapsed_seconds', 0))
    image_root = (run_dir/'figures').as_posix()
    text = f'''# 水声 BFSK 保护时间：有限范围复现报告

生成时间：{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}。运行状态已核验为 `COMPLETE_LIMITED_SCOPE`。

本轮完成两种实测时变信道上的 BFSK 保护时间实验。固定噪声参考、标签 20 dB 的 8 组配置对比中，加入等长保护时间后观察 BER 降低的有 **{improved}/8** 组，相同的有 **{equal}/8** 组；理想无错载荷率提高的有 **{rate_up}/8** 组。误码改善与数据率改善应分别判断，不能仅凭 BER 降低认定吞吐量上升。

这是带明确时序及噪声假设的物理层实验，尚未精确复现原文所有数值，也不是完成了可部署的水声 Modem。

## 实验范围

- 论文：[Schmidt 等，2024，Performance of low-power underwater modem for shallow water communications](https://doi.org/10.21008/j.0860-6897.2024.2.07)，本轮对应表 2、表 3 的 BFSK 部分。
- 数据：Watermark 的 NOF1、NCS1，各取连续的前 **{n}/60** 个记录；不是全库，也不是独立抽取的海域样本。
- 波形：采样率 48 kHz，双音 13,875/14,125 Hz；有效发射时长 4/8/16/32 ms，静默保护为 0 或等于发射时长；128 bit/包；每符号相位重置为 0。
- 接收：已知包边界，以每条实测记录的平均功率时延谱峰值确定固定接收时延，再做双音正交匹配能量检测。没有自动捕获前导，也没有根据载荷真值搜索最低误码率的时延。
- 信道：保留 Watermark 时变冲激响应的残余变化；整体速度 `V0` 没有重新引入，因此不等于官方 `sfetch` 默认恢复整体多普勒后的完整链路。
- 实验量：144 个汇总配置，即 48 个原表标签配置乘以 3 种噪声定义；每配置包含 {n} 条信道记录。每配置实际测试比特数为 {int(data.bits.min()):,} 至 {int(data.bits.max()):,}，由包长和记录长度决定。回放及接收计算约 {elapsed:.1f} 秒，不包括下载、前置验证与报告生成。

## 保护时间与速率

下面的主比较在每条信道记录内固定噪声方差，参考值由同一记录的 8 ms、无保护波形能量按有效发射样本数定标。该噪声在本记录的不同发射时长、保护时间之间保持一致；不同记录会分别定标，不能理解为全实验共用一种环境噪声强度。

{chr(10).join(observations)}

![保护时间与观察 BER]({image_root}/01_guard_time_ber.png)

图中阴影是各记录观察 BER 的最小值到最大值，属于描述性范围，不是置信区间。不同配置按各自包长铺满记录，包时刻与最终覆盖位置不同；上述比较不属于严格逐包配对的因果实验。

![误码率与理想无错载荷率]({image_root}/03_ber_payload_tradeoff.png)

无错载荷率按“完全无误码的 128 bit 包所含信息量 / 所有已发送包的总时长”计算。分母包含符号保护、前端 2 ms 补零及信道尾部包间隔；未包括前导、CRC、确认、重传和收发切换开销。是否误码由仿真真值判定，所以该数值是理想无错载荷率，不是实装协议已实现的应用吞吐量。

## 噪声定义会改变实验

三种模式独立列出：

|模式|噪声方差定义|含义|
|---|---|---|
|接收端 Eb/N0|接收整包能量 / (2 × 信息比特数 × 线性 Eb/N0)|与官方 sfetch 的能量口径对应；本轮未引入其随机包偏移和整体多普勒恢复|
|整包采样 SNR|包含静默、补零、信道尾部的整包平均功率 / 线性 SNR|保护时间会影响平均功率，从而也影响噪声方差|
|同记录固定噪声|该记录 8 ms 无保护参考功率 / 线性标签值|用于同记录内比较保护时间与时长；不是跨记录统一噪声|

{chr(10).join(delta_sentences)}

![三种噪声定义的差异]({image_root}/02_noise_conventions.png)

所有曲线保留独立模式名称。相同的 0/10/20 dB 数字不代表相同的噪声强度，因此不能横向把三种模式视为同一个物理条件。

## 与原论文 8 ms 关键配置对照

原文表格值如下，均为作者结果：

{chr(10).join(paper_table)}

本轮结果如下：

{chr(10).join(case_table)}

论文只写 SNR，未充分交代噪声归一化、双音位置、FFT 检测窗、随机数据和具体信道记录。上表按相同数字标签列出，**不表示实验条件已经对齐**，也不根据接近原文的程度挑选噪声口径。接收时序、相位、检测器与信道处理等差异都可能影响数值。完整 144 行对照在 [paper_comparison_all_conditions.csv](paper_comparison_all_conditions.csv)，原文三条 `<0.0001` 仍按界限字符串保存，未替换为零。

## 统计与结论边界

{zero_sentence}所有图均使用线性观察 BER 轴，零错误点显示真实样本观察值并标注样本数，没有偷偷换成一个正数供对数作图。

相邻记录、同一包中的位错误和连续包可能相关；本轮未将它们当作独立伯努利样本给出虚假的精密置信区间。前 {n} 条连续记录不能证明对全部 60 条记录或其他水域普遍有效。低误码配置的样本量尤其限制了对原文 `<0.0001` 的判断。

当前可支持的结论是：在明确记录的波形、时序、信道片段和噪声定义下，保护时间改变了 BER 与有效速率的权衡，三种噪声定义也产生了不同结果。FFH、在线同步、真实多普勒补偿、低功耗 MCU 运行及水池通信均不属于本轮完成范围。

## 可核查文件

- [汇总数据](summary.csv)、[逐信道记录数据](per_sounding.csv)、[逐包数据](per_packet.csv)、[运行元数据](run_metadata.json)。
- [原文全部条件对照](paper_comparison_all_conditions.csv)；图的可编辑 SVG 与 PNG 同存于 `figures/`。
- 原始参考值：`references/paper_bfsk_reference.csv`；符号歧义与公式核验：`docs/paper_audit.md`。
'''
    (run_dir/'report.md').write_text(text, encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', default='results/watermark/first_10_soundings')
    args = parser.parse_args()
    run_dir = Path(args.run_dir)
    if not run_dir.is_absolute():
        run_dir = ROOT/run_dir
    metadata, data, per, reference = load_completed(run_dir)
    figures(data, run_dir/'figures', int(metadata['n_soundings_per_channel']))
    comparison_csv(data, reference, run_dir)
    # A concurrent restart must not result in a newly written completion report.
    current = json.loads((run_dir/'run_metadata.json').read_text(encoding='utf-8'))
    if current != metadata or current.get('status') != 'COMPLETE_LIMITED_SCOPE':
        raise RuntimeError('Run metadata changed during report generation; refusing completion report.')
    make_report(metadata, data, per, reference, run_dir)
    print(f'Validated 144 conditions against per-sounding totals; wrote report and 3 PNG/SVG figures: {run_dir}')


if __name__ == '__main__':
    main()
