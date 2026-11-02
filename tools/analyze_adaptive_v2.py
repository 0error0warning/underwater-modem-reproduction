"""Fit on development records, freeze, then evaluate a causal guard policy."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/'docs/adaptive_protocol_v2.md'
GUARDS = [0, 4, 8, 16, 32]
FEATURES = ['delay_spread_ms', 'sync_score', 'abs_velocity_hat']
CONTROL_S = 0.16
PRIMARY_SNR = 20.0
NOISE_MODE = 'fixed_noise_active_reference'
RECEIVER = 'sync_compensated'
REQUIRED = ['channel', 'record', 'slot', 'on_ms', 'guard_ms', 'snr_db',
            'noise_mode', 'receiver_mode', 'detected', 'payload_bits',
            'bit_errors', 'bit_erasures', 'packet_error', 'occupied_s',
            'slot_s', 'delay_spread_ms', 'sync_score', 'velocity_hat']


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def boolean(series):
    values = series.astype(str).str.strip().str.lower()
    if not values.isin(['true', 'false', '1', '0', '1.0', '0.0']).all():
        raise ValueError('Unrecognized detected flag')
    return values.isin(['true', '1', '1.0'])


def load_records(path, stage):
    data = pd.read_csv(path)
    missing = set(REQUIRED)-set(data.columns)
    if missing:
        raise ValueError(f'Missing columns: {sorted(missing)}')
    numeric = [c for c in REQUIRED if c not in ['channel', 'noise_mode', 'receiver_mode', 'detected']]
    for col in numeric:
        data[col] = pd.to_numeric(data[col], errors='raise')
    data['detected'] = boolean(data.detected)
    allowed = set(range(1, 11)) if stage == 'fit' else set(range(31, 41))
    if not set(data.record.astype(int)) <= allowed:
        raise ValueError(f'{stage} input contains records outside {sorted(allowed)}. Keep development/test files separate.')
    data = data[(data.noise_mode == NOISE_MODE) & (data.receiver_mode == RECEIVER) & (data.on_ms == 8)].copy()
    if 'scenario' in data.columns:
        data = data[data.scenario == 'natural'].copy()
    if data.empty:
        raise ValueError('No natural 8 ms sync_compensated fixed-reference-noise data')
    if not np.isfinite(data[['occupied_s', 'slot_s', 'payload_bits', 'bit_errors', 'bit_erasures', 'packet_error']].to_numpy()).all():
        raise ValueError('Non-finite timing/error count')
    if (data.occupied_s <= 0).any() or (data.slot_s <= 0).any():
        raise ValueError('Non-positive time')
    if not data.packet_error.isin([0, 1]).all():
        raise ValueError('packet_error must be 0/1')
    if ((data.bit_errors < 0) | (data.bit_erasures < 0) | (data.bit_errors+data.bit_erasures > data.payload_bits)).any():
        raise ValueError('Invalid error/erasure count')
    failed = ~data.detected
    if ((data.loc[failed, 'packet_error'] != 1) | (data.loc[failed, 'bit_erasures'] != data.loc[failed, 'payload_bits'])).any():
        raise ValueError('Capture failure must remain a packet failure with erased payload')
    if ((data.packet_error == 0) & ((data.bit_errors != 0) | (data.bit_erasures != 0))).any():
        raise ValueError('Error-free packet flag conflicts with bit counts')
    data['abs_velocity_hat'] = data.velocity_hat.abs()
    extra = [c for c in ['trial', 'replicate'] if c in data.columns]
    keys = ['channel', 'record', 'snr_db', 'noise_mode']+extra+['slot']
    trajectory_sizes=data.groupby(keys[:-1],dropna=False).slot.nunique()
    if not (trajectory_sizes==5).all():
        raise ValueError('The locked paired protocol requires five complete slots per trajectory')
    if data.duplicated(keys+['guard_ms']).any():
        raise ValueError('Duplicate guard outcome within event')
    for key, group in data.groupby(keys, dropna=False):
        if set(group.guard_ms.astype(int)) != set(GUARDS):
            raise ValueError(f'Incomplete guard options at {key}')
        if group.payload_bits.nunique() != 1 or group.slot_s.nunique() != 1:
            raise ValueError(f'Unmatched payload size/measurement time at {key}')
    if set(data.channel) != {'NOF1', 'NCS1'}:
        raise ValueError('Expected both measured channels')
    for channel in ['NOF1', 'NCS1']:
        if set(data.loc[data.channel == channel, 'record'].astype(int)) != allowed:
            raise ValueError(f'Incomplete {stage} record set for {channel}')
        if set(data.loc[(data.channel == channel) & (data.snr_db == PRIMARY_SNR), 'record'].astype(int)) != allowed:
            raise ValueError(f'Incomplete primary-SNR {stage} record set for {channel}')
    return data, extra


def trajectories(data, extra):
    groups = ['channel', 'record', 'snr_db', 'noise_mode']+extra
    output = []
    for key, trajectory in data.groupby(groups, sort=True, dropna=False):
        events = []
        for slot, candidates in trajectory.groupby('slot', sort=True):
            options = {int(row['guard_ms']): row for row in candidates.to_dict('records')}
            events.append(options)
        output.append(events)
    return output


def fixed(guard):
    return {'kind': 'fixed', 'guard_ms': int(guard)}


def choose(model, previous):
    if model['kind'] == 'fixed':
        return model['guard_ms'], 'fixed', np.nan
    if previous is None or not previous['detected']:
        return model['fallback_guard_ms'], 'fallback', np.nan
    value = previous[model['feature']]
    if not np.isfinite(value):
        return model['fallback_guard_ms'], 'fallback', np.nan
    if value <= model['threshold']:
        return model['left_guard_ms'], 'left', value
    return model['right_guard_ms'], 'right', value


def replay(traces, model, policy_name='', extra_control_s=0.0, retain=True):
    selected = []
    good_bits = total_time = total_slot = 0.0
    branch_counts = {'left': 0, 'right': 0, 'fallback': 0, 'fixed': 0}
    for events in traces:
        previous = None
        for options in events:
            guard, branch, feature_value = choose(model, previous)
            row = options[guard]
            branch_counts[branch] += 1
            success_bits = row['payload_bits'] if row['packet_error'] == 0 else 0
            duration = row['occupied_s']+CONTROL_S+extra_control_s
            good_bits += success_bits
            total_time += duration
            total_slot += row['slot_s']
            if retain:
                result = dict(row)
                result.update(policy=policy_name, selected_guard_ms=guard,
                              decision_branch=branch, previous_feature=feature_value,
                              previous_slot=np.nan if previous is None else previous['slot'],
                              success_bits=success_bits, accounted_s=duration,
                              base_control_s=CONTROL_S, extra_control_s=extra_control_s)
                selected.append(result)
            previous = row
    score = good_bits/total_time if total_time else 0.0
    return score, selected, branch_counts, good_bits/total_slot if total_slot else 0.0


def fit_stage(data, extra, source, out, force=False):
    policy_path = out/'frozen_policy.json'
    if policy_path.exists():
        if (out/'evaluation.json').exists():
            raise RuntimeError('Final evaluation already exists; refusing to refit after test exposure.')
        if not force:
            raise RuntimeError('Frozen policy exists. Use --force-refit only before final evaluation if a development bug was fixed.')
    primary = data[data.snr_db == PRIMARY_SNR]
    train = primary[primary.record.between(1, 5)]
    valid = primary[primary.record.between(6, 10)]
    if train.empty or valid.empty:
        raise ValueError('Missing 20 dB training/validation data')
    train_traces = trajectories(train, extra)
    valid_traces = trajectories(valid, extra)
    fixed_train = [(replay(train_traces, fixed(g), retain=False)[0], g) for g in GUARDS]
    fallback = sorted(fixed_train, key=lambda v: (-v[0], v[1]))[0][1]
    candidates = []
    for guard in GUARDS:
        model = fixed(guard)
        tr = replay(train_traces, model, retain=False)
        va = replay(valid_traces, model, retain=False)
        candidates.append({'name': f'fixed_{guard}', 'model': model,
                           'train_efficiency_bps': tr[0], 'validation_efficiency_bps': va[0],
                           'validation_slot_rate_bps': va[3]})
    searched = []
    for feature in FEATURES:
        values = train.loc[train.detected, feature].to_numpy()
        values = values[np.isfinite(values)]
        if not len(values):
            continue
        thresholds = np.unique(np.quantile(values, np.arange(.1, 1.0, .1)))
        feature_candidates = []
        for threshold in thresholds:
            for left in GUARDS:
                for right in GUARDS:
                    if left == right:
                        continue
                    model = {'kind': 'stump', 'feature': feature, 'threshold': float(threshold),
                             'left_guard_ms': left, 'right_guard_ms': right,
                             'fallback_guard_ms': fallback}
                    score, _, branches, slot_rate = replay(train_traces, model, retain=False)
                    if branches['left'] < 8 or branches['right'] < 8:
                        continue
                    feature_candidates.append((score, threshold, left, right, model, branches))
        searched.append({'feature': feature, 'threshold_count': len(thresholds), 'eligible_rules': len(feature_candidates)})
        if feature_candidates:
            best = sorted(feature_candidates, key=lambda v: (-v[0], v[1], v[2], v[3]))[0]
            va = replay(valid_traces, best[4], retain=False)
            candidates.append({'name': f'stump_{feature}', 'model': best[4],
                               'train_efficiency_bps': best[0], 'validation_efficiency_bps': va[0],
                               'validation_slot_rate_bps': va[3]})
    constants = [c for c in candidates if c['model']['kind'] == 'fixed']
    stumps = [c for c in candidates if c['model']['kind'] == 'stump']
    rank = lambda c: (-c['validation_efficiency_bps'], -c['train_efficiency_bps'], c['name'])
    best_fixed = sorted(constants, key=rank)[0]
    adaptive = sorted(stumps, key=rank)[0] if stumps else None
    selection = best_fixed
    if adaptive and adaptive['validation_efficiency_bps'] > best_fixed['validation_efficiency_bps']+1e-12:
        selection = adaptive
    model = {'status': 'FROZEN_BEFORE_FINAL_TEST', 'created_utc': utc(),
             'protocol_sha256': sha256(PROTOCOL), 'development_sha256': sha256(source),
             'development_source': str(source.resolve()), 'primary_snr_db': PRIMARY_SNR,
             'noise_mode': NOISE_MODE, 'receiver_mode': RECEIVER, 'guards_ms': GUARDS,
             'feature_names': FEATURES, 'fallback_guard_ms': fallback,
             'base_control_s': CONTROL_S, 'objective': 'success_bits / sum(occupied_s + base_control_s)',
             'train_records': list(range(1, 6)), 'validation_records': list(range(6, 11)),
             'final_records': list(range(31, 41)), 'extra_trajectory_keys': extra,
             'best_fixed': best_fixed['model'], 'adaptive_stump': adaptive['model'] if adaptive else None,
             'validation_selected': selection['model'], 'validation_selected_name': selection['name'],
             'candidate_scores': candidates, 'search_counts': searched,
             'causality': 'Previous selected action preamble only; reset at each record/condition; no extra probe'}
    write_json(policy_path, model)
    pd.DataFrame([{**{k:v for k,v in c.items() if k != 'model'}, 'model': json.dumps(c['model'])} for c in candidates]).to_csv(out/'development_scores.csv', index=False)
    note = f"冻结完成。训练回退 guard={fallback} ms；验证选中的 best_fixed={best_fixed['model']['guard_ms']} ms；选择={selection['name']}。尚未读取或评测 031–040。"
    (out/'fit_status.txt').write_text(note, encoding='utf-8')
    print(note)
    print(json.dumps({'best_fixed': model['best_fixed'], 'adaptive_stump': model['adaptive_stump'],
                      'validation_selected': model['validation_selected']}, ensure_ascii=False))


def oracle_replay(traces):
    # Fractional optimization on final outcomes, for an explicitly impossible
    # reference only. A greedy per-event ratio is not globally ratio-optimal.
    events = [event for trace in traces for event in trace]
    ratio = 0.0
    chosen = []
    for _ in range(100):
        chosen = []
        for options in events:
            values = []
            for guard in GUARDS:
                row = options[guard]
                bits = row['payload_bits'] if row['packet_error'] == 0 else 0
                duration = row['occupied_s']+CONTROL_S
                values.append((bits-ratio*duration, -duration, -guard, row))
            chosen.append(max(values, key=lambda v: v[:3])[3])
        good = sum(r['payload_bits'] if r['packet_error'] == 0 else 0 for r in chosen)
        duration = sum(r['occupied_s']+CONTROL_S for r in chosen)
        updated = good/duration
        if abs(updated-ratio) < 1e-12:
            break
        ratio = updated
    results = []
    for row in chosen:
        result = dict(row)
        result.update(policy='oracle_current_outcomes', selected_guard_ms=row['guard_ms'],
                      decision_branch='oracle', previous_feature=np.nan, previous_slot=np.nan,
                      success_bits=row['payload_bits'] if row['packet_error'] == 0 else 0,
                      accounted_s=row['occupied_s']+CONTROL_S, base_control_s=CONTROL_S, extra_control_s=0.)
        results.append(result)
    return results


def summarize(rows, group_columns):
    results = []
    for key, frame in rows.groupby(group_columns, dropna=False, sort=True):
        if not isinstance(key, tuple): key = (key,)
        result = dict(zip(group_columns, key))
        payload = frame.payload_bits.sum()
        erasures = frame.bit_erasures.sum()
        decoded = payload-erasures
        result.update(packets=len(frame), payload_bits=int(payload), bit_errors=int(frame.bit_errors.sum()),
                      bit_erasures=int(erasures), conditional_ber=frame.bit_errors.sum()/decoded if decoded else np.nan,
                      capture_loss_rate=1-frame.detected.mean(), packet_error_rate=frame.packet_error.mean(),
                      success_bits=int(frame.success_bits.sum()), occupied_plus_control_s=frame.accounted_s.sum(),
                      measurement_slot_s=frame.slot_s.sum(),
                      counterfactual_efficiency_bps=frame.success_bits.sum()/frame.accounted_s.sum(),
                      measurement_slot_rate_bps=frame.success_bits.sum()/frame.slot_s.sum())
        results.append(result)
    return pd.DataFrame(results)


def bootstrap_difference(primary_rows, reps=2000):
    rng = np.random.default_rng(20261102)
    totals = primary_rows.groupby(['policy', 'channel', 'record'])[['success_bits', 'accounted_s', 'slot_s']].sum()
    policies = [x for x in ['adaptive_stump', 'validation_selected'] if x in primary_rows.policy.unique()]
    channels = sorted(primary_rows.channel.unique())
    records = {c: sorted(primary_rows.loc[primary_rows.channel == c, 'record'].unique()) for c in channels}
    keys = [(c, int(r)) for c in channels for r in records[c]]
    offset=0
    draws=[]
    for channel in channels:
        count=len(records[channel])
        draws.append(offset+rng.integers(0,count,size=(reps,count)))
        offset+=count
    draws=np.concatenate(draws,axis=1)
    baseline=np.array([totals.loc[('best_fixed',c,r)].to_numpy() for c,r in keys],dtype=float)
    b=baseline[draws].sum(axis=1)
    output = []
    for policy in policies:
        values=np.array([totals.loc[(policy,c,r)].to_numpy() for c,r in keys],dtype=float)
        a=values[draws].sum(axis=1)
        diffs=a[:,0]/a[:,1]-b[:,0]/b[:,1]
        slot_diffs=a[:,0]/a[:,2]-b[:,0]/b[:,2]
        observed_a = primary_rows[primary_rows.policy == policy]
        observed_b = primary_rows[primary_rows.policy == 'best_fixed']
        output.append({'policy': policy, 'baseline': 'best_fixed', 'repeats': reps,
                       'delta_efficiency_bps': observed_a.success_bits.sum()/observed_a.accounted_s.sum()-observed_b.success_bits.sum()/observed_b.accounted_s.sum(),
                       'record_resample_efficiency_p025': float(np.quantile(diffs,.025)),
                       'record_resample_efficiency_p975': float(np.quantile(diffs,.975)),
                       'delta_slot_rate_bps': observed_a.success_bits.sum()/observed_a.slot_s.sum()-observed_b.success_bits.sum()/observed_b.slot_s.sum(),
                       'record_resample_slot_p025': float(np.quantile(slot_diffs,.025)),
                       'record_resample_slot_p975': float(np.quantile(slot_diffs,.975))})
    return pd.DataFrame(output)


def make_figures(summary, overhead, out):
    shown = summary[(summary.snr_db == PRIMARY_SNR) & (summary.policy != 'oracle_current_outcomes')]
    names = shown.policy.tolist()
    labels = [n.replace('fixed_', 'Fixed ').replace('_', ' ') for n in names]
    fig, axes = plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    axes[0].bar(np.arange(len(shown)), shown.counterfactual_efficiency_bps, color='#2378a8')
    axes[1].bar(np.arange(len(shown)), shown.measurement_slot_rate_bps, color='#e16a30')
    for ax in axes:
        ax.set_xticks(np.arange(len(shown)), labels, rotation=30, ha='right')
        ax.grid(axis='y',alpha=.2); ax.set_axisbelow(True); ax.set_ylabel('Correct payload bit/s')
    axes[0].set_title('Counterfactual occupied-time efficiency')
    axes[1].set_title('Common 6.4 s measurement-slot rate')
    fig.suptitle('Held-out records 31–40 | 20 dB fixed-reference noise')
    for suffix in ['png','svg']: fig.savefig(out/f'01_policy_comparison.{suffix}',dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.5,4.5),layout='constrained')
    for name,group in overhead.groupby('policy'):
        ax.plot(group.extra_control_s,group.counterfactual_efficiency_bps,'o-',label=name.replace('_',' '))
    ax.set(xlabel='Additional adaptive control / waiting time (s per packet)',ylabel='Counterfactual correct payload bit/s',
           title='Sensitivity to control overhead | held-out 20 dB')
    ax.grid(alpha=.2); ax.legend()
    for suffix in ['png','svg']: fig.savefig(out/f'02_control_overhead.{suffix}',dpi=180)
    plt.close(fig)


def format_model(model):
    if model is None: return '训练中没有满足最小叶子样本数的非恒定规则。'
    if model['kind'] == 'fixed': return f"始终选择 {model['guard_ms']} ms。"
    return (f"上一包 {model['feature']} ≤ {model['threshold']:.6g} 时选 {model['left_guard_ms']} ms，"
            f"否则选 {model['right_guard_ms']} ms；首包、上一包捕获失败或特征缺失时回退 {model['fallback_guard_ms']} ms。")


def report(model, summary, by_channel, intervals, overhead, out, metadata):
    main = summary[summary.snr_db == PRIMARY_SNR]
    table = ['|策略|反事实占用效率 bit/s|共同槽位载荷率 bit/s|条件 BER|捕获丢失率|包失败率|',
             '|---|---:|---:|---:|---:|---:|']
    for row in main.itertuples():
        ber = '无已解码载荷' if pd.isna(row.conditional_ber) else f'{row.conditional_ber:.3%}'
        table.append(f'|{row.policy}|{row.counterfactual_efficiency_bps:.3f}|{row.measurement_slot_rate_bps:.3f}|{ber}|{row.capture_loss_rate:.2%}|{row.packet_error_rate:.2%}|')
    channel_table=['|信道|策略|占用效率 bit/s|槽位载荷率 bit/s|包失败率|',
                   '|---|---|---:|---:|---:|']
    channel_main=by_channel[(by_channel.snr_db==PRIMARY_SNR)&by_channel.policy.isin(['best_fixed','adaptive_stump'])]
    for row in channel_main.itertuples():
        channel_table.append(f'|{row.channel}|{row.policy}|{row.counterfactual_efficiency_bps:.3f}|{row.measurement_slot_rate_bps:.3f}|{row.packet_error_rate:.2%}|')
    robustness=['|噪声标签 dB|策略|占用效率 bit/s|槽位载荷率 bit/s|包失败率|',
                '|---:|---|---:|---:|---:|']
    for row in summary[summary.policy.isin(['best_fixed','adaptive_stump'])].itertuples():
        robustness.append(f'|{row.snr_db:g}|{row.policy}|{row.counterfactual_efficiency_bps:.3f}|{row.measurement_slot_rate_bps:.3f}|{row.packet_error_rate:.2%}|')
    conclusions = []
    for row in intervals.itertuples():
        conclusions.append(f"- `{row.policy}` 相对开发集选定的 `best_fixed`，占用效率差值 {row.delta_efficiency_bps:+.3f} bit/s；按记录配对重采样的 2.5%–97.5% 范围为 [{row.record_resample_efficiency_p025:+.3f}, {row.record_resample_efficiency_p975:+.3f}] bit/s。共同槽位载荷率差值 {row.delta_slot_rate_bps:+.3f} bit/s。")
    adaptive = main[main.policy == 'adaptive_stump']
    baseline = main[main.policy == 'best_fixed']
    if len(adaptive):
        difference = adaptive.iloc[0].counterfactual_efficiency_bps-baseline.iloc[0].counterfactual_efficiency_bps
        lead = ('观察到自适应规则的反事实占用效率高于最佳固定策略' if difference > 0 else
                '本轮自适应规则没有超过开发集选定的最佳固定策略')
    else:
        lead = '训练数据没有产生满足预设样本要求的非恒定规则，本轮采用固定回退策略'
    text = f'''# 自适应保护时间 V2：冻结规则后的留出评测

在预设主条件 20 dB 下，{lead}。结论基于同一测量序列中预先留出的 031–040 号记录；不能据此声称任意海域上的泛化或真实反馈链路已完成。

## 已锁定的方法

- 有效发射时间 8 ms，保护时间候选 0/4/8/16/32 ms；主条件为 20 dB 的同记录固定噪声参考。
- 训练 001–005、验证 006–010；最终 031–040。冻结时间：{model['created_utc']}。
- 固定基线：{format_model(model['best_fixed'])}
- 最佳非恒定候选：{format_model(model['adaptive_stump'])}
- 验证选择：`{model['validation_selected_name']}`；{format_model(model['validation_selected'])}
- 使用上一包实际选择的保护时间对应的前导估计，没有免费的额外 probe，没有用当前包真值选择可部署策略。真实信道名称和 V0 不作为策略特征。

## 最终结果

{chr(10).join(table)}

![策略比较]({(out/'01_policy_comparison.png').resolve().as_posix()})

反事实占用效率分母为发射前导/导频、载荷、保护时间、包尾加统一 0.16 s 理想控制预算；共同槽位载荷率保留每个 6.4 s 实验槽位的全部时间。两者分别回答不同问题。规则选择基于前者，后者保证实验空闲没有被隐去。连续链路在选择不同保护时间后会经历不同信道时刻，需要另行连续时序验证。

`oracle_current_outcomes` 预先知道当前包所有候选的结果，仅作为不可实现参考。条件 BER 只统计已捕获并得到判决的载荷，捕获失败以擦除和包失败单独计入。无错包的判断来自仿真真值，未声称已实现 CRC。

### 分信道检查

{chr(10).join(channel_table)}

总体效率不能替代困难信道的可靠性：包失败率为 100% 表示这批包没有完整正确接收，不能因为同步已捕获就称通信已经可靠。

### 冻结模型的噪声鲁棒性

{chr(10).join(robustness)}

这些信噪比均使用相同冻结模型，未分别重新调参。

## 差值与开销

{chr(10).join(conclusions)}

重采样在两个信道内分别按整条记录抽取，两个策略使用相同抽样记录，重复 2,000 次。相邻记录仍可能相关，这些范围是记录重采样敏感性结果，不是跨环境保证。

![控制开销敏感性]({(out/'02_control_overhead.png').resolve().as_posix()})

所有方案已有 0.16 s/包的理想 ACK/模式通知预算；图中再仅给自适应增加 0、0.16、0.5、1.0 s。该图在相同包结果上重算时间分母，没有重新模拟多等待后的信道变化。这里没有测试回传物理链路，也没有假定这些数值覆盖任意水声传播距离。完整数值在 `control_overhead.csv`。

## 可核查结果

- `frozen_policy.json`：规则、候选评分、开发数据与协议哈希；拟合不读取最终文件。
- `development_scores.csv`：训练和验证候选表现。
- `selected_packets.csv`：每包动作、上一包特征、错误/擦除与时间预算。
- `test_summary.csv`、`test_by_channel.csv`、`test_by_record.csv`：总体、分信道与分记录结果。
- `paired_record_resampling.csv`：配对记录重采样差值。
- `evaluation.json`：最终数据源、哈希、冻结模型哈希和前导跨动作差异审计。

本报告固定使用冻结规则，不因最终结果而修改阈值或选取更有利的记录。两条测量序列中的估计速度、时延特征与信道身份可能相关，学到的阈值不能解释成普适的声学规律。若改进未成立，当前结果就是该预设策略在本轮条件下的结论。
'''
    (out/'report.md').write_text(text,encoding='utf-8')


def evaluate_stage(data, extra, source, out):
    policy_path=out/'frozen_policy.json'
    if not policy_path.exists(): raise RuntimeError('Fit/freeze first; final evaluation requires frozen_policy.json')
    model=json.loads(policy_path.read_text(encoding='utf-8'))
    if model.get('status') != 'FROZEN_BEFORE_FINAL_TEST' or model['protocol_sha256'] != sha256(PROTOCOL):
        raise RuntimeError('Frozen policy/protocol mismatch')
    if model['extra_trajectory_keys'] != extra: raise RuntimeError('Development/test trajectory structure differs')
    models = {'fixed_0':fixed(0),'fixed_8':fixed(8),'fixed_32':fixed(32),'best_fixed':model['best_fixed'],
              'validation_selected':model['validation_selected']}
    if model['adaptive_stump'] is not None: models['adaptive_stump']=model['adaptive_stump']
    traces=trajectories(data,extra)
    selected=[]
    for name,rule in models.items(): selected.extend(replay(traces,rule,name)[1])
    # Separate SNR objectives so the impossible oracle does not trade one
    # noise condition against another.
    for snr in sorted(data.snr_db.unique()):
        selected.extend(oracle_replay(trajectories(data[data.snr_db==snr],extra)))
    rows=pd.DataFrame(selected)
    rows.to_csv(out/'selected_packets.csv',index=False)
    summary=summarize(rows,['policy','snr_db'])
    by_channel=summarize(rows,['policy','snr_db','channel'])
    by_record=summarize(rows,['policy','snr_db','channel','record'])
    summary.to_csv(out/'test_summary.csv',index=False)
    by_channel.to_csv(out/'test_by_channel.csv',index=False)
    by_record.to_csv(out/'test_by_record.csv',index=False)
    primary_rows=rows[rows.snr_db==PRIMARY_SNR]
    intervals=bootstrap_difference(primary_rows)
    intervals.to_csv(out/'paired_record_resampling.csv',index=False)
    overhead=[]
    for name in ['best_fixed','adaptive_stump','validation_selected']:
        if name not in models: continue
        base=primary_rows[primary_rows.policy==name]
        is_adaptive=models[name]['kind']=='stump'
        for extra_time in [0.,.16,.5,1.]:
            overhead.append({'policy':name,'extra_control_s':extra_time,
                             'applied_extra_s':extra_time if is_adaptive else 0.,
                             'counterfactual_efficiency_bps':base.success_bits.sum()/(base.accounted_s.sum()+len(base)*(extra_time if is_adaptive else 0.))})
    overhead=pd.DataFrame(overhead)
    overhead.to_csv(out/'control_overhead.csv',index=False)
    event_keys=['channel','record','snr_db','slot']+extra
    spread={}
    for feature in FEATURES:
        ranges=data.groupby(event_keys)[feature].agg(lambda v: float(v.max()-v.min()) if np.isfinite(v).any() else np.nan)
        spread[feature]={'max_finite_cross_guard_range':float(ranges.max()) if ranges.notna().any() else None,
                         'events_with_finite_values':int(ranges.notna().sum())}
    audit={'status':'COMPLETE_FROZEN_POLICY_EVALUATION','created_utc':utc(),
           'test_source':str(source.resolve()),'test_sha256':sha256(source),
           'policy_sha256':sha256(policy_path),'protocol_sha256':sha256(PROTOCOL),
           'record_ids':sorted(int(v) for v in data.record.unique()),'preamble_feature_cross_guard_ranges':spread,
           'note':'Features may differ across guard outcomes; recursive evaluation uses only the preceding chosen action, never a free reference probe.',
           'bootstrap':'2000 paired whole-record resamples, stratified by channel; adjacent-record dependence remains a limitation'}
    write_json(out/'evaluation.json',audit)
    make_figures(summary,overhead,out)
    report(model,summary,by_channel,intervals,overhead,out,audit)
    print(summary[summary.snr_db==PRIMARY_SNR].to_string(index=False))
    print(f'Completed frozen-policy evaluation: {out}')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['fit','evaluate'],required=True)
    parser.add_argument('--records',required=True)
    parser.add_argument('--output',default='results/v2/adaptive')
    parser.add_argument('--force-refit',action='store_true')
    args=parser.parse_args()
    source=Path(args.records); out=Path(args.output)
    if not source.is_absolute(): source=ROOT/source
    if not out.is_absolute(): out=ROOT/out
    out.mkdir(parents=True,exist_ok=True)
    # Check the freeze BEFORE opening the held-out file.
    if args.stage=='evaluate':
        frozen=out/'frozen_policy.json'
        if not frozen.exists(): raise RuntimeError('Cannot open final data before fitting/freezing.')
        model=json.loads(frozen.read_text(encoding='utf-8'))
        if model['protocol_sha256']!=sha256(PROTOCOL): raise RuntimeError('Protocol changed since freeze.')
    data,extra=load_records(source,args.stage)
    if args.stage=='fit': fit_stage(data,extra,source,out,args.force_refit)
    else: evaluate_stage(data,extra,source,out)


if __name__=='__main__':
    main()
