"""Audit and report the phase-only bugfix rerun; never alter V2 artifacts."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from publication import verified_source_hash

from report_link_v2 import (
    BOOTSTRAP_REPEATS, BOOTSTRAP_SEED, continuous_per_record,
    continuous_summary, file_link, md_table, paired_continuous, read_json,
    require_complete, style,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/v3/continuous_phase_corrected"
V2 = ROOT / "results/v2"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def with_overall(per_record: pd.DataFrame) -> pd.DataFrame:
    return continuous_summary(pd.concat([
        per_record, per_record.assign(channel="ALL")
    ], ignore_index=True))


def verify_frozen_sources(preflight: dict) -> dict:
    policy_path = V2 / "adaptive/frozen_policy.json"
    policy_hash = verified_source_hash(policy_path, preflight["policy_sha256"])
    assert policy_hash == preflight["policy_sha256"] == preflight["expected_original_policy_sha256"]
    old_path = ROOT / "matlab/run_continuous_v2.m"
    new_path = ROOT / "matlab/run_continuous_phasecheck_v3.m"
    assert sha256(old_path) == preflight["old_driver_sha256"]
    assert sha256(new_path) == preflight["new_driver_sha256"]
    old = old_path.read_text(encoding="utf-8")
    new = new_path.read_text(encoding="utf-8")
    old_loop = old.split("names=", 1)[1].split("metadata.status=", 1)[0]
    new_loop = new.split("names=", 1)[1].split("metadata.analysis_type=", 1)[0]
    expected = old_loop.replace("for pi=", "for policyIndex=").replace("{pi}", "{policyIndex}")
    actual = new_loop.replace("                rows(rowCount).carrier_phase_rad=cfg.phase_rad;\n", "")
    assert expected == actual, "Computational loop changed beyond loop-variable rename and phase logging"
    for expression in preflight["random_expressions_unchanged"]:
        assert expression in new and expression in old
    receiver_checks = {}
    freeze = read_json(V2 / "receiver_code_freeze.json")
    for item in freeze["files"]:
        path = ROOT / item["path"]
        assert sha256(path) == item["sha256"], f"Frozen source changed: {path}"
        receiver_checks[item["path"]] = item["sha256"]
    for name, expected_hash in preflight["receiver_hashes"].items():
        assert sha256(ROOT / "matlab/+uac2" / name) == expected_hash
    preserved = {}
    report_notice_changes = {}
    old_receipt = read_json(V2 / "link_report_validation.json")
    for path_text, expected_hash in old_receipt["sha256"].items():
        path = ROOT / Path(path_text.replace("\\", "/"))
        actual_hash = sha256(path)
        if actual_hash != expected_hash and path == ROOT / "tools/report_link_v2.py":
            report_notice_changes[str(path.relative_to(ROOT))] = dict(
                old_sha256=expected_hash, current_sha256=actual_hash,
                reason="Parent task intentionally added a V2 continuous-result correction notice; receiver, data, driver and policy frozen hashes unaffected")
            continue
        verified_source_hash(path, expected_hash)
        preserved[str(path.relative_to(ROOT))] = expected_hash
    return dict(policy_sha256=policy_hash, frozen_sources=receiver_checks,
                preserved_v2_sources=preserved,
                intentional_report_notice_changes=report_notice_changes,
                computational_loop_only_phase_shadow_fix=True)


def verify_phases(data: pd.DataFrame) -> dict:
    assert np.isfinite(data.carrier_phase_rad).all()
    assert data.carrier_phase_rad.between(0, 2*np.pi, inclusive="left").all()
    assert data.carrier_phase_rad.max() > 4, "Unexpected retained V2 truncated phase range"
    wide = data.pivot(index=["channel", "record", "packet"], columns="policy", values="carrier_phase_rad")
    common = wide.dropna()
    delta = common.adaptive_stump - common.best_fixed
    assert (delta == 0).all(), "Common packet seed has different phase across policies"
    phases = data.groupby("policy").carrier_phase_rad.agg(["min", "max", "size"])
    return dict(common_packet_pairs=len(common), identical_phase_pairs=int((delta == 0).sum()),
                maximum_pair_phase_difference_rad=float(delta.abs().max()),
                phase_range_by_policy=phases.to_dict(orient="index"),
                unpaired_extra_packet_ids=len(wide)-len(common))


def verify_frozen_actions(data: pd.DataFrame, policy: dict) -> None:
    """Re-evaluate the frozen rule on each trajectory's previous observed packet."""
    for (name, _, _), group in data.groupby(["policy", "channel", "record"]):
        model = policy[name]
        previous = None
        for row in group.sort_values("packet").itertuples():
            if model["kind"] == "fixed":
                expected = model["guard_ms"]
            else:
                expected = model["fallback_guard_ms"]
                if previous is not None and previous.detected:
                    value = abs(previous.velocity_hat) if model["feature"] == "abs_velocity_hat" else getattr(previous, model["feature"])
                    if np.isfinite(value):
                        expected = model["left_guard_ms"] if value <= model["threshold"] else model["right_guard_ms"]
            assert row.guard_ms == expected
            previous = row


def plot_correction(comparison: pd.DataFrame) -> Path:
    style()
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True, sharey=True)
    for ax, channel in zip(axes, ["NOF1", "NCS1"]):
        rows = comparison[comparison.channel == channel].set_index("policy").loc[["best_fixed", "adaptive_stump"]]
        x = np.arange(2)
        ax.bar(x-.17, rows.full_horizon_bps_old, width=.34, color="#acb7be", label="V2 phase bug (superseded)")
        ax.bar(x+.17, rows.full_horizon_bps_corrected, width=.34, color="#16798b", label="V3 phase correction")
        for index, row in enumerate(rows.itertuples()):
            ax.text(index-.17, row.full_horizon_bps_old+.4, f"{row.full_horizon_bps_old:.1f}", ha="center", fontsize=9)
            ax.text(index+.17, row.full_horizon_bps_corrected+.4, f"{row.full_horizon_bps_corrected:.1f}", ha="center", fontsize=9)
        ax.set_xticks(x, ["Frozen fixed", "Frozen adaptive"])
        ax.set_title(channel)
        ax.set_ylabel("Zero-error payload / full 32 s (bit/s)")
        ax.set_ylim(0, max(3, comparison.full_horizon_bps_old.max(), comparison.full_horizon_bps_corrected.max())*1.18)
        ax.grid(axis="y", alpha=.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, -.08), ncol=2)
    fig.suptitle("Phase-only bugfix reanalysis | same revealed records 031-040\nNo receiver or policy refit | sample-power SNR 20 dB", fontsize=11)
    path = OUT / "phase_correction_comparison.png"
    fig.savefig(path, dpi=170, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> None:
    metadata = require_complete(OUT / "metadata_continuous.json")
    assert metadata["analysis_type"] == "BUGFIX_REANALYSIS" and metadata["not_new_holdout"]
    preflight = read_json(OUT / "preflight.json")
    frozen = verify_frozen_sources(preflight)
    policy = read_json(V2 / "adaptive/frozen_policy.json")
    data = pd.read_csv(OUT / "records_continuous.csv")
    assert set(data.snr_db) == {20} and (data.payload_bits == 128).all()
    assert np.allclose(data.control_s, policy["base_control_s"])
    assert (data.groupby(["channel", "record"]).noise_variance.nunique() == 1).all()
    phases = verify_phases(data)
    verify_frozen_actions(data, policy)
    records = continuous_per_record(data)
    expected = {(p, c, r) for p in ["best_fixed", "adaptive_stump"] for c in ["NOF1", "NCS1"] for r in range(31, 41)}
    assert set(records[["policy", "channel", "record"]].itertuples(index=False, name=None)) == expected
    summary = with_overall(records)
    paired, resampling = paired_continuous(records)
    old_records = continuous_per_record(pd.read_csv(V2 / "records_continuous.csv"))
    old_summary = with_overall(old_records)
    comparison = old_summary.merge(summary, on=["policy", "channel"], suffixes=("_old", "_corrected"), validate="one_to_one")
    comparison["correction_delta_bps"] = comparison.full_horizon_bps_corrected - comparison.full_horizon_bps_old
    compare_columns = ["policy", "channel", "packets_old", "packets_corrected", "zero_error_packets_old", "zero_error_packets_corrected", "full_horizon_bps_old", "full_horizon_bps_corrected", "correction_delta_bps"]
    comparison = comparison[compare_columns]
    for name, table in {"summary.csv": summary, "per_record.csv": records,
                        "paired_records.csv": paired, "paired_record_resampling.csv": resampling,
                        "v2_vs_corrected.csv": comparison}.items():
        table.to_csv(OUT / name, index=False)
    figure = plot_correction(comparison)
    shown = summary[["policy", "channel", "records", "packets", "zero_error_packets", "conditional_ber", "per", "horizon_s", "full_horizon_bps"]]
    difference_shown = resampling[["channel", "records", "fixed_mean_bps", "adaptive_mean_bps", "mean_difference_bps", "resample_p025_bps", "resample_p975_bps"]]
    nof1 = resampling[resampling.channel == "NOF1"].iloc[0]
    report = f"""# 连续排程初相错误：修正复跑报告

状态：**BUGFIX_REANALYSIS / COMPLETE**。这是原已查看记录 NOF1/NCS1 031–040 的错误修正重分析，**不是新留出测试**。共 {len(data)} 个发送帧、40 个策略/信道/记录单元；每单元均使用完整 32 s 分母。

修正后 NOF1 固定策略为 **{nof1.fixed_mean_bps:.1f} bit/s**，冻结自适应策略为 **{nof1.adaptive_mean_bps:.1f} bit/s**，记录配对均值差（自适应 − 固定）为 **{nof1.mean_difference_bps:+.1f} bit/s**。这批记录没有显示自适应优势；原开发验证选中的方案仍为 `{policy['validation_selected_name']}`。NCS1 的两策略均未得到整包无错载荷。

## 错误与修复范围

原 `run_continuous_v2.m` 中 `for pi=1:numel(models)` 把数学常量 `pi` 遮蔽。随后 `cfg.phase_rad=2*pi*rand()` 使固定/自适应初相分别落在 [0,2) 与 [0,4) rad，同一 packet 种子也得到不同初相，破坏了原声称的公平随机相位对照。

新驱动仅将策略循环变量改为 `policyIndex`，另增加 `carrier_phase_rad` 日志与修正标识。随机种子、随机抽样顺序、比特、背景窗口生成、公共绝对时间噪声、反馈预算、接收器和冻结策略均未改变。已用代码正文比较验证：计算主循环只存在变量重命名与初相字段记录差异。

策略 SHA-256：`{frozen['policy_sha256']}`。9 个 uac2 文件及原两个 V2 驱动均与原冻结清单一致；V2 既有原始数据和元数据也通过原报告 SHA-256 清单复核。主任务有意为旧报告及生成脚本添加连续结果更正提示，因此旧统计脚本整文件哈希变化另行记录；未改统计函数。原 V2 数据保留，连续部分数值应由本报告替代引用。

## 修正结果与原值对照

{md_table(comparison)}

![初相修正前后对照](<{figure.resolve().as_posix()}>)

初相修正会改变接收特征和后续策略选择，因此自适应轨迹的包数与时间可能随之变化。这些是相同冻结规则对修正波形的响应，并非重新拟合或挑选策略。

## 完整分母与记录配对统计

{md_table(shown)}

无错载荷速率 = 完全无错包的载荷 bit / 完整记录时长；每个策略在每信道使用 10 × 32 = 320 s，总体各 640 s。启动、保护间隔、公共尾段、0.16 s 反馈预算、时间网格取整和末尾空闲均在分母中。条件 BER 只以完整判决载荷为分母；捕获失败按整包擦除和包错误计入。

{md_table(difference_shown)}

以记录作配对单位，种子 {BOOTSTRAP_SEED}、{BOOTSTRAP_REPEATS} 次重抽样；总体保持各信道权重。2.5%–97.5% 分位范围仅描述当前记录集合，记录之间可能相关，不能解释为独立海况总体置信保证。不同策略的发包时刻不同，因此不把同 packet 编号视为相同物理信道事件的配对性能样本。

## 初相及因果策略核验

- 两策略共有 **{phases['common_packet_pairs']}** 个 `(channel, record, packet)` 键，初相全部完全一致，最大差 **{phases['maximum_pair_phase_difference_rad']:.1f} rad**；额外单侧 packet 键 {phases['unpaired_extra_packet_ids']} 个。
- 全部初相有限且属于 [0, 2π)，观测最大值超过 4 rad，已排除旧的截断范围。
- 每包实际 guard 均重新用同一冻结规则、该轨迹前一包已记录的接收特征和检测标记核验；首包/漏捕按回退规则处理。
- 噪声方差在同一信道记录的两策略间固定；连续排程无重叠，全部预算末端不超过 32 s。

数据与验证：{file_link('逐包原始记录', OUT / 'records_continuous.csv')}、{file_link('逐记录分母', OUT / 'per_record.csv')}、{file_link('配对记录差', OUT / 'paired_records.csv')}、{file_link('核验 JSON', OUT / 'correction_validation.json')}。

反馈仍是假设 0.16 s 内成功的理想一包延迟模型，没有模拟返回链路传播、可靠性或 ARQ/CRC。无错性由离线载荷真值计数，不能称为已实现确认吞吐量。这里的 20 dB 是固定参考接收功率的样本 SNR，不能与同数值的 Eb/N0 等同。此次未更改任何接收器参数，也未使用本次结果优化策略。
"""
    (OUT / "correction_report.md").write_text(report, encoding="utf-8")
    receipt = dict(status="PASS", analysis_type="BUGFIX_REANALYSIS", not_new_holdout=True,
                   packets=len(data), record_policy_units=len(records),
                   full_horizon_s_per_unit=32, full_horizon_denominator_verified=True,
                   frozen_causal_policy_actions_verified=True,
                   frozen=frozen, phase=phases,
                   resampling=dict(repeats=BOOTSTRAP_REPEATS, seed=BOOTSTRAP_SEED, unit="paired record within channel"),
                   sha256={str(p.relative_to(ROOT)): sha256(p) for p in [
                       OUT / "records_continuous.csv", OUT / "metadata_continuous.json",
                       OUT / "preflight.json", Path(__file__), ROOT / "matlab/run_continuous_phasecheck_v3.m",
                   ]})
    (OUT / "correction_validation.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False), encoding="utf-8")
    print(summary[["policy", "channel", "packets", "zero_error_packets", "full_horizon_bps"]].to_string(index=False))
    print(resampling.to_string(index=False))
    print(json.dumps({"status": receipt["status"], "phase": phases}, ensure_ascii=False))


if __name__ == "__main__":
    main()
