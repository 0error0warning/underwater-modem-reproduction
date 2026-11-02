"""Report held-out receiver and genuinely sequential link experiments.

Formal generation refuses incomplete runs. --self-check exercises statistics
using completed smoke data and small synthetic records without publishing a
formal report or touching the adaptive-policy report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MODES = ["sync_compensated", "sync_uncompensated", "oracle_timing_scale"]
MODE_LABELS = {
    "sync_compensated": "Acquired + scale correction",
    "sync_uncompensated": "Same acquisition, no correction",
    "oracle_timing_scale": "Oracle diagnostic",
}
MODE_COLORS = {
    "sync_compensated": "#16798b",
    "sync_uncompensated": "#c76a37",
    "oracle_timing_scale": "#6859a1",
}
POLICY_LABELS = {"best_fixed": "Frozen fixed", "adaptive_stump": "Frozen adaptive"}
EVENT = ["channel", "record", "slot", "on_ms", "guard_ms", "snr_db", "noise_mode"]
BOOTSTRAP_REPEATS = 2000
BOOTSTRAP_SEED = 2026110207


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def as_list(value):
    return value if isinstance(value, list) else [value]


def require_complete(path: Path) -> dict:
    metadata = read_json(path)
    if metadata.get("status") != "COMPLETE":
        raise ValueError(f"Refusing a formal report from incomplete metadata: {path}")
    return metadata


def validate_counts(data: pd.DataFrame) -> None:
    columns = ["payload_bits", "bit_errors", "bit_erasures", "packet_error", "detected"]
    if not np.isfinite(data[columns].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite packet/error counts")
    assert (data.payload_bits > 0).all()
    assert data.detected.isin([0, 1]).all()
    assert data.packet_error.isin([0, 1]).all()
    assert (data.bit_errors >= 0).all() and (data.bit_erasures >= 0).all()
    assert (data.bit_errors <= data.payload_bits - data.bit_erasures).all()
    assert (data.bit_erasures == np.where(data.detected == 1, 0, data.payload_bits)).all()
    expected_error = ((data.detected == 0) | (data.bit_errors > 0)).astype(int)
    assert (data.packet_error == expected_error).all()


def validate_slot_data(data: pd.DataFrame, metadata: dict, stage: str) -> None:
    validate_counts(data)
    assert metadata.get("stage") == stage
    assert set(data.receiver_mode) == set(MODES)
    assert not data.duplicated(EVENT + ["receiver_mode"]).any()
    assert (data.groupby(EVENT, dropna=False).receiver_mode.nunique() == len(MODES)).all()
    assert set(data.channel) == set(metadata["channels"])
    assert set(data.record) == set(as_list(metadata["record_ids"]))
    assert set(data.guard_ms) == set(as_list(metadata["guards_ms"]))
    assert set(data.snr_db) == set(as_list(metadata["snr_db"]))
    assert np.allclose(data.slot_s, metadata["slot_s"])
    expected = (len(metadata["channels"]) * len(as_list(metadata["record_ids"]))
                * len(as_list(metadata["guards_ms"])) * metadata["slots_per_record"]
                * len(as_list(metadata["snr_db"])) * len(MODES))
    assert len(data) == expected, f"Incomplete {stage} packet grid: {len(data)} != {expected}"
    if stage == "test":
        assert set(data.record) == set(range(31, 41))
        assert set(data.snr_db) == {0, 10, 20}


def aggregate_packets(data: pd.DataFrame, keys: list[str], slot_denominator: bool = False) -> pd.DataFrame:
    work = data.copy()
    work["zero_error_payload_bits"] = np.where(work.packet_error == 0, work.payload_bits, 0)
    work["zero_error_packets"] = (work.packet_error == 0).astype(int)
    specs = dict(
        packets=("payload_bits", "size"), captured_packets=("detected", "sum"),
        payload_bits=("payload_bits", "sum"), bit_errors=("bit_errors", "sum"),
        bit_erasures=("bit_erasures", "sum"), packet_errors=("packet_error", "sum"),
        zero_error_packets=("zero_error_packets", "sum"),
        zero_error_payload_bits=("zero_error_payload_bits", "sum"),
    )
    if slot_denominator:
        specs["listening_seconds"] = ("slot_s", "sum")
    result = work.groupby(keys, sort=True, dropna=False).agg(**specs).reset_index()
    result["decoded_bits"] = result.payload_bits - result.bit_erasures
    result["conditional_ber"] = result.bit_errors / result.decoded_bits.replace(0, np.nan)
    result["capture_loss"] = 1 - result.captured_packets / result.packets
    result["per"] = result.packet_errors / result.packets
    if slot_denominator:
        result["slot_zero_error_bps"] = result.zero_error_payload_bits / result.listening_seconds
    return result


def paired_receivers(data: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    """Compare modes only on the exact same channel/slot/data/noise event."""
    wide = data.pivot(index=EVENT, columns="receiver_mode",
                      values=["packet_error", "detected", "bit_errors", "payload_bits"])
    assert not wide.isna().any().any()
    output = []
    comparisons = [(MODES[0], MODES[1]), (MODES[0], MODES[2]), (MODES[1], MODES[2])]
    for mode_a, mode_b in comparisons:
        event = wide.index.to_frame(index=False)
        good_a = wide[("packet_error", mode_a)].to_numpy() == 0
        good_b = wide[("packet_error", mode_b)].to_numpy() == 0
        cap_a = wide[("detected", mode_a)].to_numpy() == 1
        cap_b = wide[("detected", mode_b)].to_numpy() == 1
        common = cap_a & cap_b
        bits_a = wide[("payload_bits", mode_a)].to_numpy()
        assert np.array_equal(bits_a, wide[("payload_bits", mode_b)].to_numpy())
        event["a_only_zero_error"] = (good_a & ~good_b).astype(int)
        event["b_only_zero_error"] = (~good_a & good_b).astype(int)
        event["both_zero_error"] = (good_a & good_b).astype(int)
        event["neither_zero_error"] = (~good_a & ~good_b).astype(int)
        event["common_decoded_packets"] = common.astype(int)
        event["common_decoded_bits"] = np.where(common, bits_a, 0)
        event["a_errors_common"] = np.where(common, wide[("bit_errors", mode_a)].to_numpy(), 0)
        event["b_errors_common"] = np.where(common, wide[("bit_errors", mode_b)].to_numpy(), 0)
        event["events"] = 1
        numbers = [c for c in event if c not in EVENT]
        summary = event.groupby(groups, sort=True, dropna=False)[numbers].sum().reset_index()
        summary["mode_a"] = mode_a
        summary["mode_b"] = mode_b
        summary["paired_success_fraction_a_minus_b"] = (
            summary.a_only_zero_error - summary.b_only_zero_error) / summary.events
        summary["ber_a_common_decoded"] = summary.a_errors_common / summary.common_decoded_bits.replace(0, np.nan)
        summary["ber_b_common_decoded"] = summary.b_errors_common / summary.common_decoded_bits.replace(0, np.nan)
        summary["ber_difference_common_decoded"] = summary.ber_a_common_decoded - summary.ber_b_common_decoded
        output.append(summary)
    return pd.concat(output, ignore_index=True)


def continuous_per_record(data: pd.DataFrame) -> pd.DataFrame:
    """The horizon occurs ONCE per policy/channel/record, never per packet."""
    validate_counts(data)
    assert not data.duplicated(["policy", "channel", "record", "packet"]).any()
    keys = ["policy", "channel", "record"]
    assert (data.groupby(keys).horizon_s.nunique() == 1).all()
    assert np.allclose(data.horizon_s, 32)
    assert (data.control_s >= 0).all()
    assert np.allclose(data.next_start_s, data.start_s + data.occupied_s + data.control_s)
    assert (data.next_start_s <= data.horizon_s + 1e-9).all()
    for _, group in data.groupby(keys):
        group = group.sort_values("packet")
        assert np.array_equal(group.packet.to_numpy(), np.arange(1, len(group) + 1))
        assert (group.start_s.to_numpy()[1:] >= group.next_start_s.to_numpy()[:-1] - 1e-9).all()
    rates = aggregate_packets(data, keys)
    times = data.groupby(keys).agg(
        horizon_s=("horizon_s", "first"), occupied_seconds=("occupied_s", "sum"),
        control_seconds=("control_s", "sum"), first_tx_s=("start_s", "min"),
        final_budget_end_s=("next_start_s", "max"),
    ).reset_index()
    rates = rates.merge(times, on=keys, validate="one_to_one")
    rates["full_horizon_bps"] = rates.zero_error_payload_bits / rates.horizon_s
    rates["control_fraction_of_horizon"] = rates.control_seconds / rates.horizon_s
    rates["idle_rounding_unused_seconds"] = rates.horizon_s - rates.occupied_seconds - rates.control_seconds
    assert (rates.idle_rounding_unused_seconds >= -1e-8).all()
    return rates


def continuous_summary(per_record: pd.DataFrame) -> pd.DataFrame:
    sums = ["packets", "captured_packets", "payload_bits", "bit_errors", "bit_erasures", "packet_errors",
            "zero_error_packets", "zero_error_payload_bits", "decoded_bits", "horizon_s",
            "occupied_seconds", "control_seconds", "idle_rounding_unused_seconds"]
    result = per_record.groupby(["policy", "channel"], sort=True)[sums].sum().reset_index()
    counts = per_record.groupby(["policy", "channel"]).size().rename("records").reset_index()
    result = result.merge(counts, on=["policy", "channel"], validate="one_to_one")
    result["conditional_ber"] = result.bit_errors / result.decoded_bits.replace(0, np.nan)
    result["capture_loss"] = 1 - result.captured_packets / result.packets
    result["per"] = result.packet_errors / result.packets
    result["full_horizon_bps"] = result.zero_error_payload_bits / result.horizon_s
    result["control_fraction_of_horizon"] = result.control_seconds / result.horizon_s
    return result


def paired_continuous(per_record: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = ["channel", "records", "fixed_mean_bps", "adaptive_mean_bps", "mean_difference_bps",
               "resample_p025_bps", "resample_p975_bps", "resamples", "resample_seed"]
    if "adaptive_stump" not in set(per_record.policy):
        return pd.DataFrame(), pd.DataFrame(columns=columns)
    wide = per_record.pivot(index=["channel", "record"], columns="policy", values="full_horizon_bps")
    assert {"best_fixed", "adaptive_stump"} <= set(wide.columns)
    assert not wide[["best_fixed", "adaptive_stump"]].isna().any().any()
    wide["difference_bps"] = wide.adaptive_stump - wide.best_fixed
    paired = wide.reset_index()
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    rows = []
    channel_draw_means = []
    channel_sizes = []
    for channel, group in paired.groupby("channel", sort=True):
        delta = group.difference_bps.to_numpy()
        draws = rng.integers(0, len(delta), size=(BOOTSTRAP_REPEATS, len(delta)))
        means = delta[draws].mean(axis=1)
        channel_draw_means.append(means)
        channel_sizes.append(len(delta))
        lo, hi = np.quantile(means, [0.025, 0.975])
        rows.append(dict(channel=channel, records=len(group), fixed_mean_bps=group.best_fixed.mean(),
                         adaptive_mean_bps=group.adaptive_stump.mean(), mean_difference_bps=delta.mean(),
                         resample_p025_bps=lo, resample_p975_bps=hi,
                         resamples=BOOTSTRAP_REPEATS, resample_seed=BOOTSTRAP_SEED))
    combined = np.average(np.stack(channel_draw_means), axis=0, weights=channel_sizes)
    lo, hi = np.quantile(combined, [0.025, 0.975])
    rows.append(dict(channel="ALL (stratified)", records=len(paired), fixed_mean_bps=paired.best_fixed.mean(),
                     adaptive_mean_bps=paired.adaptive_stump.mean(), mean_difference_bps=paired.difference_bps.mean(),
                     resample_p025_bps=lo, resample_p975_bps=hi,
                     resamples=BOOTSTRAP_REPEATS, resample_seed=BOOTSTRAP_SEED))
    return paired, pd.DataFrame(rows, columns=columns)


def diagnostic_summary(test: pd.DataFrame) -> pd.DataFrame:
    data = test[(test.receiver_mode == "sync_compensated") & (test.detected == 1)].copy()
    data["bulk_velocity_deviation_mps"] = data.velocity_hat - data.velocity_true
    rows = []
    for (channel, snr), group in data.groupby(["channel", "snr_db"]):
        velocity = group.bulk_velocity_deviation_mps.dropna()
        timing = group.packet_start_error_ms_vs_mean_pdp.dropna()
        rows.append(dict(channel=channel, sample_snr_db=snr, decoded_packets=len(group),
                         bulk_velocity_deviation_median_mps=velocity.median(),
                         bulk_velocity_deviation_p10_mps=velocity.quantile(.1),
                         bulk_velocity_deviation_p90_mps=velocity.quantile(.9),
                         mean_pdp_peak_offset_median_ms=timing.median(),
                         mean_pdp_peak_offset_p10_ms=timing.quantile(.1),
                         mean_pdp_peak_offset_p90_ms=timing.quantile(.9)))
    return pd.DataFrame(rows)


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titleweight": "semibold",
                         "figure.facecolor": "white", "axes.facecolor": "white"})


def save_figure(fig, directory: Path, stem: str) -> Path:
    png = directory / f"{stem}.png"
    fig.savefig(png, dpi=170, bbox_inches="tight")
    fig.savefig(directory / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)
    return png.resolve()


def plot_receivers(summary: pd.DataFrame, directory: Path) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    metrics = [("conditional_ber", "Conditional payload BER"), ("per", "Packet error fraction"),
               ("slot_zero_error_bps", "Zero-error payload / full slot (bit/s)")]
    for row, channel in enumerate(["NOF1", "NCS1"]):
        for col, (metric, label) in enumerate(metrics):
            ax = axes[row, col]
            for mode in MODES:
                data = summary[(summary.channel == channel) & (summary.snr_db == 20)
                               & (summary.receiver_mode == mode)].sort_values("guard_ms")
                ax.plot(data.guard_ms, data[metric], "o-", ms=4, lw=1.6,
                        color=MODE_COLORS[mode], label=MODE_LABELS[mode])
            ax.set_title(channel)
            ax.set_xlabel("Guard duration (ms)")
            ax.set_ylabel(label)
            ax.set_xticks([0, 4, 8, 16, 32])
            ax.grid(alpha=.2)
            if metric == "conditional_ber":
                ax.set_yscale("symlog", linthresh=1e-4)
            elif metric == "per":
                ax.set_ylim(-.03, 1.03)
            else:
                ax.set_ylim(-.5, 20.5)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(.5, -.045))
    fig.suptitle("Held-out records 031-040 | sample-power SNR 20 dB | 128-bit payload\n"
                 "Same events; correction modes share acquisition; oracle is diagnostic", fontsize=12)
    return save_figure(fig, directory, "01_receiver_comparison_20db")


def plot_full_grid(summary: pd.DataFrame, directory: Path) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7), constrained_layout=True)
    metrics = [("capture_loss", "No complete captured payload", 1, "magma_r"),
               ("per", "Packet error fraction", 1, "magma_r"),
               ("slot_zero_error_bps", "Full-slot zero-error rate (bit/s)", 20, "viridis")]
    for row, channel in enumerate(["NOF1", "NCS1"]):
        for col, (metric, label, maximum, cmap) in enumerate(metrics):
            ax = axes[row, col]
            data = summary[(summary.channel == channel) & (summary.receiver_mode == "sync_compensated")]
            values = data.pivot(index="guard_ms", columns="snr_db", values=metric).reindex(index=[0, 4, 8, 16, 32], columns=[0, 10, 20])
            im = ax.imshow(values.to_numpy(), vmin=0, vmax=maximum, cmap=cmap, aspect="auto")
            for i in range(5):
                for j in range(3):
                    value = values.iloc[i, j]
                    text = "NA" if pd.isna(value) else (f"{value:.2f}" if maximum == 1 else f"{value:.1f}")
                    rgba = plt.get_cmap(cmap)(0 if pd.isna(value) else value/maximum)
                    luminance = .2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
                    ax.text(j, i, text, ha="center", va="center", color="white" if luminance < .5 else "black", fontsize=9)
            ax.set_xticks(range(3), [0, 10, 20])
            ax.set_yticks(range(5), [0, 4, 8, 16, 32])
            ax.set_xlabel("Reference sample-power SNR (dB)")
            ax.set_ylabel("Guard duration (ms)")
            ax.set_title(f"{channel} | {label}", fontsize=10)
            fig.colorbar(im, ax=ax, shrink=.8)
    fig.suptitle("Held-out records 031-040 | acquired + scale-corrected receiver\n"
                 "Each cell: 10 records x 5 matched 6.4 s listening slots", fontsize=12)
    return save_figure(fig, directory, "02_full_test_grid")


def plot_continuous(per_record: pd.DataFrame, differences: pd.DataFrame, directory: Path) -> Path:
    policies = [p for p in ["best_fixed", "adaptive_stump"] if p in set(per_record.policy)]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.7), constrained_layout=True, sharey=True)
    colors = {"best_fixed": "#60758d", "adaptive_stump": "#16798b"}
    for ax, channel in zip(axes, ["NOF1", "NCS1"]):
        channel_data = per_record[per_record.channel == channel]
        if len(policies) == 2:
            pivot = channel_data.pivot(index="record", columns="policy", values="full_horizon_bps")
            for _, row in pivot.iterrows():
                ax.plot([0, 1], [row[policies[0]], row[policies[1]]], "-", color="#aeb7bd", alpha=.65, lw=.9)
        for index, policy in enumerate(policies):
            values = channel_data[channel_data.policy == policy].full_horizon_bps.to_numpy()
            offsets = np.linspace(-.045, .045, len(values))
            ax.scatter(index+offsets, values, color=colors[policy], alpha=.8, s=28, zorder=3)
            ax.plot([index-.15, index+.15], [values.mean()]*2, color="#18252b", lw=3, zorder=4)
        ax.set_xticks(range(len(policies)), [POLICY_LABELS[p] for p in policies])
        ax.set_xlim(-.4, len(policies)-.6)
        ax.set_ylim(-.5, max(1, per_record.full_horizon_bps.max())*1.3+1)
        ax.set_ylabel("Zero-error payload / full 32 s (bit/s)")
        ax.set_title(channel)
        ax.grid(axis="y", alpha=.2)
        if len(policies) == 2:
            row = differences[differences.channel == channel].iloc[0]
            note = (f"Mean paired delta {row.mean_difference_bps:+.2f} bit/s\n"
                    f"2k record-resample range [{row.resample_p025_bps:+.2f}, {row.resample_p975_bps:+.2f}]")
        else:
            note = "Only the frozen fixed policy was simulated"
        ax.text(.5, .98, note, transform=ax.transAxes, ha="center", va="top", fontsize=9,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.85))
    fig.suptitle("Actual sequential scheduling | records 031-040 | sample-power SNR 20 dB\n"
                 "Dots = full records; lines = paired records; bars = means", fontsize=12)
    return save_figure(fig, directory, "03_continuous_record_rates")


def plot_schedule(data: pd.DataFrame, directory: Path) -> Path:
    from matplotlib.patches import Patch
    policies = [p for p in ["best_fixed", "adaptive_stump"] if p in set(data.policy)]
    fig, axes = plt.subplots(2, 1, figsize=(12, 4.8), constrained_layout=True, sharex=True)
    for ax, channel in zip(axes, ["NOF1", "NCS1"]):
        for index, policy in enumerate(policies):
            subset = data[(data.channel == channel) & (data.record == 31) & (data.policy == policy)]
            for row in subset.itertuples():
                color = "#258579" if row.packet_error == 0 else "#bd6953"
                ax.broken_barh([(row.start_s, row.occupied_s)], (index-.28, .56), facecolors=color)
                ax.broken_barh([(row.start_s+row.occupied_s, row.control_s)], (index-.28, .56), facecolors="#9ca6ac")
        ax.set_yticks(range(len(policies)), [POLICY_LABELS[p] for p in policies])
        ax.set_ylim(-.6, len(policies)-.4)
        ax.set_title(f"{channel} | record 031", loc="left")
        ax.set_xlim(0, 32)
        ax.grid(axis="x", alpha=.2)
    axes[-1].set_xlabel("Absolute transmit / measured-channel time (s)")
    fig.legend(handles=[Patch(color="#258579", label="Zero-error payload; frame + tail"),
                        Patch(color="#bd6953", label="Failed payload; frame + tail"),
                        Patch(color="#9ca6ac", label="Ideal feedback budget")],
               loc="lower center", bbox_to_anchor=(.5, -.06), ncol=3, fontsize=9)
    fig.suptitle("Different policies actually change subsequent packet times; unused time stays in the denominator", fontsize=11)
    return save_figure(fig, directory, "04_continuous_schedule_record031")


def md_table(data: pd.DataFrame) -> str:
    def fmt(value):
        if isinstance(value, (float, np.floating)):
            if not np.isfinite(value):
                return "NA"
            return f"{value:.5g}"
        return str(value)
    lines = ["|"+"|".join(data.columns)+"|", "|"+"|".join(["---"]*len(data.columns))+"|"]
    lines.extend("|"+"|".join(fmt(v) for v in row)+"|" for row in data.itertuples(index=False, name=None))
    return "\n".join(lines)


def file_link(label: str, path: Path) -> str:
    return f"[{label}](<{path.resolve().as_posix()}>)"


def image_link(label: str, path: Path) -> str:
    return f"![{label}](<{path.resolve().as_posix()}>)"


def generate(out: Path) -> dict:
    inputs = {name: out / f"metadata_{name}.json" for name in ["development", "test", "continuous"]}
    metadata = {name: require_complete(path) for name, path in inputs.items()}
    files = {name: out / f"records_{name}.csv" for name in ["development", "test", "continuous"]}
    data = {name: pd.read_csv(path) for name, path in files.items()}
    for stage in ["development", "test"]:
        validate_slot_data(data[stage], metadata[stage], stage)
    assert set(data["development"].record) == set(range(1, 11))
    assert set(data["continuous"].record) == set(range(31, 41))
    assert set(data["continuous"].channel) == {"NOF1", "NCS1"}
    assert set(data["continuous"].snr_db) == {20}
    assert np.allclose(data["continuous"].control_s, metadata["continuous"]["control_s"])
    validation_path = ROOT / "results/receiver_v2/validation.json"
    validation_csv = ROOT / "results/receiver_v2/validation_cases.csv"
    receiver_validation = read_json(validation_path)
    assert receiver_validation["status"] == "PASS"
    selftests = pd.read_csv(validation_csv)
    assert len(selftests) == receiver_validation["total_cases"]
    assert selftests.loc[selftests.required == 1, "passed"].astype(bool).all()
    tables_dir = out / "link_tables"
    plots_dir = out / "link_plots"
    tables_dir.mkdir(exist_ok=True)
    plots_dir.mkdir(exist_ok=True)

    summary_keys = ["channel", "guard_ms", "snr_db", "receiver_mode"]
    development = data["development"][data["development"].snr_db == 20]
    dev_summary = aggregate_packets(development, summary_keys, slot_denominator=True)
    summary = aggregate_packets(data["test"], summary_keys, slot_denominator=True)
    paired = paired_receivers(data["test"], ["channel", "guard_ms", "snr_db"])
    paired_balanced = paired_receivers(data["test"], ["channel", "snr_db"])
    per_record = continuous_per_record(data["continuous"])
    policies = set(per_record.policy)
    assert policies <= {"best_fixed", "adaptive_stump"} and "best_fixed" in policies
    expected = {(p, c, r) for p in policies for c in ["NOF1", "NCS1"] for r in range(31, 41)}
    assert set(per_record[["policy", "channel", "record"]].itertuples(index=False, name=None)) == expected
    continuous = continuous_summary(per_record)
    paired_records, differences = paired_continuous(per_record)
    diagnostics = diagnostic_summary(data["test"])
    overall = continuous.groupby("policy")[["zero_error_payload_bits", "horizon_s"]].sum()
    overall["full_horizon_bps"] = overall.zero_error_payload_bits/overall.horizon_s
    generated_tables = {
        "development_20db.csv": dev_summary, "heldout_slot_summary.csv": summary,
        "heldout_paired_receivers.csv": paired, "heldout_paired_receivers_balanced_guards.csv": paired_balanced,
        "continuous_per_record.csv": per_record, "continuous_summary.csv": continuous,
        "continuous_paired_records.csv": paired_records, "continuous_record_resampling.csv": differences,
        "heldout_reference_deviations.csv": diagnostics,
    }
    for name, table in generated_tables.items():
        table.to_csv(tables_dir / name, index=False)
    style()
    figures = [plot_receivers(summary, plots_dir), plot_full_grid(summary, plots_dir),
               plot_continuous(per_record, differences, plots_dir), plot_schedule(data["continuous"], plots_dir)]

    held20 = summary[(summary.snr_db == 20) & (summary.receiver_mode == "sync_compensated")]
    held_shown = held20[["channel", "guard_ms", "packets", "conditional_ber", "capture_loss", "per", "slot_zero_error_bps"]]
    continuous_shown = continuous[["channel", "policy", "records", "packets", "zero_error_packets", "conditional_ber",
                                  "capture_loss", "per", "full_horizon_bps", "control_seconds", "control_fraction_of_horizon"]]
    paired_shown = paired_balanced[(paired_balanced.mode_a == MODES[0]) & (paired_balanced.mode_b == MODES[1])][
        ["channel", "snr_db", "events", "a_only_zero_error", "b_only_zero_error", "both_zero_error", "neither_zero_error",
         "paired_success_fraction_a_minus_b", "ber_difference_common_decoded"]]
    oracle_shown = paired_balanced[(paired_balanced.mode_b == MODES[2]) & (paired_balanced.snr_db == 20)][
        ["channel", "mode_a", "events", "a_only_zero_error", "b_only_zero_error", "paired_success_fraction_a_minus_b"]]
    has_adaptive = "adaptive_stump" in policies
    if has_adaptive:
        fixed_rate = overall.loc["best_fixed", "full_horizon_bps"]
        adaptive_rate = overall.loc["adaptive_stump", "full_horizon_bps"]
        outcome = (f"完整 32 秒连续实验的跨信道合计速率为：冻结固定策略 **{fixed_rate:.3g} bit/s**，"
                   f"冻结自适应策略 **{adaptive_rate:.3g} bit/s**。")
        if adaptive_rate <= fixed_rate:
            outcome += "这次运行没有观察到自适应的吞吐优势；不能仅凭条件 BER 或 PER 改善宣称吞吐改进。"
        no_delivery = [channel for channel, group in continuous.groupby("channel")
                       if group.zero_error_packets.sum() == 0]
        if no_delivery:
            outcome += ("这两种冻结策略在 " + "/".join(no_delivery)
                        + " 的连续运行均未得到完全无错的整包；这与前导是否被捕获是不同指标。")
        adaptive_text = ("已在相同留出记录上运行冻结的固定策略和自适应策略；每种策略真的改变后续发送时刻。"
                         "下表以每条完整 32 秒记录为配对单位，自适应减固定；负值表示本次自适应更差。\n\n"
                         + md_table(differences))
    else:
        outcome = (f"完整 32 秒连续实验中，冻结固定策略的跨信道合计速率为 "
                   f"**{overall.loc['best_fixed', 'full_horizon_bps']:.3g} bit/s**；没有自适应连续结果可比较。")
        adaptive_text = ("本次连续结果只有冻结的固定策略，没有可报告的自适应连续运行。"
                         "因此不生成自适应相对固定的收益、配对差或重抽样范围。")
    fs = metadata["test"]["fs"]
    on_s = metadata["test"]["config"]["symbol_s"]
    db_gap = 10*np.log10(round(fs*on_s)/2)
    control_s = metadata["continuous"]["control_s"]
    noise_upper = receiver_validation["noise_only_zero_count_one_sided95_upper"]
    final_count = len(data["test"])
    event_count = len(data["test"].drop_duplicates(EVENT))
    baseline_noise = metadata["test"].get("noise", "See metadata_test.json")
    long_path = ROOT / "results/receiver_v2/long_packet_doppler.json"
    long_section = ""
    if long_path.exists():
        long_validation = read_json(long_path)
        assert long_validation["status"] == "PASS"
        long_table = pd.DataFrame(long_validation["results"])
        assert (long_table.payload_bits == 512).all()
        assert (long_table.errors_corrected == 0).all() and long_table.detected.all()
        long_section = ("另有一组**独立合成单径长包压力测试**，不参与策略拟合或真实信道统计："
                        "512 bit、未知相位与起点、20 dB 有效发射期样本 SNR，分别使用 −1.5/0/+1.5 m/s。"
                        "它演示较长数据块累积漂移时补偿器的能力，不代表实测水声多径中的性能。\n\n"
                        + md_table(long_table[["velocity_true", "velocity_hat", "payload_bits", "errors_corrected", "errors_uncorrected"]]))
    report = f"""# V2 接收链路与连续时序实验

> 2026-11-02 更正：本文件的连续实验部分保留原始运行。驱动的 `pi` 循环变量遮蔽数学常量，造成两个策略的初相分布不同，原连续策略比较不再作为正式结论。请使用 [修正复跑报告](../v3/continuous_phase_corrected/correction_report.md)。固定槽位结果及计数统计方法不受此问题影响。

正式结果已完成：固定槽实验使用 NOF1/NCS1 的留出记录 **031–040**，共 {event_count} 个匹配事件、{final_count} 行接收模式结果；连续实验也使用这 20 条记录，共 **{len(data['continuous'])} 个发送帧**，每个策略每条记录的计时分母固定为 **完整 32 秒**。开发记录 001–010 仅以 20 dB 结果作背景，未与留出结果混合。smoke 未进入本报告。

{outcome}

## 统计口径与范围

- 条件 BER = 实际比特错误数 /（发送载荷 bit 数 − 擦除 bit 数）。未得到完整载荷的包保留为捕获损失、128 bit 擦除和包错误；不会把它算作零误码。若完全没有判决，条件 BER 为 NA。
- `capture_loss` 是未输出完整载荷的比例，并非通过真值确认的独立定时关联失败率。PER 包含捕获失败和至少一个载荷错误的包。
- 无错载荷速率只计整包完全正确的载荷 bit。固定槽分母包含所有 6.4 s 监听槽；连续分母每个策略/信道/记录只计一次 32 s，含启动、反馈、取整和末尾空闲。这里不使用“占用时间效率”代替真实排程速率。
- 尚无 CRC/ARQ 协议。整包是否正确由离线载荷真值统计，故速率是**理想无错载荷统计值**，不是已实现的接收确认吞吐量。策略和接收器不读取这份正确性真值。

零误码观测也不意味着真实 BER 为零。这里保留计数与分母，未把相关的比特或包机械当成独立海况样本。

本轮横轴是固定参考接收功率的 **样本 SNR**：`variance = calibrationPower / 10^(SNR/10)`。8 ms、48 kHz 的 BFSK 参考换算与载荷有效能量的 `Eb/N0` 相差约 **{db_gap:.2f} dB**。因此本报告的 0/10/20 dB 不能直接当作 AWGN 自检或官方 `sfetch` 的同值 `Eb/N0`，也不能据此宣称精确复现原文同名 SNR 条件。每条记录单独校准一次，并在保护时间和策略之间固定噪声功率。元数据原说明：`{baseline_noise}`。

## 留出记录上的接收模式对照

下表是 20 dB 的完整同步及尺度补偿接收机，每行 10 条记录 × 5 个槽。所有保护时间均列出，没有按留出结果挑出一个最接近原论文的配置。

{md_table(held_shown)}

{image_link("20 dB 接收模式对照", figures[0])}

补偿与不补偿使用**同一次前导捕获和同一个估计起点**，只改变载荷的尺度/频率及时序校正。下面按信道与 SNR 汇总共同事件，五个预设保护时间等权出现；它是实验条件对照，不是一条实际策略轨迹。`a_only_zero_error` 表示仅补偿包完全正确，`b_only_zero_error` 表示仅不补偿包完全正确；BER 差只在双方都有完整判决的共同事件上计算。

{md_table(paired_shown)}

20 dB oracle 诊断使用已知发射起点、全记录平均功率时延轮廓的峰位置及恢复的 bulk 尺度。它既不是可部署接收机，也不保证是性能上界；该参考峰未必对应每包最佳判决时刻。以下 `mode_a` 与 oracle 在相同事件上比较，不能把差值全部归因于捕获算法。

{md_table(oracle_shown)}

{image_link("全部 SNR 与保护时间测试", figures[1])}

所有 0/10/20 dB、五种保护时间和三种接收模式的分母与计数保存在 {file_link("完整留出统计表", tables_dir / "heldout_slot_summary.csv")}；{file_link("逐保护时间的共同事件对照", tables_dir / "heldout_paired_receivers.csv")} 保留了模式间配对关系。

## 真正连续排程的结果

每发一个包，下一次发射时刻前进“所选完整帧长 + 128 ms 公共尾段 + {control_s:.3f} s 控制预算”，并在新的绝对 TVIR 时间回放信道。随机接收窗口前置背景只用于未知起点，不作为发送空闲或重复加入时间分母。不同策略共享按绝对观察时间索引的噪声序列；包数和发送时刻可以不同。

{md_table(continuous_shown)}

`control_seconds` 是该信道全部记录实际安排的控制预算总和，`control_fraction_of_horizon` 除以相应完整记录时间总和；它不代表已测得的控制信道耗时。连续结果只覆盖 20 dB，不能外推到 0/10 dB。

{adaptive_text}

重抽样以配对记录为单位，固定随机种子 {BOOTSTRAP_SEED}，共 {BOOTSTRAP_REPEATS} 次，显示均值差的 2.5%–97.5% 分位范围。总体行先在每个信道内重抽记录，再按原记录数合并，保持两信道的原有权重。同一记录内的包相关；不同记录也可能来自相关的连续测量，因此这里是对当前记录集合的**描述性重抽样范围**，不是独立海况的总体置信保证，也不据此单独宣称显著优越。

{image_link("连续记录的实际无错载荷速率", figures[2])}

{image_link("连续发送轨迹示例", figures[3])}

{file_link("每条完整记录的计数和速率", tables_dir / "continuous_per_record.csv")} 包含每条记录唯一的 32 s 分母；{file_link("配对记录差值", tables_dir / "continuous_paired_records.csv")} 与 {file_link("两千次记录重抽样", tables_dir / "continuous_record_resampling.csv")} 可复核上面的比较。固定槽的占用时间推算和连续实际速率未混合。

## 接收机验证与解释边界

接收机自检 `{receiver_validation['status']}`：必需项 {receiver_validation['required_passed']}/{receiver_validation['required_cases']} 通过；纯噪声误检 {receiver_validation['noise_only_false_alarms']}/{receiver_validation['noise_only_cases']}。零次观察的单侧 95% 上界为 {noise_upper:.2%}，不能把配置中的 0.1% 设计目标说成实测误检率。另有 {receiver_validation['ebn0_zero_stress_cases']} 个 `Eb/N0=0 dB` 压力样本，仅捕获 {receiver_validation['ebn0_zero_stress_detected']} 个；它们未被隐藏，也未作为要求零误码的自检。该 dB 定义与本报告样本 SNR 不同。

{long_section}

HFM 两段各 64 ms、中间 64 ms 间隔；10–18 kHz 扫频、加窗、两个额外 64 ms 已知单音导频及保护段是公开声明的本项目选择。载荷标称偏移 528 ms。接收器只收到波形、公共配置、载荷长度和补偿开关，没有真实信道、V0、发射起点或载荷。

连续信道使用绝对时间 TVIR 回放，并逐包恢复 bulk 时间尺度；它没有完整模拟整段链路累计传播时延及往返水声反馈。反馈被假设在 {control_s:.3f} s 预算内理想成功、延迟一个包生效；数据捕获失败触发回退。它不包含真实返回链路的可靠性、传播耗时或控制误码，也没有将离线 MATLAB 运行耗时当作实际设备处理时延。

`velocity_hat − V0` 只称为**与恢复的 bulk 速度之差**：实测 TVIR 的残余多普勒也进入导频估计，不能据此直接称为纯速度估计误差。`packet_start_error_ms_vs_mean_pdp` 只称为**相对平均 PDP 参考峰的偏移**，不是无歧义的真实首径起点误差。按已捕获包汇总如下，五种保护时间均包含在内。

{md_table(diagnostics)}

本报告的实验范围、噪声定义和反馈假设来自各阶段完成后的元数据；输入文件 SHA-256、计数与分母检查记录在 {file_link("报告核验记录", out / "link_report_validation.json")}。开发背景统计单独保存在 {file_link("开发集 20 dB 统计", tables_dir / "development_20db.csv")}。
"""
    (out / "link_report.md").write_text(report, encoding="utf-8")
    source_files = list(inputs.values()) + list(files.values()) + [validation_path, validation_csv, Path(__file__)]
    if long_path.exists():
        source_files.append(long_path)
    receipt = {
        "status": "PASS", "formal": True, "test_records": list(range(31, 41)),
        "test_event_count": event_count, "test_rows": final_count,
        "continuous_packets": len(data["continuous"]), "continuous_record_policy_units": len(per_record),
        "continuous_policy_names": sorted(policies), "continuous_horizon_seconds_per_unit": 32,
        "denominator_once_per_record_verified": True, "paired_receiver_events_verified": True,
        "resampling": {"repeats": BOOTSTRAP_REPEATS, "seed": BOOTSTRAP_SEED, "unit": "paired record within channel; overall stratified by channel"},
        "figures": [p.as_posix() for p in figures],
        "sha256": {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files},
    }
    (out / "link_report_validation.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8")
    return receipt


def self_check(out: Path) -> dict:
    smoke_metadata = require_complete(out / "metadata_smoke.json")
    smoke = pd.read_csv(out / "records_smoke.csv")
    validate_slot_data(smoke, smoke_metadata, "smoke")
    summary = aggregate_packets(smoke, ["channel", "guard_ms", "snr_db", "receiver_mode"], True)
    assert (summary.packets == smoke_metadata["slots_per_record"]).all()
    assert np.allclose(summary.listening_seconds, summary.packets*smoke_metadata["slot_s"])
    paired = paired_receivers(smoke, ["channel", "snr_db"])
    assert (paired.events == len(as_list(smoke_metadata["guards_ms"]))*smoke_metadata["slots_per_record"]).all()
    assert (paired.a_only_zero_error+paired.b_only_zero_error+paired.both_zero_error+paired.neither_zero_error == paired.events).all()
    synthetic = []
    for channel in ["NOF1", "NCS1"]:
        for record in [31, 32]:
            for policy, count in [("best_fixed", 2), ("adaptive_stump", 3)]:
                for packet in range(1, count+1):
                    start = .05+(packet-1)*4.16
                    synthetic.append(dict(channel=channel, record=record, policy=policy, packet=packet,
                                          payload_bits=128, detected=1, bit_errors=0, bit_erasures=0, packet_error=0,
                                          horizon_s=32, occupied_s=4, control_s=.16, start_s=start, next_start_s=start+4.16))
    per_record = continuous_per_record(pd.DataFrame(synthetic))
    assert np.allclose(per_record[per_record.policy == "best_fixed"].full_horizon_bps, 8)
    assert np.allclose(per_record[per_record.policy == "adaptive_stump"].full_horizon_bps, 12)
    records, difference = paired_continuous(per_record)
    assert np.allclose(records.difference_bps, 4)
    assert np.allclose(difference[["mean_difference_bps", "resample_p025_bps", "resample_p975_bps"]], 4)
    fixed_only = per_record[per_record.policy == "best_fixed"]
    assert paired_continuous(fixed_only)[0].empty
    return {"status": "PASS", "scope": "statistics self-check only; no formal report generated",
            "smoke_rows_checked": len(smoke), "unique_horizon_denominator_checked": True,
            "matched_event_counts_checked": True, "paired_resampling_checked": True,
            "fixed_only_case_checked": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT / "results/v2")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    out = args.results.resolve()
    result = self_check(out) if args.self_check else generate(out)
    print(json.dumps({k: v for k, v in result.items() if k != "sha256"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
