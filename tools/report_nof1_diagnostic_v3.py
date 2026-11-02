"""Paired, development-only NOF1 timing/phase/detector diagnosis.

Refuses incomplete results, baseline regressions, or broken common-random
conditions. No held-out records or invalid first-run artifacts are consumed.
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
CONDITION = ["design", "on_ms", "guard_ms", "f0_hz", "f1_hz", "phase_mode", "noise_mode"]
EVENT = ["channel", "record", "slot"] + CONDITION
METHODS = ["first_est", "strong_est", "first_true_scale", "strong_true_scale",
           "pdp_true", "pdp_est_scale", "first_est_fft_nearest"]
LABELS = {
    "first_est": "First-arrival\nest. scale",
    "strong_est": "HFM main peak\nest. scale",
    "first_true_scale": "First-arrival\ntrue bulk scale",
    "strong_true_scale": "HFM main peak\ntrue bulk scale",
    "pdp_true": "Mean-PDP peak\ntrue bulk scale",
    "pdp_est_scale": "Mean-PDP peak\nest. scale",
    "first_est_fft_nearest": "First-arrival\nnearest FFT bin",
}
COLORS = ["#bd6748", "#187c8c", "#8395a3", "#599b9a", "#6b62a0", "#a7a0c0", "#b09249"]


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def aggregate(data: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    specs = dict(packets=("payload_bits", "size"), bits=("payload_bits", "sum"),
                 errors=("bit_errors", "sum"))
    if "packet_error" in data:
        specs["packet_errors"] = ("packet_error", "sum")
    result = data.groupby(keys, sort=True, dropna=False).agg(**specs).reset_index()
    result["ber"] = result.errors/result.bits
    if "packet_errors" in result:
        result["per"] = result.packet_errors/result.packets
    return result


def paired_against_first(decoders: pd.DataFrame) -> pd.DataFrame:
    wide = decoders.pivot(index=EVENT, columns="method", values="bit_errors")
    assert not wide.isna().any().any()
    output = []
    for method in METHODS:
        comparison = wide.index.to_frame(index=False)
        baseline = wide.first_est.to_numpy()
        current = wide[method].to_numpy()
        comparison["events"] = 1
        comparison["baseline_errors"] = baseline
        comparison["method_errors"] = current
        comparison["fewer_error_packets"] = (current < baseline).astype(int)
        comparison["more_error_packets"] = (current > baseline).astype(int)
        comparison["equal_error_packets"] = (current == baseline).astype(int)
        comparison["zero_error_packet_delta"] = (current == 0).astype(int)-(baseline == 0).astype(int)
        numbers = [column for column in comparison if column not in EVENT]
        grouped = comparison.groupby(CONDITION)[numbers].sum().reset_index()
        grouped["method"] = method
        grouped["delta_errors_vs_first"] = grouped.method_errors-grouped.baseline_errors
        grouped["delta_ber_vs_first"] = grouped.delta_errors_vs_first/(grouped.events*64)
        output.append(grouped)
    return pd.concat(output, ignore_index=True)


def validate(metadata, decoders, events, sweep, old):
    assert metadata.get("status") == "COMPLETE", "Refuse incomplete diagnostic output"
    assert metadata.get("stage") == "full", "Formal report requires the full 001-003 run"
    assert set(metadata["record_ids"]) == {1, 2, 3}
    for data in [decoders, events, sweep]:
        assert set(data.channel) == {"NOF1"}
        assert set(data.record) == {1, 2, 3}, "Only development records 001-003 are allowed"
        assert set(data.slot) == {1, 2, 3, 4, 5}
    assert len(events) == metadata["event_rows"] == 390
    assert len(decoders) == metadata["decoder_rows"] == 2730
    assert len(sweep) == metadata["timing_sweep_rows"] == 990
    assert not events.duplicated(EVENT).any()
    assert not decoders.duplicated(EVENT+["method"]).any()
    assert not sweep.duplicated(EVENT+["shift_ms"]).any()
    assert set(decoders.method) == set(METHODS)
    assert (decoders.groupby(EVENT).method.nunique() == 7).all()
    event_tuples = set(events[EVENT].itertuples(index=False, name=None))
    assert set(decoders[EVENT].drop_duplicates().itertuples(index=False, name=None)) == event_tuples
    assert (events.detected == 1).all() and (events.original_decoder_match == 1).all()
    assert (decoders.payload_bits == 64).all() and (sweep.payload_bits == 64).all()
    assert decoders.bit_errors.between(0, 64).all() and sweep.bit_errors.between(0, 64).all()
    assert (decoders.packet_error == (decoders.bit_errors > 0).astype(int)).all()
    assert set(events.noise_mode) == {"clean", "paired_sample20"}
    assert (events[events.noise_mode == "clean"].noise_variance == 0).all()

    # Check the whole prespecified factorial grid, including the extra zero
    # payload-phase diagnostic only in design 1.
    design = np.asarray(metadata["design"], dtype=float)
    expected = set()
    for record in [1, 2, 3]:
        for slot in range(1, 6):
            for index, (on, guard, f0, f1) in enumerate(design, start=1):
                phases = ["reset", "continuous"] + (["payload_reset_zero"] if index == 1 else [])
                for phase in phases:
                    for noise in ["clean", "paired_sample20"]:
                        expected.add(("NOF1", record, slot, index, on, guard, f0, f1, phase, noise))
    assert event_tuples == expected, "Missing or extra diagnostic condition"

    # Prefix phase/start, standard-noise seeds and noise variance must not
    # change between waveform variants. The MATLAB driver fixes random seeds;
    # emitted phase/start/variance values provide an independent cross-check.
    for name in ["carrier_phase_rad", "injected_start_s"]:
        assert (events.groupby(["record", "slot"])[name].nunique() == 1).all(), name
    assert events.carrier_phase_rad.between(0, 2*np.pi).all()
    same_noise = events.groupby(["record", "slot", "on_ms", "guard_ms", "noise_mode"]).noise_variance
    assert (same_noise.nunique() == 1).all(), "Noise variance changed across phase/frequency variants"
    assert np.allclose(events[events.phase_mode != "payload_reset_zero"].payload_initial_phase_rad,
                       events[events.phase_mode != "payload_reset_zero"].carrier_phase_rad, rtol=0, atol=1e-12)
    assert (events[events.phase_mode == "payload_reset_zero"].payload_initial_phase_rad == 0).all()

    # Same-waveform ablations may replace only the declared timing/scale or
    # detector. In particular FFT and exact-tone decisions share sample blocks.
    wide_start = decoders.pivot(index=EVENT, columns="method", values="preamble_start_sample")
    wide_scale = decoders.pivot(index=EVENT, columns="method", values="alpha_used")
    for first, second in [("first_est", "first_true_scale"), ("first_est", "first_est_fft_nearest"),
                          ("strong_est", "strong_true_scale"), ("pdp_true", "pdp_est_scale")]:
        assert np.allclose(wide_start[first], wide_start[second], rtol=0, atol=1e-7)
    for first, second in [("first_est", "strong_est"), ("first_est", "pdp_est_scale"),
                          ("first_est", "first_est_fft_nearest"), ("first_true_scale", "strong_true_scale"),
                          ("first_true_scale", "pdp_true")]:
        assert np.allclose(wide_scale[first], wide_scale[second], rtol=0, atol=1e-12)
    assert (decoders[decoders.method == "first_est_fft_nearest"].detector == "fft_nearest").all()
    assert (decoders[decoders.method != "first_est_fft_nearest"].detector == "exact").all()

    # Original V2 regression: every packet in all four unchanged controls.
    reference = old[(old.channel == "NOF1") & (old.snr_label_db == 20)
                    & (old.noise_mode == "sample_snr_packet")]
    baseline = decoders[(decoders.design <= 4) & (decoders.phase_mode == "reset")
                        & (decoders.noise_mode == "paired_sample20") & (decoders.method == "first_est")]
    baseline = baseline.merge(events[EVENT+["noise_variance", "carrier_phase_rad", "injected_start_s"]],
                              on=EVENT, validate="one_to_one")
    regression = baseline.merge(reference, on=["channel", "record", "slot", "on_ms", "guard_ms"],
                                suffixes=("_new", "_old"), validate="one_to_one")
    assert len(regression) == 60
    for field in ["bit_errors", "payload_bits", "packet_error"]:
        assert np.array_equal(regression[field+"_new"], regression[field+"_old"]), field
    assert np.allclose(regression.noise_variance_new, regression.noise_variance_old, rtol=1e-9, atol=1e-10)
    assert np.allclose(regression.carrier_phase_rad_new, regression.carrier_phase_rad_old, rtol=0, atol=1e-12)
    assert np.allclose(regression.injected_start_s_new, regression.injected_start_s_old, rtol=0, atol=1e-12)
    primary = regression[(regression.on_ms == 4) & (regression.guard_ms == 0)]
    assert len(primary) == 15 and primary.bit_errors_new.sum() == 373 and primary.payload_bits_new.sum() == 960

    # Designs 2-6 have integer cycles over active+guard stride for both tones;
    # reset and continuous are equivalent waveforms apart from float roundoff.
    control_keys = ["channel", "record", "slot", "design", "on_ms", "guard_ms", "f0_hz", "f1_hz", "noise_mode"]
    controls = events[events.design >= 2]
    reset = controls[controls.phase_mode == "reset"]
    continuous = controls[controls.phase_mode == "continuous"]
    invariants = reset.merge(continuous, on=control_keys, suffixes=("_reset", "_continuous"), validate="one_to_one")
    assert len(invariants) == 150
    for field in ["noise_variance", "carrier_phase_rad", "injected_start_s"]:
        assert np.array_equal(invariants[field+"_reset"], invariants[field+"_continuous"]), field
    assert np.allclose(invariants.received_energy_reset, invariants.received_energy_continuous, rtol=1e-8, atol=1e-10)
    assert np.allclose(invariants.alpha_hat_reset, invariants.alpha_hat_continuous, rtol=0, atol=1e-9)
    control_bits = decoders[decoders.design >= 2].pivot(index=control_keys+["method"], columns="phase_mode", values="bit_errors")
    assert np.array_equal(control_bits.reset, control_bits.continuous), "Equivalent phase controls changed decisions"

    assert set(sweep.design) == {1, 3}
    assert set(sweep.phase_mode) == {"reset"} and set(sweep.noise_mode) == {"paired_sample20"}
    assert (sweep.guard_ms == 0).all()
    assert np.array_equal(np.sort(sweep.shift_ms.unique()), np.arange(-4, 4.001, .25))
    assert (sweep.groupby(EVENT).shift_ms.nunique() == 33).all()
    return regression, invariants


def waveform_label(row) -> str:
    return f"{int(row.on_ms)}/{int(row.guard_ms)} ms | {int(row.f0_hz)}/{int(row.f1_hz)} | {row.phase_mode}"


def save_figure(fig, out, name):
    path = out/f"{name}.png"
    fig.savefig(path, dpi=175, bbox_inches="tight")
    fig.savefig(out/f"{name}.svg", bbox_inches="tight")
    plt.close(fig)
    return path.resolve()


def plot_same_wave(summary, per_record, out):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), sharey=True, constrained_layout=True)
    for ax, noise in zip(axes, ["clean", "paired_sample20"]):
        subset = summary[(summary.design == 1) & (summary.phase_mode == "reset") & (summary.noise_mode == noise)]
        subset = subset.set_index("method").reindex(METHODS)
        ax.bar(np.arange(7), subset.ber, color=COLORS, width=.65)
        records = per_record[(per_record.design == 1) & (per_record.phase_mode == "reset") & (per_record.noise_mode == noise)]
        for index, method in enumerate(METHODS):
            values = records[records.method == method].sort_values("record").ber.to_numpy()
            ax.scatter(index+np.linspace(-.12, .12, len(values)), values, marker="o", s=22,
                       facecolor="white", edgecolor="#26373f", linewidth=.8, zorder=3)
            ax.text(index, subset.loc[method, "ber"]+.012,
                    f"{int(subset.loc[method, 'errors'])}/960", ha="center", va="bottom", fontsize=8)
        if noise == "paired_sample20":
            ax.axhline(.1137, color="#525d65", ls="--", lw=1, label="Paper 20 dB label: 0.1137 (different assumptions)")
            ax.legend(fontsize=7, loc="upper right")
        ax.set_xticks(range(7), [LABELS[m] for m in METHODS], rotation=38, ha="right", fontsize=8)
        ax.set_ylabel("Payload bit error rate")
        ax.set_ylim(0, .65)
        ax.set_title("No added noise" if noise == "clean" else "Shared baseline sample-SNR 20 dB noise")
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("NOF1 001-003 | 4 ms active, no guard | SAME received waveform\n"
                 "Bars: all 15 packets; dots: three record BERs; HFM main peak is not a physical-path oracle", fontsize=12)
    return save_figure(fig, out, "01_same_waveform_ablation")


def plot_waveform(summary, out):
    cases = [(1, "reset"), (1, "continuous"), (1, "payload_reset_zero"),
             (5, "reset"), (5, "continuous"), (6, "reset"), (6, "continuous")]
    labels = ["13875/14125\nreset", "13875/14125\ncontinuous", "13875/14125\npayload phase=0",
              "13750/14000\nreset", "13750/14000\ncontinuous", "14000/14250\nreset", "14000/14250\ncontinuous"]
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5), sharey=True, constrained_layout=True)
    for ax, noise in zip(axes, ["clean", "paired_sample20"]):
        subset = summary[summary.noise_mode == noise]
        for method, offset, color, label in [("first_est", -.18, COLORS[0], "First-arrival timing"),
                                             ("strong_est", .18, COLORS[1], "HFM main-peak timing")]:
            values = [subset[(subset.design == design) & (subset.phase_mode == phase) & (subset.method == method)].iloc[0].ber
                      for design, phase in cases]
            ax.bar(np.arange(len(cases))+offset, values, .34, color=color, label=label)
        ax.axvline(2.5, color="#a5afb4", lw=.8)
        ax.axvline(4.5, color="#a5afb4", lw=.8)
        ax.set_xticks(range(7), labels, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel("Payload bit error rate")
        ax.set_ylim(0, .65)
        ax.set_title("No added noise" if noise == "clean" else "Identical noise variance and realization")
        ax.grid(axis="y", alpha=.2)
    axes[0].legend(fontsize=9, loc="upper right")
    fig.suptitle("4 ms / no guard | Prespecified waveform diagnostics, not frequency or policy tuning\n"
                 "Changing tones also changes measured-channel fading and known-tone pilots", fontsize=12)
    return save_figure(fig, out, "02_phase_frequency_interaction")


def plot_sweep(sweep, out):
    combined = aggregate(sweep, ["on_ms", "shift_ms"])
    by_record = aggregate(sweep, ["on_ms", "record", "shift_ms"])
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), sharey=True, constrained_layout=True)
    for ax, duration in zip(axes, [4, 8]):
        for record, color in [(1, "#8bb0c0"), (2, "#b5a5be"), (3, "#b9ad88")]:
            data = by_record[(by_record.on_ms == duration) & (by_record.record == record)].sort_values("shift_ms")
            ax.plot(data.shift_ms, data.ber, color=color, lw=1.2, label=f"Record {record:03d}")
        data = combined[combined.on_ms == duration].sort_values("shift_ms")
        ax.plot(data.shift_ms, data.ber, color="#1f505c", lw=2.5, marker="o", ms=3, label="All 15 packets")
        ax.axvline(0, ls="--", color="#545d62", lw=1)
        ax.set_title(f"Active {duration} ms; no guard")
        ax.set_xlabel("Timing shift from full-record mean-PDP peak (ms)")
        ax.set_ylabel("Payload bit error rate")
        ax.set_xticks(np.arange(-4, 5, 1))
        ax.set_ylim(0, 1)
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Diagnostic sweep only | true bulk scale | baseline phase-reset waveform\n"
                 "All 33 offsets shown; no BER-minimizing offset is selected for deployment", fontsize=12)
    return save_figure(fig, out, "03_timing_sensitivity_sweep")


def md_table(data):
    def fmt(value):
        if isinstance(value, (float, np.floating)):
            return "NA" if not np.isfinite(value) else f"{value:.5g}"
        return str(value)
    return "\n".join(["|"+"|".join(data.columns)+"|", "|"+"|".join(["---"]*len(data.columns))+"|"]
                     + ["|"+"|".join(fmt(value) for value in row)+"|" for row in data.itertuples(index=False, name=None)])


def link(label, path):
    return f"[{label}](<{path.resolve().as_posix()}>)"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT/"results/v3/nof1_diagnostic")
    args = parser.parse_args()
    out = args.results.resolve()
    metadata_path = out/"metadata.json"
    metadata = read_json(metadata_path)
    if metadata.get("status") != "COMPLETE":
        raise ValueError("Diagnostic run has not completed; no formal report generated")
    paths = [out/"decoder_results.csv", out/"events.csv", out/"timing_sweep.csv"]
    decoders, events, sweep = [pd.read_csv(path) for path in paths]
    old_path = ROOT/"results/v2/paper_alignment/per_packet.csv"
    old = pd.read_csv(old_path)
    regression, invariants = validate(metadata, decoders, events, sweep, old)
    summary = aggregate(decoders, CONDITION+["method"])
    per_record = aggregate(decoders, CONDITION+["method", "record"])
    paired = paired_against_first(decoders)
    sweep_summary = aggregate(sweep, ["design", "on_ms", "shift_ms"])
    sweep_record = aggregate(sweep, ["design", "on_ms", "record", "shift_ms"])
    plots_dir = out/"plots"
    plots_dir.mkdir(exist_ok=True)
    summary.to_csv(out/"summary_all_conditions.csv", index=False)
    per_record.to_csv(out/"per_record_summary.csv", index=False)
    paired.to_csv(out/"paired_decoder_changes.csv", index=False)
    regression.to_csv(out/"baseline_regression.csv", index=False)
    invariants.to_csv(out/"integer_stride_phase_invariants.csv", index=False)
    sweep_summary.to_csv(out/"timing_sweep_summary.csv", index=False)
    sweep_record.to_csv(out/"timing_sweep_per_record.csv", index=False)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titleweight": "semibold", "figure.facecolor": "white"})
    figures = [plot_same_wave(summary, per_record, plots_dir), plot_waveform(summary, plots_dir), plot_sweep(sweep, plots_dir)]
    primary = summary[(summary.design == 1) & (summary.phase_mode == "reset") & (summary.noise_mode == "paired_sample20")].set_index("method")
    first = primary.loc["first_est"]
    strong = primary.loc["strong_est"]
    scale = primary.loc["first_true_scale"]
    fft = primary.loc["first_est_fft_nearest"]
    pdp = primary.loc["pdp_true"]
    timing = events[(events.design == 1) & (events.phase_mode == "reset") & (events.noise_mode == "paired_sample20")]
    timing_shown = timing.groupby("record").agg(
        packets=("slot", "size"), first_minus_pdp_median_ms=("first_minus_pdp_ms", "median"),
        mainpeak_minus_pdp_median_ms=("strong_minus_pdp_ms", "median"),
        mainpeak_minus_first_median_ms=("strong_minus_first_ms", "median"),
        alpha_hat_median=("alpha_hat", "median"), alpha_true_median=("alpha_true", "median"),
    ).reset_index()
    primary_table = summary[(summary.design == 1) & (summary.phase_mode == "reset")][
        ["noise_mode", "method", "packets", "bits", "errors", "ber", "packet_errors", "per"]]
    primary_records = per_record[(per_record.design == 1) & (per_record.phase_mode == "reset")
                                 & (per_record.noise_mode == "paired_sample20")].pivot(index="record", columns="method", values="ber").reindex(columns=METHODS).reset_index()
    phase_table = summary[(summary.design == 1) & (summary.noise_mode == "paired_sample20")
                          & summary.method.isin(["first_est", "strong_est", "pdp_true"])][
                              ["phase_mode", "method", "errors", "bits", "ber", "per"]]
    variant_tables = []
    for noise in ["clean", "paired_sample20"]:
        table = summary[summary.noise_mode == noise].copy()
        table["waveform"] = [waveform_label(row) for row in table.itertuples()]
        pivot = table.pivot(index="waveform", columns="method", values="ber").reindex(columns=METHODS).reset_index()
        variant_tables.append(f"**{noise}**\n\n"+md_table(pivot))
    primary_change = paired[(paired.design == 1) & (paired.phase_mode == "reset") & (paired.noise_mode == "paired_sample20")][
        ["method", "events", "delta_errors_vs_first", "delta_ber_vs_first", "fewer_error_packets", "more_error_packets", "equal_error_packets"]]
    clean_errors = int(summary[(summary.design == 1) & (summary.phase_mode == "reset") & (summary.noise_mode == "clean") & (summary.method == "first_est")].iloc[0].errors)
    clean_scale_errors = int(summary[(summary.design == 1) & (summary.phase_mode == "reset") & (summary.noise_mode == "clean") & (summary.method == "first_true_scale")].iloc[0].errors)
    continuous_phase = summary[(summary.design == 1) & (summary.phase_mode == "continuous") & (summary.noise_mode == "paired_sample20")].set_index("method")
    mechanism = (f"前导主峰平均比首到达估计晚 **{timing.strong_minus_first_ms.mean():.6f} ms**。"
                 f"无加噪声时原首到达窗口仍有 {clean_errors}/960 错误，同一无加噪声条件换真实bulk尺度后为 {clean_scale_errors}/960；20 dB主条件换真实bulk尺度后为 {int(scale.errors)}/960，"
                 "所以噪声和bulk尺度不是足以解释这一现象的主因。4 ms下两音均为半整数周期，逐符号同相reset对应符号边界的π反转；"
                 "过早的积分窗可能混入相邻符号成分，并在特定比特/路径相位下相消。"
                 f"连续相位在首到达窗口为 {int(continuous_phase.loc['first_est', 'errors'])}/960，"
                 f"前导主峰窗口为 {int(continuous_phase.loc['strong_est', 'errors'])}/960；"
                 f"后者并未优于reset的 {int(strong.errors)}/960，因此不能把连续相位说成普遍更好。"
                 "这支持“相位边界×提前积分窗口”的交互解释，而不是单一调制方式排名。")
    toy_path = out/"toy/README.md"
    if toy_path.exists():
        mechanism += ("独立的 " + link("双径toy与受限相消公式", toy_path)
                      + " 只提供机制佐证；公式要求同一比特、单路径、矩形窗等条件，不能直接预测真实NOF1的BER。")
    confirmation_section = ""
    confirmation_sources = []
    confirmation_path = out/"confirmation"
    if (confirmation_path/"metadata.json").exists():
        cm = read_json(confirmation_path/"metadata.json")
        assert cm["status"] == "COMPLETE" and set(cm["record_ids"]) == {4, 5}
        cd = pd.read_csv(confirmation_path/"decoder_results.csv")
        ce = pd.read_csv(confirmation_path/"events.csv")
        assert len(cd) == 280 and len(ce) == 40
        assert set(cd.record) == {4, 5} and set(cd.design) == {1, 3}
        assert set(cd.phase_mode) == {"reset"} and (cd.guard_ms == 0).all()
        assert not cd.duplicated(EVENT+["method"]).any() and not ce.duplicated(EVENT).any()
        assert (cd.groupby(EVENT).method.nunique() == 7).all()
        assert (ce.detected == 1).all() and (ce.original_decoder_match == 1).all()
        assert (cd.payload_bits == 64).all() and cd.bit_errors.between(0, 64).all()
        cs = aggregate(cd, ["on_ms", "noise_mode", "method"])
        cr = aggregate(cd, ["on_ms", "record", "noise_mode", "method"])
        cs.to_csv(out/"confirmation_summary.csv", index=False)
        cr.to_csv(out/"confirmation_per_record.csv", index=False)
        shown = cr[cr.noise_mode == "paired_sample20"].pivot(index=["on_ms", "record"], columns="method", values="errors")
        shown = shown[["first_est", "strong_est", "first_true_scale", "pdp_true"]].reset_index()
        shown["bits_per_row"] = 320
        cp = cs[(cs.on_ms == 4) & (cs.noise_mode == "paired_sample20")].set_index("method")
        record4 = cr[(cr.on_ms == 4) & (cr.record == 4) & (cr.noise_mode == "paired_sample20") & (cr.method == "strong_est")].iloc[0]
        confirmation_section = ("## 已开发记录004–005的固定确认\n\n"
            "这两条也是此前使用过的开发/训练记录，不是新独立留出。只固定比较reset波形、4/8 ms无保护和既定定时方案，没有调参数。\n\n"
            + md_table(shown) + "\n\n"
            + f"4 ms的原首到达窗口为 {int(cp.loc['first_est','errors'])}/640（{cp.loc['first_est','ber']:.2%}），"
            + f"前导主峰窗口为 {int(cp.loc['strong_est','errors'])}/640（{cp.loc['strong_est','ber']:.2%}）。"
            + f"但记录004的主峰窗口仍有 {int(record4.errors)}/320（{record4.ber:.2%}）错误，"
            + "说明改定时没有消除多径等残余失真，不能宣称问题全部解决。两种噪声及全部方法仍保留在确认统计CSV中。\n")
        confirmation_sources = [confirmation_path/name for name in ["metadata.json", "decoder_results.csv", "events.csv"]]
    main_conclusion = (f"保持收到的波形、估计尺度和精确双音检测不变，只把原来的首到达估计换为前导主峰定时，"
                       f"总错误由 **{int(first.errors)}/960（{first.ber:.6f}）** 变为 **{int(strong.errors)}/960（{strong.ber:.6f}）**。"
                       f"仅换为真实 bulk 尺度、仍用原首到达估计时则为 **{int(scale.errors)}/960（{scale.ber:.6f}）**。")
    if strong.errors < first.errors and first.errors-strong.errors > max(0, first.errors-scale.errors):
        main_conclusion += "在这些开发样本上，前导定时选择造成的差异明显大于单独更换 bulk 尺度；这是该配置高 BER 的重要来源，不能当作唯一原因或原论文实现的证明。"
    report = f"""# NOF1 4 ms 无保护偏差诊断

本报告主诊断使用开发记录 **001–003**，追加确认使用开发记录 **004–005**，没有读取留出 031–040；冻结接收器和策略没有修改。旧 `paper_alignment` 的主基线已逐包复核：15 包、960 bit、**373 个错误，BER=0.388542**；原论文表 2 对应 NOF1、4 ms、无保护、20 dB 标签为 **0.1137**。两者噪声、前导和接收实现并未完全对齐，不能仅据数值距离判断是否正确。

{main_conclusion}

{mechanism}

## 数据和回归检查

- 修正后的运行共 {len(events)} 个波形/噪声事件、{len(decoders)} 个解码对照和 {len(sweep)} 个定时扫描点。每个汇总条件均为 3 条记录 × 5 包 × 64 bit = 960 bit。
- 原始四组时长/保护控制的 **60 包**逐包错误数、包错误、载荷长度、噪声方差、初相和起点均与旧记录一致；其中 4 ms/无保护的 15 包恢复为 373/960。
- 同一事件的七种解码方法使用同一个收到的波形。`first_est` 的逐比特输出在 MATLAB 中与冻结接收器硬断言一致；精确双音和 FFT 方法使用相同插值后的判决采样块。
- reset/continuous/频率变体保留相同载荷随机种子、前导初相、未知起点、标准噪声种子及**同一噪声方差**。新增 `payload_reset_zero` 仅把载荷符号重置相位置零，前导与导频保持原随机初相。
- design 2–6 的有效符号加保护间隔均对应两音整数周期，reset 和 continuous 理应等价。本报告已核对150组事件的初相/噪声完全相同、接收能量/尺度近似相同，以及所有方法错误数相同。

先前 `pi` 循环变量遮蔽数学常量导致的运行已作废归档，本报告不使用其数据。当前运行用 `phaseIndex` 并通过上述旧基线及相位不变性核验。

## 同一波形上的定时、尺度、检测器

`strong_est` 在本报告称为**前导主峰定时**：它是接收器已有的 HFM 联合检测统计量局部峰，不保证是物理真最强路径，也不保证对两个载荷频率都是最优路径。首到达和主峰对照均只使用接收前导；本轮只是诊断替换，没有部署到冻结接收器。

`first_true_scale`/`strong_true_scale` 使用真实恢复的 bulk 时间尺度，仅作 oracle 诊断。实测 TVIR 还含残余多普勒，因此 `alpha_hat` 与 bulk 尺度不同不自动等于估计器错误。`pdp_true` 更使用已知发射起点和全记录平均 PDP 峰，无法直接部署，也不是保证的性能上界。

{md_table(primary_table)}

![同波形消融](<{figures[0].as_posix()}>)

在相同15包上，相对首到达基线的逐包变化如下。负的 `delta_errors_vs_first` 表示错误减少；同时列出改善、恶化和不变的包数，避免只看总和。

{md_table(primary_change)}

逐记录 BER 保留下表；只有3条开发记录，不把960个相关 bit 当作独立海况样本计算总体显著性。

{md_table(primary_records)}

这三条记录的定时位置和尺度概况如下。这里的参考是**平均 PDP 峰**，并非无歧义的真首径起点。

{md_table(timing_shown)}

FFT 对照使用无补零、矩形窗 N 点 FFT 的最近一个频点，非论文未披露的窄带积分实现。本组频率13,875/14,125 Hz在4 ms下为55.5/56.5个周期，最近正频率 bin 为56/57；精确双音方法仍可使用正交的四个实数正交模板，**半整数FFT bin不等于信号不正交**。同采样块的最近-bin结果为 {int(fft.errors)}/960（{fft.ber:.6f}），其差异只描述当前检测器假设，不能反推作者代码。

## 符号相位与固定频率诊断

连续相位按当前频率在“有效符号+保护段”期间累积，保护段仍发静默；reset则每个符号从该包固定初相重新开始。对13,875/14,125 Hz、4 ms且无保护，每步是半整数周期，因此两种策略的符号边界不同；其余整数周期控制应保持等价。

{md_table(phase_table)}

![相位和频率诊断](<{figures[1].as_posix()}>)

另两组频率13,750/14,000 Hz及14,000/14,250 Hz事先固定并全部列出。改变频率也改变实测信道的频率选择性衰落和已知单音导频，不是纯FFT网格实验。所有变体的噪声方差沿用同一时长/保护设置的原始reset参考，因此实际接收功率SNR可随波形变化；没有重新归一化噪声去制造公平性假象，也没有按与原论文接近程度选频率或相位。

主基线的噪声是旧 `sample_snr_packet` 20 dB定义：整帧接收能量除以含128 ms尾段的帧时长所得样本平均功率，再除以100。它不是 `Eb/N0`，不能与原论文未披露的SNR实现直接等同。无加噪声对照也完整保留。

## 定时扫描仅用于解释敏感性

对原始reset波形的4/8 ms无保护条件，以全记录平均 PDP峰为中心扫描±4 ms、步长0.25 ms，使用真实bulk尺度及共同20 dB噪声。图中展示全部33个偏移；没有用载荷BER挑选一个最低点加入接收器或策略。扫描需要真值评分及oracle参考位置，只能说明当前链路的定时敏感性。

![定时扫描](<{figures[2].as_posix()}>)

{confirmation_section}

## 全部预设配置

以下列出全部13个波形配置、两种噪声情况和七种方法，表中为BER；每格960 bit。`on/guard ms | f0/f1 | phase`定义行标签。完整计数见 {link('全部条件统计', out/'summary_all_conditions.csv')}；逐记录计数见 {link('逐记录统计', out/'per_record_summary.csv')}。

{chr(10).join(variant_tables)}

这次诊断能说明实现选择如何改变当前开发样本的误码，不能证明已经严格复现论文，也不构成留出验证后的接收机改进。所有冻结源码与策略保持不变；已知范围是NOF1三条主诊断及两条开发确认记录、短载荷和本项目公开帧格式。原文来源：[2024论文](https://vibsys.put.poznan.pl/_journal/2024-35-2/articles/vps_2024207.pdf)。

{link('旧基线逐包核对', out/'baseline_regression.csv')}、{link('整数周期相位不变性', out/'integer_stride_phase_invariants.csv')}、{link('全量定时扫描统计', out/'timing_sweep_summary.csv')} 和 {link('报告核验收据', out/'report_validation.json')}保留检查依据。
"""
    (out/"report.md").write_text(report, encoding="utf-8")
    source_files = paths+[metadata_path, old_path, ROOT/"matlab/run_nof1_diagnostic_v3.m",
                          ROOT/"matlab/+uac2/receive_packet.m", ROOT/"matlab/+uac2/default_config.m", Path(__file__)]
    source_files += confirmation_sources
    if toy_path.exists(): source_files.append(toy_path)
    receipt = {
        "status": "PASS", "records": [1, 2, 3], "heldout_records_used": [],
        "events": len(events), "decoder_rows": len(decoders), "sweep_rows": len(sweep),
        "original_regression_packets": 60, "primary_baseline_errors": 373, "primary_baseline_bits": 960,
        "primary_main_peak_errors": int(strong.errors), "primary_true_scale_errors": int(scale.errors),
        "primary_pdp_oracle_errors": int(pdp.errors), "primary_fft_errors": int(fft.errors),
        "integer_stride_phase_pairs": len(invariants), "shared_variance_phase_start_verified": True,
        "confirmation_records": [4, 5] if confirmation_sources else [], "confirmation_is_independent_holdout": False,
        "all_prespecified_conditions_reported": True, "timing_sweep_not_selected_for_deployment": True,
        "figures": [path.as_posix() for path in figures],
        "sha256": {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files},
    }
    (out/"report_validation.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in receipt.items() if key != "sha256"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
