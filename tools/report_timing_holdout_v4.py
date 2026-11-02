"""Verify and summarize the preregistered, paired V4 timing experiment.

No raw channel files are required for report verification. This reports
within-record comparisons, not an independent-environment significance test.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/v4/timing_holdout"
CHANNELS = ("NOF1", "NCS1")
DESIGN = {1: (4, 0), 2: (4, 4), 3: (8, 0), 4: (8, 8)}
METHODS = ("first", "strong")
INTS = ("record", "slot", "design", "on_ms", "guard_ms", "payload_bits", "decoded_bits",
        "bit_errors", "bit_erasures", "packet_error", "detected", "received_samples")
FLOATS = ("preamble_start_sample", "alpha_hat", "sync_score", "delay_spread_ms",
          "carrier_phase_rad", "injected_start_s", "alpha_true_simulation_only", "noise_variance")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_freeze():
    freeze_path = OUT / "freeze.json"
    freeze = read_json(freeze_path)
    require(bool(freeze.get("files")), "Missing freeze source list")
    for entry in freeze["files"]:
        path = (ROOT / entry["path"]).resolve()
        require(path.is_relative_to(ROOT), "Frozen path escapes repository")
        require(sha256(path) == entry["sha256"], f"Frozen file changed: {entry['path']}")
    return sha256(freeze_path), len(freeze["files"])


def read_packets(out):
    with (out / "per_packet.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for name in INTS:
            number = float(row[name])
            require(math.isfinite(number) and number == int(number), f"Noninteger {name}")
            row[name] = int(number)
        for name in FLOATS:
            row[name] = float(row[name])
    return rows


def validate(out, development):
    metadata = read_json(out / "metadata.json")
    records = (1,) if development else tuple(range(11, 21))
    require(metadata["status"] == "COMPLETE", "Refuse incomplete experiment")
    require(metadata["stage"] == ("development" if development else "final"), "Wrong stage")
    actual_records = metadata["record_ids"]
    if isinstance(actual_records, int):
        actual_records = [actual_records]
    require(tuple(actual_records) == records, "Wrong predeclared records")
    require(metadata["channels"] == list(CHANNELS), "Wrong channels")
    require(metadata["design"] == [list(DESIGN[k]) for k in DESIGN], "Wrong design")
    require(metadata["payload_bits"] == 64 and metadata["slots_per_record"] == 5, "Packet design")
    require(metadata["horizon_s"] == 32 and metadata["slot_duration_s"] == 6.4, "Window design")
    require(metadata["seed"] == 2026110205, "Random seed")
    require(metadata["snr_label_db"] == 20 and metadata["noise_mode"] == "sample_snr_packet", "Noise design")
    require(metadata["phase_mode"] == "reset" and metadata["methods"] == list(METHODS), "Phase/method design")
    require(metadata["reserved_unread_record_ids"] == list(range(21, 31)), "Reserved records")
    freeze_count = 0
    if not development:
        freeze_hash, freeze_count = check_freeze()
        require(metadata["freeze_sha256"] == freeze_hash, "Run used a different freeze manifest")
        require(metadata["frozen_sources_unchanged"] is True, "Frozen code not verified after run")
    else:
        require(metadata["development_v3_checks"] == 40, "Missing V3 regressions")
    rows = read_packets(out)
    expected = {(c, r, s, d, m) for c in CHANNELS for r in records
                for s in range(1, 6) for d in DESIGN for m in METHODS}
    keys = [(r["channel"], r["record"], r["slot"], r["design"], r["method"]) for r in rows]
    require(len(rows) == len(expected) == len(set(keys)), "Incomplete or duplicate packet grid")
    require(set(keys) == expected, "Unexpected packet grid")
    require(metadata["method_rows"] == len(rows) and metadata["waveform_count"] == len(rows)//2, "Metadata count")
    pairs = defaultdict(dict)
    sources = {entry["path"]: entry["sha256"] for entry in metadata["sources"]}
    require(len(sources) == len(CHANNELS)*len(records), "Channel source inventory")
    for row in rows:
        require((row["on_ms"], row["guard_ms"]) == DESIGN[row["design"]], "Wrong timing condition")
        require(row["payload_bits"] == 64, "Payload size")
        require(row["detected"] in (0, 1), "Detection flag")
        require(row["decoded_bits"] == 64*row["detected"], "Missing payload not erased")
        require(row["bit_erasures"] == 64-row["decoded_bits"], "Erasure accounting")
        require(0 <= row["bit_errors"] <= row["decoded_bits"], "Error accounting")
        require(row["packet_error"] == int(row["bit_errors"]+row["bit_erasures"] > 0), "PER accounting")
        require(0 <= row["carrier_phase_rad"] < 2*math.pi, "Invalid initial phase")
        require(math.isfinite(row["noise_variance"]) and row["noise_variance"] >= 0, "Noise variance")
        require(row["received_samples"] > 0, "Empty waveform")
        require(len(row["rx_sha256"]) == 64 and all(c in "0123456789abcdef" for c in row["rx_sha256"]), "Waveform digest")
        path = f"vendor/Watermark/input/channels/{row['channel']}/mat/{row['channel']}_{row['record']:03}.mat"
        require(sources.get(path) == row["source_sha256"], "Channel source digest")
        if row["detected"]:
            require(math.isfinite(row["preamble_start_sample"]) and math.isfinite(row["alpha_hat"]), "Acquired timing/scale")
        pairs[(row["channel"], row["record"], row["slot"], row["design"])][row["method"]] = row
    for pair in pairs.values():
        a, b = pair["first"], pair["strong"]
        for name in ("rx_sha256", "source_sha256", "carrier_phase_rad", "noise_variance",
                     "injected_start_s", "alpha_true_simulation_only", "received_samples"):
            require(a[name] == b[name], f"Broken same-waveform pairing: {name}")
        require(a["alpha_hat"] == b["alpha_hat"] or
                (math.isnan(a["alpha_hat"]) and math.isnan(b["alpha_hat"])), "Different estimated scales")
    return metadata, rows, freeze_count


def aggregate(rows, keys):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in keys)].append(row)
    output = []
    for values, group in sorted(groups.items()):
        result = dict(zip(keys, values))
        for output_name, source in (("bits", "payload_bits"), ("decoded_bits", "decoded_bits"),
                                    ("errors", "bit_errors"), ("erasures", "bit_erasures"),
                                    ("packet_errors", "packet_error"), ("detected_packets", "detected")):
            result[output_name] = sum(row[source] for row in group)
        result["packets"] = len(group)
        result["bit_loss_rate"] = (result["errors"]+result["erasures"])/result["bits"]
        result["ber_decoded"] = result["errors"]/result["decoded_bits"] if result["decoded_bits"] else None
        result["erasure_rate"] = result["erasures"]/result["bits"]
        result["per"] = result["packet_errors"]/result["packets"]
        output.append(result)
    return output


def pair_records(records):
    groups = defaultdict(dict)
    for row in records:
        groups[(row["channel"], row["record"], row["on_ms"], row["guard_ms"])][row["method"]] = row
    output = []
    for (channel, record, on_ms, guard_ms), pair in sorted(groups.items()):
        a, b = pair["first"], pair["strong"]
        output.append(dict(channel=channel, record=record, on_ms=on_ms, guard_ms=guard_ms,
                           bits=a["bits"], first_errors=a["errors"], strong_errors=b["errors"],
                           first_erasures=a["erasures"], strong_erasures=b["erasures"],
                           first_loss=a["bit_loss_rate"], strong_loss=b["bit_loss_rate"],
                           delta_loss=b["bit_loss_rate"]-a["bit_loss_rate"],
                           delta_packet_errors=b["packet_errors"]-a["packet_errors"]))
    return output


def paired_summary(pairs):
    groups = defaultdict(list)
    for row in pairs:
        groups[(row["channel"], row["on_ms"], row["guard_ms"])].append(row)
    return [dict(channel=key[0], on_ms=key[1], guard_ms=key[2], records=len(group),
                 improved_records=sum(r["delta_loss"] < 0 for r in group),
                 tied_records=sum(r["delta_loss"] == 0 for r in group),
                 worse_records=sum(r["delta_loss"] > 0 for r in group),
                 min_delta_loss=min(r["delta_loss"] for r in group),
                 median_delta_loss=median(r["delta_loss"] for r in group),
                 max_delta_loss=max(r["delta_loss"] for r in group))
            for key, group in sorted(groups.items())]


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def percent(value):
    return "—" if value is None else f"{value*100:.2f}%"


def plot(out, pairs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, channel in zip(axes, CHANNELS):
        selected = [r for r in pairs if r["channel"] == channel and r["on_ms"] == 4 and r["guard_ms"] == 0]
        x = [r["record"] for r in selected]
        ax.plot(x, [r["first_loss"]*100 for r in selected], "o-", color="#bd6748", label="First arrival")
        ax.plot(x, [r["strong_loss"]*100 for r in selected], "s-", color="#187c8c", label="HFM main peak")
        ax.set_title(f"{channel} | 4 ms, no guard")
        ax.set_xlabel("Previously unused record ID")
        ax.set_xticks(x)
        ax.grid(alpha=.2)
        ax.legend(frameon=False)
    axes[0].set_ylabel("Bit loss (errors + erasures), %")
    fig.suptitle("Frozen timing rules; same waveform within each pair", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "primary_by_record.png", dpi=180)
    plt.close(fig)


def report(out, metadata, summary, pairs, comparisons, development):
    stage = "开发自检（001）" if development else "新记录验证（011–020）"
    docs_link = "../../../../docs/timing_holdout_protocol_v4.md" if development else "../../../docs/timing_holdout_protocol_v4.md"
    freeze_link = "../freeze.json" if development else "freeze.json"
    by_key = {(r["channel"], r["on_ms"], r["guard_ms"], r["method"]): r for r in summary}
    a = by_key[("NOF1", 4, 0, "first")]
    b = by_key[("NOF1", 4, 0, "strong")]
    direction = "下降" if b["bit_loss_rate"] < a["bit_loss_rate"] else "未下降"
    worse_conditions = []
    for channel in CHANNELS:
        for on_ms, guard_ms in DESIGN.values():
            first = by_key[(channel, on_ms, guard_ms, "first")]
            strong = by_key[(channel, on_ms, guard_ms, "strong")]
            if strong["bit_loss_rate"] > first["bit_loss_rate"]:
                worse_conditions.append(f"{channel} {on_ms}/{guard_ms} ms（{percent(first['bit_loss_rate'])} → {percent(strong['bit_loss_rate'])}）")
    degradation = ("合计比特损失率恶化的预定条件："+"；".join(worse_conditions)+"。主峰规则不能据本轮结果作为全条件统一改进。") if worse_conditions else "本批预定条件中，没有条件的合计比特损失率因主峰规则增加；仍应查看逐记录恶化数量，不能推断所有包或其他环境均受益。"
    lines = [f"# V4 主峰定时 {stage}", "",
             f"预定主条件 NOF1、4 ms、无保护：合计比特损失率从 **{percent(a['bit_loss_rate'])}** 到 **{percent(b['bit_loss_rate'])}**，{direction}。",
             f"首到达规则为 {a['errors']} 错误 + {a['erasures']} 擦除 / {a['bits']} bit；主峰规则为 {b['errors']} 错误 + {b['erasures']} 擦除 / {b['bits']} bit。",
             "", "本轮只改有效符号积分起点；导频估计、检测器、每个配对的波形、噪声与初相相同。",
             "", "## 全部预定条件", "",
             "|信道|发声/保护 ms|方法|错误/已解码 bit|擦除/总 bit|条件 BER|比特损失率|错误包/总包|",
             "|---|---|---|---|---|---|---|---|"]
    for row in summary:
        label = "首到达" if row["method"] == "first" else "主峰"
        lines.append(f"|{row['channel']}|{row['on_ms']}/{row['guard_ms']}|{label}|{row['errors']}/{row['decoded_bits']}|{row['erasures']}/{row['bits']}|{percent(row['ber_decoded'])}|{percent(row['bit_loss_rate'])}|{row['packet_errors']}/{row['packets']}|")
    lines += ["", "条件 BER 只使用已解码比特；比特损失率将未捕获或不完整的整包计为 64 bit 擦除。无擦除时两者相同。PER 由错误包/总包给出。", "",
              "## 逐记录配对差异", "",
              "差异定义为主峰减首到达，负数表示改善；每条记录先聚合 5 包。", "",
              "|信道|发声/保护 ms|改善/持平/恶化记录|差异最小值|差异中位数|差异最大值|",
              "|---|---|---|---|---|---|"]
    for row in comparisons:
        lines.append(f"|{row['channel']}|{row['on_ms']}/{row['guard_ms']}|{row['improved_records']}/{row['tied_records']}/{row['worse_records']}|{row['min_delta_loss']*100:+.2f} pp|{row['median_delta_loss']*100:+.2f} pp|{row['max_delta_loss']*100:+.2f} pp|")
    lines += ["", "![4 ms 无保护：逐记录比特损失率](primary_by_record.png)", "",
              "## 可支持的结论与边界", "",
              degradation, "",
              "- 这是冻结规则后的新记录评测；011–020 未参与此前开发，021–030 继续保留。001–010 和旧评测 031–040 不参与最终统计。" if not development else
              "- 这是已使用的开发记录 001 的实现自检，不是新记录泛化证据；NOF1 20 个波形、40 条方法观测与原 V3 错误数、起点和尺度逐项回归。",
              "- 同一 Watermark 信道序列的相邻记录可能相关。表中的范围和中位数描述本批记录的离散程度，不是独立海域置信区间；不把比特当独立统计重复。",
              "- 主条件改善只支持本条件、本噪声定义和这些记录；次要条件若恶化，必须保留并限定适用范围，不能据主条件平均数自动替换所有信道的接收器。",
              "- 未用本轮结果修改定时参数、记录选择或自适应策略；没有新增硬件、功耗、CRC 或连续排程结论。",
              "", "## 可追溯文件", "",
              f"- [预注册]({docs_link})；[冻结清单]({freeze_link})。",
              "- [逐包方法结果](per_packet.csv)、[逐记录聚合](per_record.csv)、[逐记录配对](paired_records.csv)、[条件汇总](summary.csv)。",
              "- [运行元数据与数据 SHA-256](metadata.json)、[结果核验](verification.json)。",
              f"- 收到波形数：{metadata['waveform_count']}；方法观测数：{metadata['method_rows']}。每个事件的两行拥有相同收到波形 SHA-256。",
              f"- 复跑：`addpath('matlab'); run_timing_holdout_v4('{metadata['stage']}');`（需要 JVM 和预先下载的数据）；随后 `python tools/report_timing_holdout_v4.py{' --development' if development else ''}`。", ""]
    # Relative docs link from results/v4/timing_holdout is three parents.
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--development", action="store_true")
    parser.add_argument("--check-freeze", action="store_true")
    args = parser.parse_args()
    if args.check_freeze:
        digest, count = check_freeze()
        print(f"PASS: {count} frozen files; manifest {digest}")
        return
    out = OUT / "development" if args.development else OUT
    metadata, rows, freeze_count = validate(out, args.development)
    summary = aggregate(rows, ("channel", "on_ms", "guard_ms", "method"))
    records = aggregate(rows, ("channel", "record", "on_ms", "guard_ms", "method"))
    pairs = pair_records(records)
    comparisons = paired_summary(pairs)
    for name, data in (("summary", summary), ("per_record", records),
                       ("paired_records", pairs), ("paired_summary", comparisons)):
        write_csv(out / f"{name}.csv", data)
    plot(out, pairs)
    verification = dict(status="PASS", stage=metadata["stage"], method_rows=len(rows),
                        waveform_pairs=len(rows)//2, frozen_files_checked=freeze_count,
                        same_waveform_pairs_verified=True, all_erasures_retained=True,
                        record_grid_verified=True, conditional_ber_and_total_loss_separate=True,
                        source_files={name: sha256(out / name) for name in ("per_packet.csv", "metadata.json")},
                        scope="Count, grid, pairing and frozen-code integrity; no claim of independent environments")
    (out / "verification.json").write_text(json.dumps(verification, indent=2)+"\n", encoding="utf-8")
    report(out, metadata, summary, pairs, comparisons, args.development)
    print(f"PASS: {len(rows)} method observations / {len(rows)//2} waveform pairs; erasures and source hashes retained.")
    for r in summary:
        if r["on_ms"] == 4 and r["guard_ms"] == 0:
            print(f"{r['channel']} {r['method']}: {r['errors']} errors + {r['erasures']} erasures / {r['bits']} bit; loss {percent(r['bit_loss_rate'])}")


if __name__ == "__main__":
    main()
