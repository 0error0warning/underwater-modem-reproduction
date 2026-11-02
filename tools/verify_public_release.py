"""Verify the published artifact using repository files only (Python stdlib)."""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(relative):
    return json.loads((ROOT / relative).read_text(encoding="utf-8-sig"))


def read_csv(relative):
    with (ROOT / relative).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row, name):
    value = float(row[name])
    require(math.isfinite(value), f"Nonfinite {name}")
    return value


def check_hash(relative, expected):
    path = (ROOT / relative).resolve()
    require(path.is_relative_to(ROOT), f"Path escapes repository: {relative}")
    require(path.is_file(), f"Missing file: {relative}")
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    require(actual == expected, f"SHA-256 mismatch: {relative}")


def check_diagnostic(relative, records, expected_errors, expected_bits):
    selected = [r for r in read_csv(relative) if
                r["channel"] == "NOF1" and number(r, "design") == 1 and
                number(r, "on_ms") == 4 and number(r, "guard_ms") == 0 and
                r["phase_mode"] == "reset" and r["noise_mode"] == "paired_sample20"]
    event_sets = []
    for method, expected in zip(("first_est", "strong_est"), expected_errors):
        rows = [r for r in selected if r["method"] == method]
        events = {(int(r["record"]), int(r["slot"])) for r in rows}
        require(len(rows) == len(records)*5 == len(events), f"{relative}: event grid")
        require({int(r["record"]) for r in rows} == records, f"{relative}: records")
        require(sum(number(r, "payload_bits") for r in rows) == expected_bits, f"{method}: bits")
        require(sum(number(r, "bit_errors") for r in rows) == expected, f"{method}: errors")
        event_sets.append(events)
    require(event_sets[0] == event_sets[1], "Diagnostic methods do not share events")


def check_continuous():
    rows = read_csv("results/v3/continuous_phase_corrected/records_continuous.csv")
    require(len(rows) == 595, "Continuous packet count")
    groups = defaultdict(list)
    for row in rows:
        groups[(row["policy"], row["channel"], int(row["record"]))].append(row)
    policies = {"best_fixed": 11.6, "adaptive_stump": 9.8}
    expected = {(p, c, r) for p in policies for c in ("NOF1", "NCS1") for r in range(31, 41)}
    require(set(groups) == expected and len(groups) == 40, "Continuous record/policy grid")
    bits, seconds, phases = defaultdict(float), defaultdict(float), defaultdict(dict)
    for (policy, channel, record), group in groups.items():
        group.sort(key=lambda r: int(r["packet"]))
        previous_end = 0
        for packet, row in enumerate(group, 1):
            start, end = number(row, "start_s"), number(row, "next_start_s")
            require(int(row["packet"]) == packet and number(row, "horizon_s") == 32, "Packet/horizon")
            require(start >= previous_end-1e-8 and end <= 32+1e-8, "Schedule outside horizon")
            require(math.isclose(end-start, number(row, "occupied_s")+number(row, "control_s"), abs_tol=1e-8), "Time budget")
            phase = number(row, "carrier_phase_rad")
            require(0 <= phase < 2*math.pi, "Invalid carrier phase")
            phases[(channel, record, packet)][policy] = phase
            good = number(row, "detected") == 1 and number(row, "bit_errors") == 0 and number(row, "bit_erasures") == 0
            require(number(row, "packet_error") == int(not good), "Packet error accounting")
            bits[policy] += number(row, "payload_bits") if good else 0
            previous_end = end
        seconds[policy] += 32  # One complete observation horizon per record, not per packet.
    for policy, expected_rate in policies.items():
        require(seconds[policy] == 640, f"{policy}: full denominator")
        require(math.isclose(bits[policy]/seconds[policy], expected_rate, abs_tol=1e-12), f"{policy}: rate")
    common = [value for value in phases.values() if len(value) == 2]
    require(len(common) == 255, "Common packet phase count")
    require(all(p["best_fixed"] == p["adaptive_stump"] for p in common), "Unequal common phases")


def main():
    files = read_json("publication_manifest.json")["files"]
    require(bool(files) and "publication_manifest.json" not in files, "Invalid manifest")
    for relative, entry in files.items():
        check_hash(relative, entry["sha256"])
    for entry in read_json("results/v2/receiver_code_freeze.json")["files"]:
        check_hash(entry["path"], entry["sha256"])
    check_diagnostic("results/v3/nof1_diagnostic/decoder_results.csv", {1, 2, 3}, (373, 24), 960)
    check_diagnostic("results/v3/nof1_diagnostic/confirmation/decoder_results.csv", {4, 5}, (163, 59), 640)
    check_continuous()
    print(f"PASS: {len(files)} manifest files; frozen MATLAB hashes; diagnostic counts; 595 packets / 40 horizons; 11.6 / 9.8 bit/s; 255 identical phase pairs.")


if __name__ == "__main__":
    main()
