"""Independently reconcile V4 P0 artifacts; Python standard library only."""
import csv
import binascii
import hashlib
import itertools
import json
import math
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8-sig"))


def read_csv(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def integer(row, field):
    value = float(row[field])
    require(math.isfinite(value) and value.is_integer(), f"Invalid {field}: {value}")
    return int(value)


def check_timing():
    base = "results/v4/timing_holdout/"
    meta = read_json(base + "metadata.json")
    freeze = read_json(base + "freeze.json")
    require(meta["status"] == "COMPLETE" and meta["stage"] == "final", "Incomplete timing run")
    require(meta["record_ids"] == list(range(11, 21)), "Changed evaluation records")
    require(meta["reserved_unread_record_ids"] == list(range(21, 31)), "Changed reserve")
    require(meta["waveform_count"] == 400 and meta["method_rows"] == 800, "Timing run size")
    require(meta["freeze_sha256"] == digest(base + "freeze.json"), "Freeze digest changed")
    require(meta["frozen_sources_unchanged"], "Frozen sources changed during execution")
    for item in freeze["files"]:
        require(digest(item["path"]) == item["sha256"], f"Frozen file changed: {item['path']}")
    require(datetime.fromisoformat(freeze["frozen_utc"].replace("Z", "+00:00")) <
            datetime.fromisoformat(meta["started_utc"].replace("Z", "+00:00")), "Freeze must precede evaluation")
    sources = {item["path"]: item["sha256"] for item in meta["sources"]}
    expected_sources = {f"vendor/Watermark/input/channels/{c}/mat/{c}_{r:03}.mat"
                        for c, r in itertools.product(("NOF1", "NCS1"), range(11, 21))}
    require(set(sources) == expected_sources, "Unexpected input record set")
    availability = freeze["input_availability_before_download"]
    require(set(availability) == expected_sources and not any(availability.values()),
            "Expected input absence receipt before download")
    downloads = {"vendor/" + item["name"]: item["sha256"]
                 for item in read_json("references/watermark_download_manifest.json")["files"]}
    require(all(downloads.get(path) == sha for path, sha in sources.items()), "Input hashes disagree")

    rows = read_csv(base + "per_packet.csv")
    expected = set(itertools.product(("NOF1", "NCS1"), range(11, 21), range(1, 6), range(1, 5), ("first", "strong")))
    actual, pairs, groups = set(), defaultdict(list), defaultdict(lambda: [0, 0, 0, 0, 0, 0])
    designs = {1: (4, 0), 2: (4, 4), 3: (8, 0), 4: (8, 8)}
    for row in rows:
        channel, record, slot, design, method = (row["channel"], integer(row, "record"),
                                                integer(row, "slot"), integer(row, "design"), row["method"])
        key = (channel, record, slot, design, method)
        require(key not in actual, f"Duplicate timing observation: {key}")
        actual.add(key)
        require(designs[design] == (integer(row, "on_ms"), integer(row, "guard_ms")), "Design changed")
        bits, decoded, errors, erased, packet_error, detected = (
            integer(row, field) for field in ("payload_bits", "decoded_bits", "bit_errors",
                                              "bit_erasures", "packet_error", "detected"))
        require(bits == 64 and decoded + erased == bits and 0 <= errors <= decoded,
                f"Invalid error/erasure accounting: {key}")
        require(detected in (0, 1) and decoded == 64 * detected, "Detection/decoded count mismatch")
        require(packet_error == int(errors + erased > 0), "Packet failure accounting mismatch")
        require(0 <= float(row["carrier_phase_rad"]) < 2 * math.pi, "Carrier phase outside range")
        source = f"vendor/Watermark/input/channels/{channel}/mat/{channel}_{record:03}.mat"
        require(row["source_sha256"] == sources[source], "Per-packet input hash mismatch")
        require(len(row["rx_sha256"]) == 64, "Missing waveform digest")
        pairs[key[:-1]].append(row)
        group = groups[(channel, design, method)]
        for index, value in enumerate((bits, decoded, errors, erased, packet_error, 1)):
            group[index] += value
    require(actual == expected and len(rows) == 800, "Missing or unexpected timing observations")
    for key, pair in pairs.items():
        require(len(pair) == 2, "Unpaired observation")
        for field in ("rx_sha256", "source_sha256", "carrier_phase_rad", "noise_variance",
                      "injected_start_s", "alpha_true_simulation_only", "received_samples", "alpha_hat"):
            require(pair[0][field] == pair[1][field], f"Unpaired {field}: {key}")
    primary = {method: groups[("NOF1", 1, method)] for method in ("first", "strong")}
    print("PASS: V4 frozen inputs/code; 400 paired waveforms; 800 observations; errors and erasures reconciled.")
    print("NOF1 4/0 primary (errors, erasures, total): " + "; ".join(
        f"{method}=({values[2]}, {values[3]}, {values[0]})" for method, values in primary.items()))
    return groups


def check_packet_link():
    base = "results/v4/packet_link/"
    validation = read_json(base + "validation.json")
    cases = read_csv(base + "validation_cases.csv")
    require(validation["status"] == "PASS" and validation["all_passed"], "Packet validation failed")
    require(len(cases) == validation["total_cases"] == validation["passed_cases"] == 35,
            "Packet validation count mismatch")
    require(len({row["name"] for row in cases}) == 35 and
            all(integer(row, "passed") == 1 for row in cases), "Missing/failed validation case")
    require(validation["crc_single_bit_trials"] == 256 and
            validation["noise_only_trials"] == 10 and validation["noise_only_false_accepts"] == 0,
            "CRC/noise test evidence mismatch")
    demo = read_json(base + "demo.json")
    require(demo["status"] == "PASS" and demo["all_passed"], "Packet demo failed")
    wire = bytes(demo["example_wire_bytes"])
    require(wire[:4] == bytes.fromhex("d3910101"), "Unexpected wire format")
    require(int.from_bytes(wire[6:8], "big") == len(wire) - 10 <= 64, "Wire length mismatch")
    require(wire[8:-2].decode("utf-8") == demo["message"], "Message was not recovered from wire bytes")
    require(binascii.crc_hqx(b"123456789", 0xFFFF) == 0x29B1, "CRC check value")
    require(binascii.crc_hqx(wire[:-2], 0xFFFF) == int.from_bytes(wire[-2:], "big"),
            "Independent Python CRC disagrees with MATLAB frame")
    scenarios = read_csv(base + "demo_scenarios.csv")
    events = read_csv(base + "demo_events.csv")
    require(len(scenarios) == demo["scenario_count"] == 5 and len(events) == demo["event_count"] == 11,
            "Demo trace size mismatch")
    expected = {"normal": (1, 1, 1), "lose_first_data": (2, 1, 1), "lose_first_ack": (2, 1, 1),
                "lose_all_data": (3, 0, 0), "lose_all_ack": (3, 0, 1)}
    require({row["scenario"] for row in scenarios} == set(expected), "Missing loss scenario")
    for row in scenarios:
        name = row["scenario"]
        attempts, confirmed, deliveries = expected[name]
        require(tuple(integer(row, field) for field in ("attempts", "sender_confirmed", "receiver_deliveries"))
                == expected[name], f"Wrong outcome: {name}")
        require(row["sender_state"] == ("idle" if confirmed else "exhausted"), "Wrong sender state")
        require(integer(row, "passed") == 1, "Demo scenario failed")
        require(row["received_text"] == (demo["message"] if deliveries else ""), "Application text mismatch")
        trace = sorted((event for event in events if event["scenario"] == name), key=lambda x: integer(x, "attempt"))
        require([integer(event, "attempt") for event in trace] == list(range(1, attempts + 1)), "Retry count mismatch")
        require(sum(integer(event, "application_delivered") for event in trace) == deliveries,
                "Duplicate application delivery or missing trace")
        require(sum(integer(event, "ack_crc_valid") for event in trace) == confirmed,
                "ACK confirmation not supported by decoded waveform")
        for event in trace:
            require(integer(event, "application_delivered") <= integer(event, "data_crc_valid")
                    <= integer(event, "data_reached_b"), "Data delivered without decoded CRC")
            require(integer(event, "ack_crc_valid") <= integer(event, "ack_reached_a"), "ACK accepted without waveform")
        if name in ("lose_first_ack", "lose_all_ack"):
            require(all(event["b_event"] == "duplicate_acknowledged" for event in trace[1:]),
                    "Retry after lost ACK not duplicate-suppressed")
    provenance = read_json(base + "provenance.json")
    for item in provenance["source_files"] + provenance["output_files"]:
        require(digest(item["path"]) == item["sha256"], f"Protocol evidence/source changed: {item['path']}")
    print("PASS: 35 protocol checks; independent UTF-8/length/CRC decoding; 5 loss scenarios / 11 transmissions; duplicate suppression and unconfirmed-delivery semantics.")


if __name__ == "__main__":
    check_timing()
    check_packet_link()
