"""Score the benchmark result files, and check the published scoring against them.

Every number in `../../experiments/real_data/README.md` comes out of this
script.  It reads the stored result files and recomputes (ADD, FAR, DR) under
two rules:

*published*
    a trial counts as a detection whenever a finite delay was recorded and the
    detector did not stop before the change point.

*protocol*
    the rule the study states in words -- a trial ends in a false alarm if the
    detector stops before the change point, in a detection if it stops after it,
    and **in a missed target if no alarm is raised before the end of the test
    segment**.

The two disagree because the detector used to report `detection_time =
len(test_segment)` when no alarm ever fired, and the harness scored that as a
detection at the last sample.  The baseline wrappers reported `None` in the same
situation and were scored correctly, so the discrepancy is confined to the GSA
rows.  The *protocol* column is the correct one.

Nothing is re-run and nothing is downloaded: the result files carry every trial.

    python experiments/real_data/score.py
    python experiments/real_data/score.py --set refit
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

HERE = Path(__file__).resolve().parent

TIER1_ORDER = [
    ("US RealInt", "realint_modern"),
    ("SKAB", "skab"),
    ("NASA IMS", "nasa_ims_3bearings"),
    ("TCPD", "tcpd"),
    ("NAB EC2", "nab"),
    ("NSL-KDD", "nslkdd_modern"),
]

REFIT_ORDER = [
    ("TCPD", "tcpd"),
    ("NAB EC2", "nab"),
    ("SKAB", "skab"),
    ("NSL-KDD", "nslkdd"),
]


def _finite_delay(record: dict[str, Any]) -> float | None:
    delay = record.get("delay")
    if delay is None:
        return None
    delay = float(delay)
    return delay if np.isfinite(delay) and delay < 1e100 else None


def classify_published(trial: dict, name: str) -> tuple[str, float | None]:
    record = trial["detectors"].get(name)
    if record is None or "error" in record:
        return "error", None
    if record.get("is_false_alarm"):
        return "fa", None
    delay = _finite_delay(record)
    return ("detect", delay) if delay is not None else ("miss", None)


def classify_protocol(trial: dict, name: str) -> tuple[str, float | None]:
    """`alarm_raised` decides where it exists; otherwise the end-of-segment tell.

    Files written before the flag existed cannot say directly whether an alarm
    fired.  For those, a recorded stop at or after the last index of the test
    segment can only have come from the end-of-segment fallback, so it is a
    missed target.
    """
    record = trial["detectors"].get(name)
    if record is None or "error" in record:
        return "error", None
    if record.get("is_false_alarm"):
        return "fa", None

    alarm = record.get("diagnostics", {}).get("alarm_raised", record.get("alarm_raised"))
    stop = record.get("detection_time")
    if alarm is False:
        return "miss", None
    if alarm is None and (stop is None or stop >= trial["test_length"] - 1):
        return "miss", None

    delay = _finite_delay(record)
    return ("detect", delay) if delay is not None else ("miss", None)


def score(trials: list[dict], name: str, classify: Callable) -> dict[str, Any]:
    delays: list[float] = []
    false_alarms = missed = errors = 0
    for trial in trials:
        kind, delay = classify(trial, name)
        if kind == "error":
            errors += 1
        elif kind == "fa":
            false_alarms += 1
        elif kind == "miss":
            missed += 1
        else:
            delays.append(delay)

    n = len(trials)
    return {
        "add": float(np.mean(delays)) if delays else float("inf"),
        "far": false_alarms / n if n else 0.0,
        "detection_rate": len(delays) / n if n else 0.0,
        "n_detected": len(delays),
        "n_false_alarms": false_alarms,
        "n_missed": missed,
        "n_errors": errors,
        "n_trials": n,
    }


def _fmt(metrics: dict[str, Any]) -> str:
    add = metrics["add"]
    add_s = f"{add:7.1f}" if np.isfinite(add) else "    inf"
    return f"({add_s}, {metrics['far']:.3f}, {metrics['detection_rate']:.2f})"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--set", choices=["published", "refit"], default="published",
                        help="'published' scores the files behind the original "
                             "table; 'refit' scores the standardised re-run")
    parser.add_argument("--out", type=Path, default=None,
                        help="Optional JSON output path")
    args = parser.parse_args()

    directory = HERE / "results" / args.set
    order = TIER1_ORDER if args.set == "published" else REFIT_ORDER
    report: dict[str, Any] = {"set": args.set, "datasets": {}}

    for label, stem in order:
        path = directory / f"{stem}.json"
        if not path.exists():
            print(f"[skip] {label}: {path.name} not present")
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        trials = payload["per_trial"]
        names = sorted(trials[0]["detectors"].keys())

        print(f"\n=== {label}  ({len(trials)} trials)")
        print(f"{'detector':<16} {'published rule':<25} {'protocol rule':<25} "
              f"{'silent':>7}")

        per_detector = {}
        for name in names:
            published = score(trials, name, classify_published)
            protocol = score(trials, name, classify_protocol)
            silent = published["n_detected"] - protocol["n_detected"]
            per_detector[name] = {"published": published, "protocol": protocol,
                                  "silent_runs": silent}
            flag = "  <--" if silent else ""
            print(f"{name:<16} {_fmt(published):<25} {_fmt(protocol):<25} "
                  f"{silent:>7}{flag}")

        report["datasets"][label] = {"n_trials": len(trials),
                                     "detectors": per_detector}

    if args.out:
        args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"\nWritten: {args.out}")


if __name__ == "__main__":
    main()
