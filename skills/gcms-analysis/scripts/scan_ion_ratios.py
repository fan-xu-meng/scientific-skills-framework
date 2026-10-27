#!/usr/bin/env python3
"""Compute raw per-scan ion-channel ratios from an exported long-form CSV."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def run(source: Path, target: Path, numerator: float, denominator: float, tolerance: float) -> None:
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"run", "scan", "time_min", "mz", "intensity"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("input requires run,scan,time_min,mz,intensity")
        rows = list(reader)
    grouped: dict[tuple[str, str, str], dict[str, float]] = defaultdict(dict)
    for row in rows:
        try:
            mz = float(row["mz"])
            intensity = float(row["intensity"])
        except ValueError as exc:
            raise ValueError(f"non-numeric m/z or intensity in {row!r}") from exc
        if not math.isfinite(mz) or not math.isfinite(intensity) or intensity < 0:
            raise ValueError(f"invalid m/z or intensity in {row!r}")
        key = (row["run"], row["scan"], row["time_min"])
        grouped[key]
        for label, target_mz in (("numerator", numerator), ("denominator", denominator)):
            if abs(mz - target_mz) <= tolerance:
                if label in grouped[key]:
                    raise ValueError(f"multiple ions match {label} for {key}")
                grouped[key][label] = intensity
    output = []
    for (run_id, scan, time_min), ions in sorted(grouped.items()):
        top, base = ions.get("numerator"), ions.get("denominator")
        status = "ok" if top is not None and base is not None and base > 0 else "missing_or_zero_channel"
        output.append({"run": run_id, "scan": scan, "time_min": time_min,
                       "numerator_mz": numerator, "denominator_mz": denominator,
                       "numerator_intensity": "" if top is None else top,
                       "denominator_intensity": "" if base is None else base,
                       "raw_ratio": top / base if status == "ok" else "", "qc": status})
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]) if output else [
            "run", "scan", "time_min", "numerator_mz", "denominator_mz",
            "numerator_intensity", "denominator_intensity", "raw_ratio", "qc"])
        writer.writeheader()
        writer.writerows(output)
    print(f"Wrote {len(output)} raw scan ratios; no isotope correction or molecular assignment applied")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--numerator-mz", required=True, type=float)
    parser.add_argument("--denominator-mz", required=True, type=float)
    parser.add_argument("--tolerance", default=0.5, type=float)
    args = parser.parse_args()
    if args.tolerance < 0 or args.numerator_mz == args.denominator_mz:
        parser.error("use distinct channels and a nonnegative tolerance")
    run(args.input, args.output, args.numerator_mz, args.denominator_mz, args.tolerance)


if __name__ == "__main__":
    main()
