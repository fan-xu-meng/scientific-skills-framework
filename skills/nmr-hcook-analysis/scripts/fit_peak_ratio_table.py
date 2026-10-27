#!/usr/bin/env python3
"""Fit reviewed NMR formate/internal-standard peak ratios to concentration."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def run(source: Path, target: Path) -> None:
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"sample", "formate_area", "internal_standard_area", "known_formate_mM", "qc"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("input requires " + ",".join(sorted(required)))
        rows = list(reader)
    standards: list[tuple[float, float]] = []
    ratios: list[float] = []
    for line, row in enumerate(rows, start=2):
        try:
            formate = float(row["formate_area"])
            reference = float(row["internal_standard_area"])
        except ValueError as exc:
            raise ValueError(f"line {line}: invalid peak area") from exc
        if not all(math.isfinite(v) and v >= 0 for v in (formate, reference)) or reference == 0:
            raise ValueError(f"line {line}: invalid target or reference area")
        ratio = formate / reference
        ratios.append(ratio)
        if row["known_formate_mM"].strip():
            if row["qc"] != "pass":
                raise ValueError(f"line {line}: calibration standard QC is not pass")
            concentration = float(row["known_formate_mM"])
            if not math.isfinite(concentration) or concentration < 0:
                raise ValueError(f"line {line}: invalid standard concentration")
            standards.append((ratio, concentration))
    if len(standards) < 3 or len({x for x, _ in standards}) < 3:
        raise ValueError("at least three distinct, reviewed standard ratios are required")
    xbar = sum(x for x, _ in standards) / len(standards)
    ybar = sum(y for _, y in standards) / len(standards)
    slope = sum((x - xbar) * (y - ybar) for x, y in standards) / sum((x - xbar) ** 2 for x, _ in standards)
    intercept = ybar - slope * xbar
    if slope <= 0:
        raise ValueError("nonpositive response slope; review integration and standards")
    low, high = min(x for x, _ in standards), max(x for x, _ in standards)
    for row, ratio in zip(rows, ratios, strict=True):
        row["peak_ratio"] = f"{ratio:.9g}"
        row["predicted_formate_mM"] = f"{slope * ratio + intercept:.9g}"
        row["calibration_status"] = "within_range" if low <= ratio <= high else "outside_range"
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*reader.fieldnames, "peak_ratio", "predicted_formate_mM", "calibration_status"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} ratios; fit concentration_mM={slope:.6g}*ratio+{intercept:.6g}")
    print("Predictions require original-spectrum QC and dilution correction before reporting")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.input, args.output)


if __name__ == "__main__":
    main()
