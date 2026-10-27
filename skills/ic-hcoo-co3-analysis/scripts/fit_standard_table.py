#!/usr/bin/env python3
"""Fit linear IC area-vs-concentration standards from a reviewed CSV."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def fit(points: list[tuple[float, float]]) -> tuple[float, float, float, float]:
    if len(points) < 3 or len({x for x, _ in points}) < 3:
        raise ValueError("each analyte needs at least three distinct concentration levels")
    n = len(points)
    xbar = sum(x for x, _ in points) / n
    ybar = sum(y for _, y in points) / n
    sxx = sum((x - xbar) ** 2 for x, _ in points)
    slope = sum((x - xbar) * (y - ybar) for x, y in points) / sxx
    intercept = ybar - slope * xbar
    residuals = [y - (slope * x + intercept) for x, y in points]
    sst = sum((y - ybar) ** 2 for _, y in points)
    r2 = 1 - sum(e * e for e in residuals) / sst if sst > 0 else math.nan
    return slope, intercept, r2, max(abs(e) for e in residuals)


def run(source: Path, target: Path) -> None:
    groups: dict[str, list[tuple[float, float]]] = defaultdict(list)
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not {"analyte", "concentration_mM", "area_uS_min"}.issubset(reader.fieldnames or ()):
            raise ValueError("input requires analyte,concentration_mM,area_uS_min")
        for line, row in enumerate(reader, start=2):
            analyte = row["analyte"].strip()
            try:
                concentration = float(row["concentration_mM"])
                area = float(row["area_uS_min"])
            except ValueError as exc:
                raise ValueError(f"line {line}: invalid numeric value") from exc
            if not analyte or not all(math.isfinite(v) and v >= 0 for v in (concentration, area)):
                raise ValueError(f"line {line}: invalid analyte, concentration, or area")
            groups[analyte].append((concentration, area))
    if not groups:
        raise ValueError("no standards found")
    output = []
    for analyte, points in sorted(groups.items()):
        slope, intercept, r2, max_residual = fit(points)
        output.append({"analyte": analyte, "n_standards": len(points),
                       "minimum_mM": min(x for x, _ in points),
                       "maximum_mM": max(x for x, _ in points),
                       "slope_area_per_mM": slope, "intercept_area": intercept,
                       "r2": r2, "max_abs_residual_area": max_residual,
                       "qc": "review" if slope <= 0 or not math.isfinite(r2) else "pending_chromatogram_review"})
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    print(f"Wrote {len(output)} analyte fits; inspect residuals and raw chromatograms")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.input, args.output)


if __name__ == "__main__":
    main()
