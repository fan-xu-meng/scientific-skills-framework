#!/usr/bin/env python3
"""Audit basic numeric integrity of an exported XPS component table."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def run(source: Path) -> int:
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"region", "component", "center_eV", "fwhm_eV", "area_counts_eV", "qc"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("input requires " + ",".join(sorted(required)))
        rows = list(reader)
    if not rows:
        raise ValueError("component table is empty")
    holds: list[str] = []
    seen: set[tuple[str, str]] = set()
    for line, row in enumerate(rows, start=2):
        key = (row["region"].strip(), row["component"].strip())
        if not all(key) or key in seen:
            holds.append(f"line {line}: missing or duplicated region/component")
        seen.add(key)
        try:
            center = float(row["center_eV"])
            width = float(row["fwhm_eV"])
            area = float(row["area_counts_eV"])
        except ValueError:
            holds.append(f"line {line}: non-numeric component value")
            continue
        if not all(math.isfinite(v) for v in (center, width, area)) or width <= 0 or area < 0:
            holds.append(f"line {line}: non-finite center or invalid width/area")
        if row["qc"].strip().lower() not in {"pass", "review", "pending_assignment", "rejected"}:
            holds.append(f"line {line}: unrecognized QC status")
    print(f"Rows: {len(rows)}; numeric/identity holds: {len(holds)}")
    for hold in holds:
        print("QC HOLD:", hold)
    print("No chemical-state, background, charge-reference, or spin-orbit validation is implied.")
    return 1 if holds else 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    args = parser.parse_args()
    raise SystemExit(run(args.input))


if __name__ == "__main__":
    main()
