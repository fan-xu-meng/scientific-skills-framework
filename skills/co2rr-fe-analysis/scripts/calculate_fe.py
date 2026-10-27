#!/usr/bin/env python3
"""Calculate Faradaic efficiency from a reviewed product-amount CSV.

Input columns: product,moles,electrons,charge_C. All rows must share one charge
interval. Gas/liquid calibration and interval alignment happen upstream.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

FARADAY_C_PER_MOL = 96485.33212
REQUIRED = ("product", "moles", "electrons", "charge_C")


def positive_number(row: dict[str, str], field: str, line: int, *, allow_zero: bool = False) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"line {line}: invalid {field}") from exc
    if not math.isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        raise ValueError(f"line {line}: {field} must be finite and {'nonnegative' if allow_zero else 'positive'}")
    return value


def calculate(source: Path, target: Path) -> None:
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or any(name not in reader.fieldnames for name in REQUIRED):
            raise ValueError("input needs columns: " + ", ".join(REQUIRED))
        rows = list(reader)
    if not rows:
        raise ValueError("input contains no product rows")

    output: list[dict[str, str]] = []
    charges: set[float] = set()
    products: set[str] = set()
    for line, row in enumerate(rows, start=2):
        name = row["product"].strip()
        if not name:
            raise ValueError(f"line {line}: product name is empty")
        if name in products:
            raise ValueError(f"line {line}: duplicate product {name!r}")
        products.add(name)
        moles = positive_number(row, "moles", line, allow_zero=True)
        electrons = positive_number(row, "electrons", line)
        charge = positive_number(row, "charge_C", line)
        charges.add(charge)
        fraction = electrons * FARADAY_C_PER_MOL * moles / charge
        output.append({**row, "fe_fraction": f"{fraction:.9g}", "fe_percent": f"{100 * fraction:.9g}"})

    if len(charges) != 1:
        raise ValueError("all product rows must use the same charge_C interval magnitude")

    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*reader.fieldnames, "fe_fraction", "fe_percent"])
        writer.writeheader()
        writer.writerows(output)
    total_percent = sum(float(row["fe_percent"]) for row in output)
    print(f"Wrote {len(output)} products to {target}; summed FE={total_percent:.3f}%")
    if total_percent > 100:
        print("QC HOLD: summed FE exceeds 100%; check charge interval, calibration, and product overlap")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    calculate(args.input, args.output)


if __name__ == "__main__":
    main()
