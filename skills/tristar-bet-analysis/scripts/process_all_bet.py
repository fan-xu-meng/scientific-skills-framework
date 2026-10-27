#!/usr/bin/env python3
"""Extract TriStar II adsorption data and reported BET surface properties."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import re
import shutil
import sqlite3
import struct

import matplotlib.pyplot as plt
import numpy as np


SAMPLE_RECORD_MARKER = (
    b"\x00\xe0\x00\x00\x00\x00\x00\x00\xf8\x7f\x02\x00"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_directory_database(path: Path) -> tuple[list[dict[str, object]], str]:
    if not path.exists():
        return [], "not_present"
    try:
        connection = sqlite3.connect(path)
        rows = [
            {"file_name": row[0], "attribute_id": row[1], "value": row[2]}
            for row in connection.execute(
                "SELECT FileName, AttributeID, Value FROM Attributes ORDER BY FileName, AttributeID"
            )
        ]
        connection.close()
        return rows, "read_ok"
    except (sqlite3.DatabaseError, sqlite3.OperationalError) as exc:
        return [], f"read_failed:{type(exc).__name__}:{exc}"


def parse_smp_measurements(path: Path) -> list[dict[str, object]]:
    data = path.read_bytes()
    start = data.index(b"SUBSET303")
    end = data.find(b"SUBSET", start + len(b"SUBSET303"))
    if end < 0:
        end = len(data)
    subset = data[start:end]
    rows: list[dict[str, object]] = []
    for position in range(len(subset) - 40):
        elapsed = struct.unpack_from("<H", subset, position)[0]
        if subset[position + 2] != 0:
            continue
        selection_count = struct.unpack_from("<H", subset, position + 3)[0]
        if selection_count > 12:
            continue
        marker_position = position + 5 + 2 * selection_count
        if subset[marker_position : marker_position + len(SAMPLE_RECORD_MARKER)] != SAMPLE_RECORD_MARKER:
            continue
        values_position = marker_position + len(SAMPLE_RECORD_MARKER)
        absolute_pressure, relative_pressure, raw_collected_value = struct.unpack_from(
            "<3d", subset, values_position
        )
        if not (
            0 <= absolute_pressure < 2000
            and 0 <= relative_pressure <= 1.1
            and 0 <= raw_collected_value < 1000
        ):
            continue
        labels = subset[position + 5 : marker_position].decode("utf-16le", "ignore")
        rows.append(
            {
                "elapsed_seconds_raw": elapsed,
                "selection_labels": labels,
                "absolute_pressure_mmHg": absolute_pressure,
                "relative_pressure": relative_pressure,
                "raw_collected_value": raw_collected_value,
            }
        )
    if not rows:
        raise ValueError(f"No TriStar measurement records found in {path}")
    return rows


def utf16_strings(data: bytes) -> list[str]:
    return [
        match.group().decode("utf-16le")
        for match in re.finditer(rb"(?:(?:[\x20-\x7e])\x00){4,}", data)
    ]


def numeric_prefix(value: str) -> float:
    match = re.match(r"\s*([-+]?\d+(?:\.\d+)?)", value)
    if not match:
        raise ValueError(f"No numeric value found in {value!r}")
    return float(match.group(1))


def parse_rep_report(path: Path, smp_rows: list[dict[str, object]]) -> dict[str, object]:
    data = path.read_bytes()
    strings = utf16_strings(data)
    point_count = len(smp_rows)
    relative_pressure = np.array([float(row["relative_pressure"]) for row in smp_rows])
    absolute_pressure = np.array([float(row["absolute_pressure_mmHg"]) for row in smp_rows])

    pressure_bytes = struct.pack(f"<{point_count}d", *relative_pressure)
    relative_offset = data.find(pressure_bytes)
    if relative_offset < 0:
        raise ValueError("The REP file does not contain the SMP relative-pressure sequence")
    array_stride = point_count * 8 + 47
    absolute_offset = relative_offset + array_stride
    quantity_offset = absolute_offset + array_stride
    report_absolute = np.array(struct.unpack_from(f"<{point_count}d", data, absolute_offset))
    if not np.allclose(report_absolute, absolute_pressure, rtol=0, atol=1e-10):
        raise ValueError("REP absolute-pressure array failed the SMP cross-check")
    quantity = np.array(struct.unpack_from(f"<{point_count}d", data, quantity_offset))
    if not np.all(np.isfinite(quantity)) or np.any(quantity <= 0):
        raise ValueError("REP quantity-adsorbed array is invalid")

    times = [value for value in strings if re.fullmatch(r"\d{2}:\d{2}", value)]
    if len(times) < point_count:
        times = [""] * point_count
    else:
        times = times[:point_count]

    summary_values = [numeric_prefix(value) for value in strings[13:20]]
    volume_values = [numeric_prefix(value) for value in strings[32:36]]
    pore_size_values = [numeric_prefix(strings[index]) for index in (39, 40, 41, 43)]
    return {
        "relative_pressure": relative_pressure,
        "absolute_pressure_mmHg": absolute_pressure,
        "quantity_adsorbed_cm3_stp_g": quantity,
        "elapsed_time": times,
        "surface_areas": {
            "Single-point surface area": summary_values[0],
            "BET multipoint surface area": summary_values[1],
            "Langmuir surface area": summary_values[2],
            "t-Plot micropore area": summary_values[3],
            "t-Plot external surface area": summary_values[4],
            "BJH adsorption cumulative pore area": summary_values[5],
            "BJH desorption cumulative pore area": summary_values[6],
        },
        "pore_properties": {
            "Single-point adsorption total pore volume": volume_values[0],
            "t-Plot micropore volume": volume_values[1],
            "BJH adsorption cumulative pore volume": volume_values[2],
            "BJH desorption cumulative pore volume": volume_values[3],
            "Adsorption average pore width": pore_size_values[0],
            "BJH adsorption average pore radius": pore_size_values[1],
            "BJH desorption average pore radius": pore_size_values[2],
            "Average particle size": pore_size_values[3],
        },
    }


def plot_isotherm(report: dict[str, object], path: Path) -> None:
    pressure = np.asarray(report["relative_pressure"], dtype=float)
    quantity = np.asarray(report["quantity_adsorbed_cm3_stp_g"], dtype=float)
    maximum_index = int(np.argmax(pressure))
    fig, axis = plt.subplots(figsize=(6.6, 4.8), dpi=180)
    axis.plot(
        pressure[: maximum_index + 1],
        quantity[: maximum_index + 1],
        "o-",
        color="#167D9A",
        markersize=4,
        linewidth=1.5,
        label="Adsorption",
    )
    axis.plot(
        pressure[maximum_index + 1 :],
        quantity[maximum_index + 1 :],
        "s-",
        color="#D66A3A",
        markersize=3.8,
        linewidth=1.5,
        label="Desorption",
    )
    axis.set_xlim(0, 1.02)
    axis.set_ylim(bottom=0)
    axis.set_xlabel(r"Relative pressure, $P/P_0$")
    axis.set_ylabel(r"Quantity adsorbed (cm$^3$ STP g$^{-1}$)")
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(direction="out")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_bet_fit(report: dict[str, object], path: Path) -> dict[str, float]:
    pressure = np.asarray(report["relative_pressure"], dtype=float)
    quantity = np.asarray(report["quantity_adsorbed_cm3_stp_g"], dtype=float)
    adsorption_end = int(np.argmax(pressure)) + 1
    selected = (pressure[:adsorption_end] >= 0.065) & (pressure[:adsorption_end] <= 0.200)
    x = pressure[:adsorption_end][selected]
    y = x / (quantity[:adsorption_end][selected] * (1 - x))
    slope, intercept = np.polyfit(x, y, 1)
    fitted = slope * x + intercept
    r_squared = 1 - float(np.sum((y - fitted) ** 2) / np.sum((y - np.mean(y)) ** 2))
    monolayer_quantity = 1 / (slope + intercept)
    bet_constant = slope / intercept + 1
    nitrogen_cross_section_nm2 = 0.162
    molar_volume_cm3_stp = 22414.0
    avogadro_constant = 6.02214076e23
    calculated_surface_area = (
        monolayer_quantity
        / molar_volume_cm3_stp
        * avogadro_constant
        * nitrogen_cross_section_nm2
        * 1e-18
    )

    fig, axis = plt.subplots(figsize=(5.8, 4.4), dpi=180)
    axis.scatter(x, y, color="#167D9A", s=28, label="Selected data")
    axis.plot(x, fitted, color="#333333", linewidth=1.4, label="Linear fit")
    axis.set_xlabel(r"Relative pressure, $P/P_0$")
    axis.set_ylabel(r"$P/[Q(P_0-P)]$ (g cm$^{-3}$ STP)")
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(direction="out")
    axis.text(
        0.04,
        0.96,
        f"$R^2$ = {r_squared:.6f}\nPoints = {len(x)}",
        transform=axis.transAxes,
        va="top",
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return {
        "selected_point_count": float(len(x)),
        "selected_pressure_min": float(x.min()),
        "selected_pressure_max": float(x.max()),
        "slope": float(slope),
        "intercept": float(intercept),
        "correlation_coefficient": float(math.sqrt(r_squared)),
        "r_squared": r_squared,
        "monolayer_quantity_cm3_stp_g": float(monolayer_quantity),
        "bet_constant": float(bet_constant),
        "calculated_bet_surface_area_m2_g": float(calculated_surface_area),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    data_root = args.data_root.resolve()
    out_dir = args.out_dir.resolve()
    if not data_root.is_dir():
        raise NotADirectoryError(f"BET input directory not found: {data_root}")
    if out_dir.exists():
        raise FileExistsError(f"Output directory already exists: {out_dir}")
    out_dir.mkdir(parents=True)

    files = sorted(path for path in data_root.rglob("*") if path.is_file() and path.name != ".DS_Store")
    inventory = [
        {
            "file": str(path.relative_to(data_root)),
            "extension": path.suffix.lower(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
            "first_16_bytes_hex": path.read_bytes()[:16].hex(" "),
            "classification": (
                "TriStar_sample_binary"
                if path.suffix.lower() == ".smp"
                else "TriStar_report_binary"
                if path.suffix.lower() == ".rep"
                else "Micromeritics_directory_database"
                if path.suffix.lower() == ".db3"
                else "unrecognized"
            ),
        }
        for path in files
    ]
    write_csv(
        out_dir / "bet_file_inventory.csv",
        inventory,
        ["file", "extension", "size_bytes", "sha256", "first_16_bytes_hex", "classification"],
    )

    attribute_rows, db_status = read_directory_database(data_root / "MICDIR.DB3")
    write_csv(
        out_dir / "micromeritics_directory_attributes.csv",
        attribute_rows,
        ["file_name", "attribute_id", "value"],
    )

    sample_path = next(iter(sorted(data_root.glob("*.SMP"))), None)
    report_path = next(iter(sorted(data_root.glob("*.rep"))), None)
    if sample_path is None or report_path is None:
        raise FileNotFoundError("A matching TriStar SMP and REP pair is required")
    smp_rows = parse_smp_measurements(sample_path)
    report = parse_rep_report(report_path, smp_rows)
    pressure = np.asarray(report["relative_pressure"], dtype=float)
    maximum_index = int(np.argmax(pressure))
    isotherm_rows = []
    for index in range(len(smp_rows)):
        isotherm_rows.append(
            {
                "point_index": index + 1,
                "branch": "adsorption" if index <= maximum_index else "desorption",
                "elapsed_time": report["elapsed_time"][index],
                "absolute_pressure_mmHg": f'{report["absolute_pressure_mmHg"][index]:.9f}',
                "relative_pressure": f'{pressure[index]:.9f}',
                "quantity_adsorbed_cm3_stp_g": f'{report["quantity_adsorbed_cm3_stp_g"][index]:.9f}',
            }
        )
    write_csv(
        out_dir / "bet_isotherm.csv",
        isotherm_rows,
        [
            "point_index",
            "branch",
            "elapsed_time",
            "absolute_pressure_mmHg",
            "relative_pressure",
            "quantity_adsorbed_cm3_stp_g",
        ],
    )

    area_rows = []
    for method, value in report["surface_areas"].items():
        details = ""
        if method.startswith("Single-point"):
            details = "P/P0 = 0.199864355"
        elif method.startswith("BJH"):
            details = "Pore radius range = 8.500-1500.000 angstrom"
        area_rows.append(
            {
                "method": method,
                "surface_area_m2_g": f"{value:.4f}",
                "source": "TriStar II 3020 report",
                "details": details,
            }
        )
    write_csv(
        out_dir / "bet_surface_area_summary.csv",
        area_rows,
        ["method", "surface_area_m2_g", "source", "details"],
    )

    pore_rows = []
    for property_name, value in report["pore_properties"].items():
        if "volume" in property_name.lower():
            unit = "cm3_g"
        elif "particle" in property_name.lower():
            unit = "angstrom"
        else:
            unit = "angstrom"
        pore_rows.append(
            {
                "property": property_name,
                "value": f"{value:.6f}",
                "unit": unit,
                "source": "TriStar II 3020 report",
            }
        )
    write_csv(
        out_dir / "bet_pore_properties_summary.csv",
        pore_rows,
        ["property", "value", "unit", "source"],
    )

    plot_isotherm(report, out_dir / "bet_adsorption_desorption_isotherm.png")
    fit = plot_bet_fit(report, out_dir / "bet_linear_fit.png")
    write_csv(
        out_dir / "bet_linear_fit_summary.csv",
        [{key: f"{value:.9g}" for key, value in fit.items()}],
        list(fit),
    )

    summary = [
        "BET processing summary",
        "======================",
        f"Input sample: {sample_path}",
        f"Input report: {report_path}",
        f"Micromeritics directory database status: {db_status}",
        f"Isotherm points: {len(isotherm_rows)}",
        f"Adsorption points: {maximum_index + 1}",
        f"Desorption points: {len(isotherm_rows) - maximum_index - 1}",
        f"Maximum quantity adsorbed: {max(report['quantity_adsorbed_cm3_stp_g']):.6f} cm3 STP/g",
        f"Reported BET surface area: {report['surface_areas']['BET multipoint surface area']:.4f} m2/g",
        "",
        "Data provenance:",
        "- Pressure arrays were parsed from the SMP measurement subset and cross-checked against the REP arrays.",
        "- Quantity adsorbed and reported surface/pore properties were parsed from the matching REP file.",
        "- The adsorption/desorption split occurs immediately after the maximum relative-pressure point.",
        "- The source SMP and REP files were not modified.",
        "",
        "Generated numeric products:",
        "- bet_isotherm.csv",
        "- bet_surface_area_summary.csv",
        "- bet_pore_properties_summary.csv",
        "- bet_linear_fit_summary.csv",
        "- bet_adsorption_desorption_isotherm.png",
        "- bet_linear_fit.png",
        "",
    ]
    (out_dir / "bet_processing_summary.txt").write_text("\n".join(summary), encoding="utf-8")
    (out_dir / "processing_config.yaml").write_text(
        "workflow: tristar_ii_bet_binary_report_extraction\n"
        f"data_root: {data_root}\n"
        f"output_root: {out_dir}\n"
        "instrument_family: TriStar II 3020\n"
        "isotherm_source: matching_smp_and_rep_arrays\n"
        "surface_area_source: tristar_report\n"
        "source_files_modified: false\n",
        encoding="utf-8",
    )
    shutil.copy2(__file__, out_dir / Path(__file__).name)
    print(f"Isotherm points: {len(isotherm_rows)}")
    print(f"Adsorption/desorption: {maximum_index + 1}/{len(isotherm_rows) - maximum_index - 1}")
    print(f"BET surface area: {report['surface_areas']['BET multipoint surface area']:.4f} m2/g")
    print(f"Output: {out_dir}")


if __name__ == "__main__":
    main()
