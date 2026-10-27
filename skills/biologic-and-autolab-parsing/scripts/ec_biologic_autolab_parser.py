"""
ec_biologic_autolab_parser.py
=============================
Reusable module: parse BioLogic .mpr/.mpt and Autolab NOVA .nox files
into a canonical (time_s, current_A, potential_V, charge_C) DataFrame,
then generate publication-style CA/CP traces.

Can also be run as a standalone script on a folder of CA files:

    python ec_biologic_autolab_parser.py \\
        --data-dir  "path/to/CA_folder" \\
        --output-dir "output/EC"
"""

from __future__ import annotations

import argparse
import hashlib
import math
import re
import struct
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ── Canonical column mapping for BioLogic ─────────────────────────────────────
_COL_MAP = {
    "time/s":       "time_s",
    "Ewe/V":        "potential_V",
    "<Ewe>/V":      "potential_V",
    "<I>/mA":       "_I_mA",
    "I/mA":         "_I_mA",
    "(Q-Qo)/C":     "charge_C",      # CP: already Coulombs
    "Q-Q0/mA.h":   "_Q_mAh",
    "(Q-Qo)/mA.h": "_Q_mAh",
    "control/V":    "ctrl_V",
    "control/mA":   "ctrl_mA",
}

# ── Plot style ─────────────────────────────────────────────────────────────────
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
          "#9467bd", "#8c564b", "#e377c2", "#17becf"]
CHARGE_COLOR = "#555555"

_SKIP_RE = re.compile(r"(?i)(cv|peis|lsv|blank|\bCV\b)")
_VOLT_RE = re.compile(r"([0-9]+\.?[0-9]*)[\s_]?[Vv](?:\s|_|$|[\s,\.])")
_NUM_AFTER_FOR = re.compile(r"(?:for|at)\s+([0-9]+\.?[0-9]*)\b", re.IGNORECASE)


# ═══════════════════════════════════════════════════════════════════════════════
# BioLogic .mpr / .mpt parser
# ═══════════════════════════════════════════════════════════════════════════════

def parse_mpr(path: str | Path) -> pd.DataFrame:
    """
    Parse a BioLogic .mpr (or .mpt) file.
    Returns DataFrame with canonical columns: time_s, current_A, potential_V, charge_C.

    Library priority: galvani → eclabfiles.
    """
    path = Path(path)
    try:
        import galvani
        raw = pd.DataFrame(galvani.MPRfile(str(path)).data)
    except Exception:
        try:
            import eclabfiles as ecf
            raw = ecf.to_df(str(path))
        except Exception as exc:
            raise RuntimeError(f"Cannot parse {path.name}: {exc}") from exc

    df = raw.rename(columns={k: v for k, v in _COL_MAP.items() if k in raw.columns})

    # Make time relative (start at 0)
    if "time_s" in df.columns:
        df["time_s"] = df["time_s"] - df["time_s"].iloc[0]

    # Current → Amperes
    if "current_A" not in df.columns and "_I_mA" in df.columns:
        df["current_A"] = df["_I_mA"] / 1000.0

    # Charge → Coulombs
    if "charge_C" not in df.columns:
        if "_Q_mAh" in df.columns:
            df["charge_C"] = df["_Q_mAh"] * 3.6          # mA·h → C
        elif "current_A" in df.columns and "time_s" in df.columns:
            dt = np.diff(df["time_s"].values, prepend=df["time_s"].iloc[0])
            df["charge_C"] = np.cumsum(df["current_A"].values * dt)

    keep = ["time_s", "current_A", "potential_V", "charge_C"]
    return df[[c for c in keep if c in df.columns]].copy()


# ═══════════════════════════════════════════════════════════════════════════════
# Autolab NOVA .nox binary parser
# ═══════════════════════════════════════════════════════════════════════════════

def _find_corr_time(data: bytes) -> tuple[int, int, float]:
    """
    Locate the 'Corrected time (s)' float64 array inside an Autolab NOVA .nox
    binary blob.

    The array starts at exactly 0.0 and increases monotonically with a
    consistent step (<10 s).  A run of ≥200 consecutive valid values is
    required before the full length is counted.

    Returns
    -------
    (byte_offset, n_points, step_s)
    """
    end = len(data) - 8 * 100
    for pos in range(0, end):
        v0, = struct.unpack_from("<d", data, pos)
        if v0 != 0.0:
            continue
        try:
            v1, = struct.unpack_from("<d", data, pos + 8)
        except struct.error:
            continue
        if not (0.0 < v1 < 10.0):
            continue
        step = v1
        prev, ok = 0.0, True
        for j in range(1, min(200, (len(data) - pos) // 8)):
            try:
                v, = struct.unpack_from("<d", data, pos + j * 8)
            except struct.error:
                ok = False; break
            if not (v > prev and abs((v - prev) - step) < step * 0.01):
                ok = False; break
            prev = v
        if not ok:
            continue
        # Count full length
        count = 200
        prev = v  # last confirmed value
        for j in range(200, (len(data) - pos) // 8):
            try:
                v, = struct.unpack_from("<d", data, pos + j * 8)
            except struct.error:
                break
            if not np.isfinite(v) or v <= prev or abs((v - prev) - step) > step * 0.05:
                count = j; break
            prev = v
        else:
            count = (len(data) - pos) // 8
        return pos, count, step
    raise ValueError("Corrected time array not found in .nox file")


def _find_current(data: bytes, ct_offset: int, n: int,
                  ct_step: float, val_max: float = 1.0) -> np.ndarray:
    """
    Locate the current float64 array of length n in the .nox binary blob.

    Searches all 8 byte-alignment offsets in the region before the Corrected
    time array.  Picks the first run of n values all in [-val_max, val_max] A
    with non-trivial standard deviation (>1e-8).

    Parameters
    ----------
    val_max : float
        Upper bound on |current| in Amperes.  Default 1 A covers most
        experiments.  Raise to ~10 for high-current (> 1 A) runs.
    """
    scan_start = max(0, ct_offset - 4 * n * 8)
    scan_end   = ct_offset + n * 8
    scan_data  = data[scan_start:scan_end]

    for align in range(8):
        chunk = scan_data[align:]
        pad   = (-len(chunk)) % 8
        raw   = np.frombuffer(chunk + b"\x00" * pad, dtype="<f8")
        i = 0
        while i <= len(raw) - n:
            window = raw[i : i + n]
            wmin, wmax = float(window.min()), float(window.max())
            if (
                -val_max <= wmin and wmax <= val_max
                and np.all(np.isfinite(window))
                and float(window.std()) > 1e-8
            ):
                return window.copy()
            i += 1

    raise ValueError(
        f"Current array not found in .nox file (n={n}, val_max={val_max} A). "
        "Try increasing val_max for high-current runs."
    )


def parse_nox(path: str | Path, current_max_A: float = 1.0) -> pd.DataFrame:
    """
    Parse an Autolab NOVA .nox binary file.

    Returns DataFrame with canonical columns: time_s, current_A, charge_C.
    Note: .nox CA files do not store potential as a variable column;
    the applied potential is embedded in the filename / folder name.

    Parameters
    ----------
    path : str | Path
    current_max_A : float
        Upper bound on expected |current| in Amperes (default 1.0 A = 1000 mA).
        Raise to 10 or higher for high-current experiments.
    """
    data = Path(path).read_bytes()
    ct_offset, n, ct_step = _find_corr_time(data)
    corr_time = np.frombuffer(data, dtype="<f8", count=n, offset=ct_offset)
    current   = _find_current(data, ct_offset, n, ct_step, val_max=current_max_A)
    dt        = np.diff(corr_time, prepend=corr_time[0] - ct_step)
    charge    = np.cumsum(current * dt)
    return pd.DataFrame({
        "time_s":   corr_time,
        "current_A": current,
        "charge_C": charge,
    })


# ═══════════════════════════════════════════════════════════════════════════════
# Unified reader
# ═══════════════════════════════════════════════════════════════════════════════

def read_ec(path: str | Path, **kwargs) -> pd.DataFrame:
    """
    Auto-detect file format and return canonical DataFrame.

    Supported: .nox (Autolab NOVA), .mpr / .mpt (BioLogic EC-Lab).
    kwargs are forwarded to the relevant parser (e.g. current_max_A for .nox).
    """
    path = Path(path)
    if path.suffix.lower() == ".nox":
        return parse_nox(path, **kwargs)
    elif path.suffix.lower() in (".mpr", ".mpt"):
        return parse_mpr(path)
    else:
        raise ValueError(f"Unsupported extension: {path.suffix}")


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def extract_potential(stem: str, folder_name: str = "") -> str:
    """Extract applied potential label from filename or folder name."""
    for text in [stem, folder_name]:
        m = _VOLT_RE.search(text)
        if m:
            return m.group(1) + " V"
    for text in [stem, folder_name]:
        m = _NUM_AFTER_FOR.search(text)
        if m:
            return m.group(1) + " V"
    return "? V"


def sha256_hex(path: Path, n: int = 16) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()[:n]


def _subplot_grid(n: int) -> tuple[int, int]:
    r = max(1, int(math.sqrt(n)))
    c = r if r * r == n else r + 1
    return r, c


# ═══════════════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════════════

def _setup_ax(ax: plt.Axes) -> None:
    ax.grid(False)
    ax.tick_params(direction="in")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_ca_traces(entries: list[dict], out_path: str | Path, title: str) -> None:
    """
    Plot current (mA) vs time (min) with cumulative charge (C) on right y-axis.

    Parameters
    ----------
    entries : list of dicts with keys:
        'label' : str
        'df'    : DataFrame with time_s, current_A, charge_C
        'color' : str (optional, default from COLORS cycle)
    """
    fig, ax = plt.subplots(figsize=(7, 4))
    ax2 = ax.twinx()
    for idx, entry in enumerate(entries):
        df    = entry["df"]
        color = entry.get("color", COLORS[idx % len(COLORS)])
        t     = df["time_s"] / 60.0
        ax.plot(t, df["current_A"] * 1000.0, color=color,
                lw=1.2, alpha=0.85, label=entry["label"])
        ax2.plot(t, df["charge_C"], color=color, lw=1.0, alpha=0.5, ls="--")
    ax.set_xlabel("Time (min)", fontsize=10)
    ax.set_ylabel("Current (mA)", fontsize=10)
    ax2.set_ylabel("Cumulative charge (C)", fontsize=10, color=CHARGE_COLOR)
    ax2.tick_params(axis="y", labelcolor=CHARGE_COLOR)
    ax2.grid(False)
    _setup_ax(ax)
    ax.spines["right"].set_visible(True)
    leg = ax.legend(fontsize=8, loc="best")
    leg.get_frame().set_visible(False)
    ax.set_title(title, fontsize=10, pad=6)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_ca_charge(entries: list[dict], out_path: str | Path, title: str) -> None:
    """
    Plot cumulative charge (C) vs time (min) for all entries.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    for idx, entry in enumerate(entries):
        df    = entry["df"]
        color = entry.get("color", COLORS[idx % len(COLORS)])
        ax.plot(df["time_s"] / 60.0, df["charge_C"],
                color=color, lw=1.4, alpha=0.85, label=entry["label"])
    ax.set_xlabel("Time (min)", fontsize=10)
    ax.set_ylabel("Cumulative charge (C)", fontsize=10)
    ax.set_title(title, fontsize=10, pad=6)
    _setup_ax(ax)
    leg = ax.legend(fontsize=8, loc="best")
    leg.get_frame().set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# Per-folder processor
# ═══════════════════════════════════════════════════════════════════════════════

def process_folder(data_folder: Path, out_folder: Path) -> list[dict]:
    """
    Find all CA .nox / .mpr files in data_folder (skip CV/PEIS/LSV/blank),
    parse, plot ca_traces.png + ca_charge.png, and save ca_summary.csv.

    Returns a list of per-file summary dicts.
    """
    out_folder.mkdir(parents=True, exist_ok=True)
    folder_label = data_folder.name

    ca_files: list[Path] = []
    for suffix in (".nox", ".mpr", ".mpt"):
        for p in sorted(data_folder.glob(f"*CA*{suffix}")):
            if not _SKIP_RE.search(p.stem):
                ca_files.append(p)

    if not ca_files:
        print(f"  [skip] no CA files in {folder_label}")
        return []

    entries: list[dict] = []
    rows: list[dict] = []

    for idx, p in enumerate(ca_files):
        color = COLORS[idx % len(COLORS)]
        pot   = extract_potential(p.stem, folder_label)
        try:
            df = read_ec(p)
        except Exception as exc:
            print(f"  WARNING {p.name}: {exc}")
            continue

        # Truncate long label for legend
        raw_label = p.stem.split("_02_CA")[0].split(", CA")[0]
        label = f"{raw_label[-40:]} ({pot})"

        entries.append({"label": label, "df": df, "color": color})
        rows.append({
            "session":          folder_label,
            "file":             p.name,
            "source":           p.suffix.lower().lstrip("."),
            "potential_label":  pot,
            "duration_min":     round(float(df["time_s"].iloc[-1]) / 60, 2),
            "n_points":         len(df),
            "mean_current_mA":  round(float(df["current_A"].mean()) * 1000, 4),
            "total_charge_C":   round(float(df["charge_C"].iloc[-1]), 4),
            "sha256":           sha256_hex(p),
        })
        print(f"    {p.name}: {rows[-1]['duration_min']:.1f} min, "
              f"I̅={rows[-1]['mean_current_mA']:.3f} mA, "
              f"Q={rows[-1]['total_charge_C']:.3f} C")

    if not entries:
        return []

    plot_ca_traces(entries,
                   out_folder / "ca_traces.png",
                   f"{folder_label}\nCurrent & charge vs time")
    plot_ca_charge(entries,
                   out_folder / "ca_charge.png",
                   f"{folder_label}\nCumulative charge vs time")
    pd.DataFrame(rows).to_csv(out_folder / "ca_summary.csv", index=False)
    print(f"    → {len(rows)} files → ca_traces.png, ca_charge.png, ca_summary.csv")
    return rows


# ═══════════════════════════════════════════════════════════════════════════════
# CLI entry point
# ═══════════════════════════════════════════════════════════════════════════════

def _cli() -> None:
    parser = argparse.ArgumentParser(
        description="Parse .nox/.mpr CA files and generate traces + summary CSV"
    )
    parser.add_argument("--data-dir",    required=True,
                        help="Root folder containing per-session subfolders")
    parser.add_argument("--output-dir",  required=True,
                        help="Root output folder (mirrored structure)")
    args = parser.parse_args()

    data_root = Path(args.data_dir)
    out_root  = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    for data_folder in sorted(data_root.iterdir()):
        if not data_folder.is_dir():
            continue
        print(f"\n{'─'*60}\n  Session: {data_folder.name}\n{'─'*60}")
        rows = process_folder(data_folder, out_root / data_folder.name)
        all_rows.extend(rows)

    if all_rows:
        out_csv = out_root / "ca_all_sessions.csv"
        pd.DataFrame(all_rows).to_csv(out_csv, index=False)
        print(f"\n✓ Saved {len(all_rows)} rows → {out_csv}")


if __name__ == "__main__":
    _cli()
