"""Agilent/Varian Cary UV-Vis CSV adaptive calibration.

Same three-pass design as nmr_adaptive_calibrate.py / ic_adaptive.py:
  Pass 1 — parse all spectra from the Cary export CSV, blank-subtract,
            auto-detect peak wavelength from calibration standards,
            read absorbance at that wavelength for every spectrum.
  Pass 2 — linear regression: absorbance → concentration from
            known-concentration standard columns.
  Pass 3 — apply regression to unknown sample columns.

CSV format (Agilent Cary 5000 multi-spectrum export):
  Row 0 : sample labels in EVEN columns (0, 2, 4, …); odd columns empty.
          Concentration encoded as "0.02mmol", "0.04mmol", etc.
          Blank/baseline encoded as "Baseline 100%T" (or similar).
  Row 1 : column sub-headers — "Wavelength (nm)", "Abs" repeated.
  Rows 2+: interleaved data — wavelength_i, abs_i, wavelength_j, abs_j, …
  Footer : non-numeric metadata lines after the last data row (ignored).

Standard detection  : column label contains e.g. "0.02mmol" or "0.04 mM"
Blank detection     : column label starts with "Base" / "Blank" / "100%"

Peak wavelength auto-detection:
  - Blank-subtract each standard spectrum.
  - Average across all standards at each wavelength.
  - Peak = wavelength of maximum average blank-subtracted absorbance
    within WL_SEARCH_RANGE.
  - Can be overridden with --cal-wl.

Usage
-----
uv run python uv_cary_csv_adaptive.py \\
    --csv "<spectra-export.csv>" \\
    --out-dir output/uv-vis \\
    --label "UV-Vis calibration"

# Fix calibration wavelength manually:
uv run python uv_cary_csv_adaptive.py --csv ... --cal-wl 660 --out-dir ...
"""

from __future__ import annotations

import argparse
import math
import re
import warnings
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

# ─────────────────────────────────────────────────────────────────────────────
#  Constants
# ─────────────────────────────────────────────────────────────────────────────

# Wavelength range searched for peak auto-detection (nm)
WL_SEARCH_RANGE: tuple[float, float] = (500.0, 800.0)

# Common wavelength grid for interpolation (descending, high → low nm)
_WL_GRID_STEP = 1.0  # nm

# Colour palette — same tableau cycle as GC/HPLC adaptive scripts.
# Standards use the blue anchor; samples cycle through the full tableau set.
_TABLEAU = [
    "#1f77b4",  # blue   (H₂ in GC)
    "#d62728",  # red    (CO  in GC)
    "#2ca02c",  # green  (CH₄ in GC)
    "#ff7f0e",  # orange (CO₂ in GC)
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#17becf",  # cyan
    "#bcbd22",  # olive
]
_STD_COLORS  = _TABLEAU[:4]   # blue, red, green, orange for up to 4 standards
_SAMP_COLORS = _TABLEAU        # full palette for samples


# ─────────────────────────────────────────────────────────────────────────────
#  CSV parser
# ─────────────────────────────────────────────────────────────────────────────

def _parse_conc_from_label(label: str) -> float | None:
    """Return concentration in mmol from a label like '0.02mmol' or '0.04 mM'."""
    m = re.search(r"([\d.]+)\s*mmol", label, re.IGNORECASE)
    if m:
        return float(m.group(1))
    m = re.search(r"([\d.]+)\s*mM", label, re.IGNORECASE)
    if m:
        return float(m.group(1))
    return None


def _is_blank(label: str) -> bool:
    return bool(re.search(r"[Bb]ase|[Bb]lank|100\s*%", label))


def parse_cary_csv(csv_path: Path) -> list[dict]:
    """Parse a Cary multi-spectrum CSV export.

    Returns a list of dicts, one per spectrum column pair:
        label     : str   — column header from row 0
        wl        : ndarray — wavelengths (nm), ascending
        abs       : ndarray — absorbance values
        conc_mmol : float | None — concentration if parseable from label
        is_blank  : bool
        is_std    : bool — True when conc_mmol is not None
    """
    raw = csv_path.read_text(encoding="utf-8", errors="replace")
    lines = raw.splitlines()

    # --- Row 0: sample labels (even column positions) ---
    header_row = next((i for i, l in enumerate(lines)
                       if "Wavelength" in l or "wavelength" in l), 1) - 1
    label_row_idx = max(header_row, 0)
    label_parts = lines[label_row_idx].split(",")

    labels: list[str] = []
    n_cols = len(label_parts)
    for col in range(0, n_cols, 2):
        lbl = label_parts[col].strip()
        if lbl:
            labels.append(lbl)

    # --- Data rows: stop at first non-numeric line ---
    data_start = header_row + 2          # skip label row + sub-header row
    data_rows: list[list[float | None]] = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith(","):
            break
        parts = stripped.split(",")
        try:
            float(parts[0])
        except (ValueError, IndexError):
            break
        row_vals: list[float | None] = []
        for p in parts:
            try:
                row_vals.append(float(p))
            except ValueError:
                row_vals.append(None)
        data_rows.append(row_vals)

    if not data_rows:
        raise ValueError(f"No numeric data found in {csv_path.name}")

    arr = np.array(
        [[v if v is not None else np.nan for v in r] for r in data_rows],
        dtype=float,
    )

    # --- Build spectrum list ---
    spectra: list[dict] = []
    for i, label in enumerate(labels):
        wl_col = i * 2
        ab_col = i * 2 + 1
        if ab_col >= arr.shape[1]:
            break
        wl = arr[:, wl_col]
        ab = arr[:, ab_col]
        # Remove NaN rows and sort ascending wavelength
        valid = np.isfinite(wl) & np.isfinite(ab)
        wl, ab = wl[valid], ab[valid]
        if len(wl) < 5:
            continue
        order = np.argsort(wl)
        wl, ab = wl[order], ab[order]
        conc = _parse_conc_from_label(label)
        spectra.append({
            "label":     label,
            "wl":        wl,
            "abs":       ab,
            "conc_mmol": conc,
            "is_blank":  _is_blank(label),
            "is_std":    conc is not None,
        })

    return spectra


# ─────────────────────────────────────────────────────────────────────────────
#  Blank subtraction + peak detection
# ─────────────────────────────────────────────────────────────────────────────

def _common_wl_grid(spectra: list[dict]) -> np.ndarray:
    lo = max(WL_SEARCH_RANGE[0], max(float(s["wl"].min()) for s in spectra))
    hi = min(WL_SEARCH_RANGE[1], min(float(s["wl"].max()) for s in spectra))
    return np.arange(lo, hi + _WL_GRID_STEP, _WL_GRID_STEP)


def blank_subtract(spectra: list[dict]) -> tuple[dict | None, list[dict]]:
    """Return (blank_spectrum, [blank-subtracted copies of all non-blank spectra])."""
    blanks = [s for s in spectra if s["is_blank"]]
    blank = blanks[0] if blanks else None

    result: list[dict] = []
    for s in spectra:
        if s["is_blank"]:
            continue
        ab_corr = s["abs"].copy()
        if blank is not None:
            ab_bl = np.interp(s["wl"], blank["wl"], blank["abs"],
                              left=np.nan, right=np.nan)
            ab_corr = ab_corr - ab_bl
        result.append({**s, "abs_corr": ab_corr})
    return blank, result


def detect_peak_wavelength(spectra_corr: list[dict]) -> float:
    """Auto-detect calibration wavelength from average blank-subtracted standards."""
    stds = [s for s in spectra_corr if s["is_std"]]
    if not stds:
        return 660.0  # fallback
    grid = _common_wl_grid(stds)
    avg = np.zeros_like(grid)
    for s in stds:
        avg += np.interp(grid, s["wl"], s["abs_corr"],
                         left=np.nan, right=np.nan)
    avg /= len(stds)
    valid = np.isfinite(avg)
    if not valid.any():
        return 660.0
    return float(grid[np.nanargmax(avg)])


def read_abs_at_wl(spectra_corr: list[dict], cal_wl: float) -> list[dict]:
    """Add 'abs_at_cal_wl' to each spectrum dict."""
    out = []
    for s in spectra_corr:
        ab = float(np.interp(cal_wl, s["wl"], s["abs_corr"],
                             left=np.nan, right=np.nan))
        out.append({**s, "abs_at_cal_wl": ab})
    return out


# ─────────────────────────────────────────────────────────────────────────────
#  DataFrame assembly + 3-pass calibration
# ─────────────────────────────────────────────────────────────────────────────

def build_dataframe(spectra_wl: list[dict]) -> pd.DataFrame:
    rows = []
    for s in spectra_wl:
        rows.append({
            "sample":       s["label"],
            "conc_mmol":    s["conc_mmol"],
            "abs_at_cal_wl": s["abs_at_cal_wl"],
            "_is_std":      s["is_std"],
        })
    return pd.DataFrame(rows)


def apply_uv_calibration(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """3-pass linear calibration: abs_at_cal_wl → conc_mmol."""
    std_df = df[df["_is_std"] & df["conc_mmol"].notna()
                & df["abs_at_cal_wl"].notna()].copy()

    cal_meta: dict = {"slope": math.nan, "intercept": math.nan, "r2": math.nan,
                      "areas": [], "concs": [], "source_samples": []}

    if len(std_df) >= 2:
        areas = std_df["abs_at_cal_wl"].values.astype(float)
        concs = std_df["conc_mmol"].values.astype(float)
        coeffs = np.polyfit(areas, concs, 1)
        slope, intercept = float(coeffs[0]), float(coeffs[1])
        pred   = slope * areas + intercept
        ss_res = float(np.sum((concs - pred) ** 2))
        ss_tot = float(np.sum((concs - np.mean(concs)) ** 2))
        r2     = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
        cal_meta.update({
            "slope": slope, "intercept": intercept, "r2": r2,
            "areas": list(areas), "concs": list(concs),
            "source_samples": std_df["sample"].tolist(),
        })
        # Pass 3: fill unknown concentrations
        for idx, row in df.iterrows():
            if math.isfinite(row.get("conc_mmol", math.nan)):
                continue
            ab = row.get("abs_at_cal_wl", math.nan)
            if math.isfinite(ab) and math.isfinite(slope) and math.isfinite(intercept):
                df.at[idx, "conc_mmol"] = max(slope * ab + intercept, 0.0)
    elif len(std_df) == 1:
        areas = std_df["abs_at_cal_wl"].values.astype(float)
        concs = std_df["conc_mmol"].values.astype(float)
        slope = float(concs[0] / areas[0]) if areas[0] != 0 else math.nan
        intercept = 0.0
        cal_meta.update({"slope": slope, "intercept": intercept,
                         "areas": list(areas), "concs": list(concs),
                         "source_samples": std_df["sample"].tolist()})

    return df, cal_meta


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting
# ─────────────────────────────────────────────────────────────────────────────

def _subplot_grid(n: int) -> tuple[int, int]:
    r = max(1, int(math.sqrt(n)))
    c = r if r * r == n else r + 1
    return r, c


def plot_uv_spectra(
    spectra_corr: list[dict],
    blank: dict | None,
    cal_wl: float,
    df: pd.DataFrame,
    out_path: Path,
    label: str,
) -> None:
    """Overlaid spectra: standards (ivory bg) / samples (light-blue bg), grouped."""
    stds  = [s for s in spectra_corr if s["is_std"]]
    samps = [s for s in spectra_corr if not s["is_std"]]

    fig, axes = plt.subplots(
        2, 1,
        figsize=(10, 8),
        gridspec_kw={"hspace": 0.38},
    )

    for ax_i, (group, sp_list, title_sfx, bg, label_color, palette) in enumerate([
        (stds,  stds,  "Calibration Standards", "#fffff5", "#7a5c00", _STD_COLORS),
        (samps, samps, "Experiment Samples",     "#f0f6ff", "#1a3f80", _SAMP_COLORS),
    ]):
        ax = axes[ax_i]
        ax.set_facecolor(bg)

        # Group label (no box)
        ax.text(0.003, 0.97, title_sfx,
                transform=ax.transAxes, fontsize=13, fontweight="bold",
                va="top", ha="left", color=label_color)

        for sp_i, sp in enumerate(sp_list):
            color = palette[sp_i % len(palette)]
            conc_row = df[df["sample"] == sp["label"]]
            conc_val = float(conc_row["conc_mmol"].iloc[0]) if not conc_row.empty else math.nan
            lbl_txt = sp["label"]
            if math.isfinite(conc_val):
                lbl_txt += f"  ({conc_val:.3f} mmol)"
            ax.plot(sp["wl"], sp["abs_corr"],
                    color=color, lw=1.5, label=lbl_txt, alpha=0.85)

        # Expand y upper limit to give headroom for the wavelength label
        y_lo, y_hi = ax.get_ylim()
        ax.set_ylim(y_lo, y_hi * 1.15)

        # Calibration wavelength marker
        ax.axvline(cal_wl, color="gray", ls="--", lw=1.0, alpha=0.6, zorder=5)
        ax.text(cal_wl + 2, ax.get_ylim()[1] * 0.97,
                f"{cal_wl:.0f} nm", fontsize=9, color="gray", va="top")

        ax.set_xlabel("Wavelength (nm)", fontsize=12)
        ax.set_ylabel("Abs (blank-subtracted)", fontsize=12)
        ax.tick_params(labelsize=11)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(True)
        leg = ax.legend(fontsize=9, loc="upper right", frameon=False)

    fig.suptitle(f"{label} — UV-Vis Spectra", fontsize=14)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    Saved: {out_path.name}")


def plot_uv_calibration(cal_meta: dict, cal_wl: float,
                        out_path: Path, label: str) -> None:
    """Scatter + gray dashed linear fit; R² in title."""
    areas = np.array(cal_meta.get("areas", []), dtype=float)
    concs = np.array(cal_meta.get("concs", []), dtype=float)
    if len(areas) < 2:
        return

    slope     = cal_meta.get("slope",     math.nan)
    intercept = cal_meta.get("intercept", math.nan)
    r2        = cal_meta.get("r2",        math.nan)

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    # GC/HPLC-style: tableau blue dots + gray dashed fit line
    ax.scatter(areas, concs, color="#1f77b4", s=70, zorder=4, alpha=0.6,
               edgecolors="white", linewidths=0.5, label="Cal points")

    if math.isfinite(slope) and math.isfinite(intercept):
        span = float(areas.max() - areas.min())
        # Fit line spans data range only — does NOT extend to x=0
        a_lo = float(areas.min()) - span * 0.05
        a_hi = float(areas.max()) + span * 0.05
        a_fit = np.linspace(a_lo, a_hi, 300)
        ax.plot(a_fit, slope * a_fit + intercept,
                color="gray", lw=1.8, ls="--", alpha=0.7, zorder=3,
                label="Linear fit")


    ax.set_xlabel(f"Absorbance at {cal_wl:.0f} nm (AU)", fontsize=10)
    ax.set_ylabel("Concentration (mmol)", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
    leg = ax.legend(fontsize=9, loc="best")
    leg.get_frame().set_visible(False)

    sl_s = f"{slope:.4e}" if math.isfinite(slope) else "N/A"
    ic_s = f"{intercept:+.4e}" if math.isfinite(intercept) else ""
    r2_s = f"R² = {r2:.4f}" if math.isfinite(r2) else "R² = N/A"
    ax.set_title(
        f"UV-Vis Calibration at {cal_wl:.0f} nm\n"
        f"y = {sl_s}·x{ic_s}    {r2_s}",
        fontsize=9,
    )
    fig.suptitle(f"{label} — UV-Vis Calibration", fontsize=11)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved: {out_path.name}")


def plot_uv_concentrations(df: pd.DataFrame, out_path: Path, label: str) -> None:
    """Concentration vs sample name (experiment samples only), PCHIP curve."""
    samp_df = df[~df["_is_std"]].copy().reset_index(drop=True)
    if samp_df.empty:
        return

    x = np.arange(1, len(samp_df) + 1, dtype=float)
    vals = samp_df["conc_mmol"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(max(6, len(samp_df) * 0.9 + 2), 4))
    ax.scatter(x, vals, color="#d62728", s=50, zorder=4, alpha=0.85,
               label="Sample concentration")

    valid = np.isfinite(vals)
    if valid.sum() >= 2:
        xv, yv = x[valid], vals[valid]
        xi = np.linspace(xv[0], xv[-1], 300)
        ax.plot(xi, PchipInterpolator(xv, yv)(xi),
                color="#d62728", lw=1.5, alpha=0.8)

    ax.set_xticks(list(x))
    ax.set_xticklabels(samp_df["sample"].tolist(),
                       rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Concentration (mmol)", fontsize=11)
    ax.set_title(f"{label} — UV-Vis Sample Concentrations", fontsize=11)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
    leg = ax.legend(fontsize=9, loc="best")
    leg.get_frame().set_visible(False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved: {out_path.name}")


# ─────────────────────────────────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cary CSV UV-Vis 3-pass adaptive calibration")
    parser.add_argument("--csv",     type=Path, required=True,
                        help="Path to Cary multi-spectrum CSV export")
    parser.add_argument("--out-dir", type=Path, default=Path("output/uv"),
                        help="Output directory")
    parser.add_argument("--label",   default="UV-Vis",
                        help="Title label for plots")
    parser.add_argument("--cal-wl",  type=float, default=None,
                        help="Calibration wavelength (nm). "
                             "Auto-detected from standards if omitted.")
    args = parser.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    import shutil
    _s = out / "scripts"; _s.mkdir(exist_ok=True)
    shutil.copy2(__file__, _s / Path(__file__).name)

    print(f"[1] Parsing {args.csv.name}")
    spectra = parse_cary_csv(args.csv)
    print(f"    {len(spectra)} spectra found:")
    for sp in spectra:
        flag = "STD" if sp["is_std"] else ("BLK" if sp["is_blank"] else "SMP")
        conc_s = f"  {sp['conc_mmol']:.4f} mmol" if sp["conc_mmol"] is not None else ""
        print(f"      [{flag}]  {sp['label']:20s}  "
              f"wl {sp['wl'].min():.0f}–{sp['wl'].max():.0f} nm{conc_s}")

    print("[2] Blank subtraction")
    blank, spectra_corr = blank_subtract(spectra)
    print(f"    Blank: {blank['label'] if blank else 'none'}")

    if args.cal_wl is not None:
        cal_wl = args.cal_wl
        print(f"    Calibration wavelength: {cal_wl:.1f} nm (user-specified)")
    else:
        cal_wl = detect_peak_wavelength(spectra_corr)
        print(f"    Calibration wavelength: {cal_wl:.1f} nm (auto-detected)")

    print("[3] Reading absorbance at calibration wavelength")
    spectra_wl = read_abs_at_wl(spectra_corr, cal_wl)
    df = build_dataframe(spectra_wl)

    cols = ["sample", "conc_mmol", "abs_at_cal_wl"]
    print(df[cols].to_string())

    print("[4] 3-pass calibration")
    df, cal_meta = apply_uv_calibration(df)
    sl   = cal_meta.get("slope",     math.nan)
    ic_  = cal_meta.get("intercept", math.nan)
    r2   = cal_meta.get("r2",        math.nan)
    n_   = len(cal_meta.get("areas", []))
    print(f"    y = {sl:.4e}·x {ic_:+.4e}   R² = {r2:.4f}   ({n_} pts)")

    print("[5] Saving plots and CSV")
    plot_uv_spectra(spectra_corr, blank, cal_wl, df,
                    out / "uv_spectra.png", args.label)
    plot_uv_calibration(cal_meta, cal_wl,
                        out / "uv_calibration.png", args.label)
    plot_uv_concentrations(df, out / "uv_concentrations.png", args.label)

    df.drop(columns=["_is_std"], errors="ignore").to_csv(
        out / "uv_raw.csv", index=False)
    print(f"    Saved: uv_raw.csv")

    print(f"\n✓ Done. Output: {out.resolve()}")

    # Print sample concentration summary
    samp_df = df[~df["_is_std"]].copy()
    if not samp_df.empty:
        print("\nSample concentrations:")
        print(samp_df[["sample", "conc_mmol", "abs_at_cal_wl"]].to_string(index=False))


if __name__ == "__main__":
    main()
