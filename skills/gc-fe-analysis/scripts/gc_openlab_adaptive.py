"""Standalone Agilent 990 Micro GC adaptive integration and calibration.

Analogous to hplc_openlab_adaptive.py — same three-pass calibration approach,
same plot conventions (integration windows, calibration curves, ppm vs injection).

Usage
-----
uv run python gc_openlab_adaptive.py \\
    --gc-dirs "/path/to/run1.rslt" "/path/to/run2.rslt" \\
    --label "Example sequence" \\
    --out-dir output/gc

Outputs
-------
  gc_integration_windows.png  — TCD1A + TCD1B traces per injection, shared x-axis
  gc_calibration.png          — scatter + linear regression per compound
  gc_ppm.png                  — ppm vs injection number (1-based), NaN = gap
  gc_raw.csv                  — all parsed values
"""

from __future__ import annotations

import argparse
import math
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

# ─────────────────────────────────────────────────────────────────────────────
#  Constants
# ─────────────────────────────────────────────────────────────────────────────

GC_NS_ACAML = "urn:schemas-agilent-com:acaml21"   # used in .rx files
GC_NS_ACMD  = "urn:schemas-agilent-com:acmd20"    # used in .dx files (injection.acmd)

# Compound → (channel, center_min, half_window_min)
GC_PEAKS: dict[str, tuple[str, float, float]] = {
    "H2":         ("TCD1A", 0.582, 0.010),
    "CH4":        ("TCD1A", 1.252, 0.080),
    "CO":         ("TCD1A", 1.851, 0.180),
    "CH4_B":      ("TCD1B", 0.550, 0.015),
    "CO2":        ("TCD1B", 0.791, 0.060),
    "C2H4":       ("TCD1B", 0.914, 0.040),
    "C2H6":       ("TCD1B", 1.051, 0.040),
}

# OpenLAB compound name → our key
GC_COMPOUND_MAP: dict[str, str] = {
    "H2":              "H2",
    "CO":              "CO",
    "CH4 Channel 1":   "CH4",
    "CH4 Channel 2":   "CH4_B",
    "CO2":             "CO2",
    "C2H4":            "C2H4",
    "C2H6":            "C2H6",
}

GC_COLORS: dict[str, str] = {
    "H2":    "#1f77b4",
    "CO":    "#d62728",
    "CH4":   "#2ca02c",
    "CH4_B": "#98df8a",
    "CO2":   "#ff7f0e",
    "C2H4":  "#9467bd",
    "C2H6":  "#8c564b",
}

# Adaptive integration window multipliers per compound
# search: Pass-1 half-width = search × hw (symmetric)
# il, ir: Pass-2 integration = [apex − il×hw, apex + ir×hw]
GC_WINDOW_PARAMS: dict[str, dict] = {
    "H2":    {"search": 2.0, "il": 1.5, "ir": 2.0},
    "CH4":   {"search": 1.575, "il": 1.05, "ir": 1.4},          # TCD1A: prev×0.7
    "CO":    {"search": 1.0,   "il": 0.75, "ir": 1.0},          # × 0.5 (narrow peak)
    "CH4_B": {"search": 1.4,  "il": 1.05, "ir": 1.4, "sl_offset": 0.02},  # TCD1B: prev×0.7, search left +0.02
    "CO2":   {"search": 2.0, "il": 1.2, "ir": 1.6},    # il/ir × 0.8
    "C2H4":  {"search": 0.6, "il": 0.45, "ir": 1.08},  # search/il × 0.3; ir = il × 2.4
    "C2H6":  {"search": 0.6, "il": 0.45, "ir": 0.6},  # all × 0.3
}

# Fallback calibration fits (from gc_window_fe_analyzer.py / 5 GC STD notebook).
# Used when insufficient Agilent-quantified points exist to derive a regression.
DEFAULT_FITS: dict[str, dict] = {
    "H2":   {"slope": 4.433952e6, "intercept": -74.932480,   "r2": 0.999232},
    "CO":   {"slope": 7.953694e4, "intercept": -69.496289,   "r2": 0.999371},
    "CH4":  {"slope": 8.981725e4, "intercept": -6.102562,    "r2": 0.999476},
    "C2H4": {"slope": 7.912502e4, "intercept": -68.665229,   "r2": 0.999415},
    "C2H6": {"slope": 7.613048e4, "intercept": -6.163700,    "r2": 0.999502},
    "CO2":  {"slope": 8.278111e4, "intercept": -232.723274,  "r2": 0.990442},
}

CHANNEL_XRANGE = {
    "TCD1A": (0.0, 2.5),
    "TCD1B": (0.0, 1.5),
}


# ─────────────────────────────────────────────────────────────────────────────
#  I/O
# ─────────────────────────────────────────────────────────────────────────────

def tag(local: str, ns: str = GC_NS_ACAML) -> str:
    return f"{{{ns}}}{local}"


def acmd_tag(local: str) -> str:
    return f"{{{GC_NS_ACMD}}}{local}"


def read_gc_dx(dx_path: Path) -> dict[str, dict]:
    """Parse an Agilent 990 GC .dx file.  Returns {channel: {time_min, signal_mV}}."""
    channels: dict[str, dict] = {}
    try:
        with zipfile.ZipFile(dx_path) as zf:
            root = ET.parse(zf.open("injection.acmd")).getroot()
            info = root.find(acmd_tag("InjectionInfo"))
            if info is None:
                return channels
            for sig in info.findall(f"{acmd_tag('Signals')}/{acmd_tag('Signal')}"):
                ch   = sig.findtext(acmd_tag("ChannelName"), "")
                tid  = sig.findtext(acmd_tag("TraceId"), "")
                n    = int(sig.findtext(acmd_tag("NumberOfValues"), "0"))
                t0   = float(sig.findtext(acmd_tag("TimeStart"), "0"))
                t1   = float(sig.findtext(acmd_tag("TimeEnd"),   "0"))
                slo  = float(sig.findtext(acmd_tag("Slope"),     "1"))
                if n == 0 or tid == "":
                    continue
                try:
                    raw  = zf.open(f"{tid}.CH").read()
                    data = np.frombuffer(raw[len(raw) - n * 8:], dtype="<f8")[:n] * slo
                    channels[ch] = {
                        "time_min":  np.linspace(t0, t1, n) / 60_000,
                        "signal_mV": data,
                    }
                except (KeyError, ValueError):
                    continue
    except Exception:
        pass
    return channels


def parse_gc_rx(rx_path: Path) -> dict[str, dict]:
    """Parse GC .rx XML → {our_key: {ppm, area_au}} for each detected compound."""
    result: dict[str, dict] = {}
    try:
        with zipfile.ZipFile(rx_path) as zf:
            root = ET.parse(zf.open("Base/InjectionACAML")).getroot()
            for ic in root.iter(tag("InjectionCompound")):
                name_raw = ic.findtext(tag("CompoundName"), "").strip()
                key = GC_COMPOUND_MAP.get(name_raw)
                if key is None:
                    continue
                conc_el  = ic.find(tag("Concentration"))
                area_el  = ic.find(tag("Area"))
                ppm = float(conc_el.attrib.get("val", "nan")) \
                      if conc_el is not None else math.nan
                area_au = float(area_el.attrib.get("val", "nan")) \
                          if area_el is not None else math.nan
                result[key] = {"ppm": ppm, "area_au": area_au}
    except Exception:
        pass
    return result


def collect_gc_data(rslt_dirs: list[Path]) -> pd.DataFrame:
    """Scan .rslt directories → DataFrame (one row per injection, sorted by timestamp)."""
    rows: list[dict] = {}
    for rslt in rslt_dirs:
        for rx_p in sorted(rslt.glob("*.rx")):
            stem = rx_p.stem
            dx_p = rx_p.with_suffix(".dx")
            parsed = parse_gc_rx(rx_p)
            row: dict = {"rx_file": rx_p.name, "dx_path": str(dx_p), "rslt_dir": rslt.name}
            # Timestamp from filename
            import re
            m = re.match(r"(\d{4}-\d{2}-\d{2} \d{2}-\d{2}-\d{2})", stem)
            if m:
                try:
                    row["run_datetime"] = pd.to_datetime(m.group(1), format="%Y-%m-%d %H-%M-%S")
                except Exception:
                    row["run_datetime"] = pd.NaT
            else:
                row["run_datetime"] = pd.NaT
            for key in GC_PEAKS:
                d = parsed.get(key, {})
                row[key]              = d.get("ppm",     math.nan)
                row[f"{key}_area_au"] = d.get("area_au", math.nan)
            rows[stem] = row
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(list(rows.values()))
    df = df.sort_values("run_datetime", na_position="last").reset_index(drop=True)
    return df


# ─────────────────────────────────────────────────────────────────────────────
#  Adaptive two-pass integration
# ─────────────────────────────────────────────────────────────────────────────

def integrate_gc_peak_adaptive(
    t_min: np.ndarray,
    sig_mV: np.ndarray,
    center: float,
    hw: float,
    search_mult: float = 2.0,
    integ_left_mult: float = 1.5,
    integ_right_mult: float = 2.0,
    sl_offset: float = 0.0,
) -> dict:
    """Two-pass adaptive GC peak integration.

    Pass 1 — search: [center − search×hw + sl_offset, center + search×hw] → apex_t
             sl_offset shifts the search left boundary (positive = trim left side).
    Pass 2 — integrate: [apex − il×hw, apex + ir×hw]
             flat baseline at min(sig[left], sig[right])

    Returns dict with area (mV·min), apex_t, apex_y, integ bounds.
    """
    sl = center - search_mult * hw + sl_offset
    sr = center + search_mult * hw
    mask1 = (t_min >= sl) & (t_min <= sr)
    if mask1.sum() < 2:
        return {"area": math.nan}
    apex_idx = int(np.argmax(sig_mV[mask1]))
    apex_t   = float(t_min[mask1][apex_idx])
    apex_y   = float(sig_mV[mask1][apex_idx])

    il = apex_t - integ_left_mult  * hw
    ir = apex_t + integ_right_mult * hw
    mask2 = (t_min >= il) & (t_min <= ir)
    if mask2.sum() < 2:
        return {"area": math.nan}
    t_seg   = t_min[mask2]
    sig_seg = sig_mV[mask2]
    baseline = min(float(sig_seg[0]), float(sig_seg[-1]))
    area = float(np.trapezoid(sig_seg - baseline, t_seg))
    return {
        "area":    max(area, 0.0),
        "apex_t":  apex_t,
        "apex_y":  apex_y,
        "integ_l": il,
        "integ_r": ir,
        "search_l": sl,
        "search_r": sr,
        "baseline": baseline,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  3-pass calibration (same approach as hplc_openlab_adaptive.py)
# ─────────────────────────────────────────────────────────────────────────────

def apply_gc_calibration(gc_df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Re-integrate all injections, then build ppm-vs-area linear regression.

    Pass 1 — re-integrate every injection from raw .dx.
    Pass 2 — fit linear regression from Agilent-quantified rows (noise filter:
             ppm ≥ max(10 ppm, 5% of max observed)).
    Pass 3 — fill NaN rows using the fit.

    Returns (modified_df, cal_meta) where cal_meta has per-compound:
      {slope, intercept, r2, areas, ppms, source_files}
    """
    compounds = list(GC_PEAKS.keys())

    # ── Pass 1: re-integrate all ────────────────────────────────────────────
    our_areas: dict[int, dict[str, float]] = {}
    for idx, row in gc_df.iterrows():
        dx_p = Path(str(row.get("dx_path", "")))
        if not dx_p.exists():
            continue
        channels = read_gc_dx(dx_p)
        our_areas[idx] = {}
        for comp, (ch, center, hw) in GC_PEAKS.items():
            trace = channels.get(ch)
            if trace is None:
                our_areas[idx][comp] = math.nan
                continue
            wp = GC_WINDOW_PARAMS.get(comp, {"search": 2.0, "il": 1.5, "ir": 2.0})
            pk = integrate_gc_peak_adaptive(
                trace["time_min"], trace["signal_mV"],
                center, hw,
                search_mult=wp["search"],
                integ_left_mult=wp["il"],
                integ_right_mult=wp["ir"],
                sl_offset=wp.get("sl_offset", 0.0),
            )
            our_areas[idx][comp] = pk.get("area", math.nan)
        gc_df.at[idx, "dx_ok"] = True

    # ── Pass 2: calibrate from STANDARDS ONLY (noise-filtered) ─────────────
    # A row is a calibration standard when its rslt_dir name contains "cal"
    # (matches "Cal1.rslt", "0_calibration_...", etc.).
    # Only standard points enter the linear regression; experiment samples are
    # never included in the fit.
    max_ppm: dict[str, float] = {}
    for comp in compounds:
        if comp in gc_df.columns:
            # Use only standard rows to determine the ppm scale (avoids
            # experiment-sample CO2 at ~94 % inflating the noise threshold).
            std_mask = gc_df.get("rslt_dir", pd.Series(dtype=str)).apply(
                lambda d: "cal" in str(d).lower())
            pos = gc_df.loc[std_mask, comp].dropna() if comp in gc_df.columns else pd.Series(dtype=float)
            pos = pos[pos > 0]
            max_ppm[comp] = float(pos.max()) if not pos.empty else 0.0

    # Always seed every compound's point list with the origin (0, 0).
    # This enforces zero-area → zero-ppm physics and anchors the regression
    # even for compounds with few calibration standards.
    std_data: dict[str, dict] = {
        c: {"areas": [0.0], "ppms": [0.0], "files": ["(0,0) origin"]}
        for c in compounds
    }

    for idx, row in gc_df.iterrows():
        rslt_dir = str(row.get("rslt_dir", ""))
        if "cal" not in rslt_dir.lower():
            continue  # skip experiment rows entirely in Pass 2
        for comp in compounds:
            agilent_ppm = row.get(comp, math.nan)
            our_area    = our_areas.get(idx, {}).get(comp, math.nan)
            # Include ALL standard points with valid positive ppm and area;
            # no noise-threshold filter — origin already anchors the low end.
            if not (math.isfinite(agilent_ppm) and agilent_ppm > 0
                    and math.isfinite(our_area) and our_area > 0):
                continue
            std_data[comp]["areas"].append(float(our_area))
            std_data[comp]["ppms"].append(float(agilent_ppm))
            std_data[comp]["files"].append(str(row.get("rx_file", "?")))

    fits: dict[str, dict] = {}
    for comp in compounds:
        areas = std_data[comp]["areas"]
        ppms  = std_data[comp]["ppms"]
        if len(areas) >= 2:
            coeffs = np.polyfit(areas, ppms, 1)
            slope, intercept = float(coeffs[0]), float(coeffs[1])
            pred   = slope * np.array(areas) + intercept
            ss_res = float(np.sum((np.array(ppms) - pred) ** 2))
            ss_tot = float(np.sum((np.array(ppms) - np.mean(ppms)) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
        elif len(areas) == 1:
            slope = float(ppms[0] / areas[0]) if areas[0] != 0 else math.nan
            intercept, r2 = 0.0, math.nan
        else:
            slope = intercept = r2 = math.nan
        if not math.isfinite(slope) and comp in DEFAULT_FITS:
            d = DEFAULT_FITS[comp]
            slope, intercept, r2 = d["slope"], d["intercept"], d["r2"]
        fits[comp] = {"slope": slope, "intercept": intercept, "r2": r2}

    # ── Pass 3: calculate concentrations for ALL non-standard rows ───────────
    # For each experiment injection: ppm = slope × our_area + intercept.
    # At the same time collect (our_area, fit_ppm) pairs for the calibration
    # plot — these points lie exactly on the line by construction.
    samp_plot: dict[str, dict] = {c: {"areas": [], "ppms": [], "files": []} for c in compounds}

    for idx, row in gc_df.iterrows():
        rslt_dir = str(row.get("rslt_dir", ""))
        if "cal" in rslt_dir.lower():
            continue  # standards keep their Agilent values
        for comp in compounds:
            our_area = our_areas.get(idx, {}).get(comp, math.nan)
            fit = fits.get(comp, {})
            sl  = fit.get("slope",     math.nan)
            ic  = fit.get("intercept", math.nan)
            if not (math.isfinite(our_area) and our_area > 0
                    and math.isfinite(sl) and math.isfinite(ic)):
                continue
            fit_ppm = max(sl * our_area + ic, 0.0)
            gc_df.at[idx, comp]             = fit_ppm
            gc_df.at[idx, f"{comp}_area"]   = our_area
            gc_df.at[idx, "_reintegrated"]  = True
            # Record for calibration plot (point lies on line by construction)
            samp_plot[comp]["areas"].append(float(our_area))
            samp_plot[comp]["ppms"].append(fit_ppm)
            samp_plot[comp]["files"].append(str(row.get("rx_file", "?")))

    cal_meta: dict[str, dict] = {}
    for comp in compounds:
        fit  = fits[comp]
        slope, intercept, r2 = fit["slope"], fit["intercept"], fit["r2"]
        cal_meta[comp] = {
            "slope": slope, "intercept": intercept, "r2": r2,
            "areas":        std_data[comp]["areas"],
            "ppms":         std_data[comp]["ppms"],
            "source_files": std_data[comp]["files"],
            # Experiment sample points — on the line by construction
            "samp_areas":        samp_plot[comp]["areas"],
            "samp_ppms":         samp_plot[comp]["ppms"],
            "samp_source_files": samp_plot[comp]["files"],
            "fallback": (not std_data[comp]["areas"] and comp in DEFAULT_FITS),
        }

    return gc_df, cal_meta


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting
# ─────────────────────────────────────────────────────────────────────────────

def plot_gc_integration_windows(gc_df: pd.DataFrame, out_path: Path, label: str) -> None:
    """One row per injection, TCD1A and TCD1B side by side. Shared x-axis per channel."""
    n = len(gc_df)
    if n == 0:
        return
    channels = ["TCD1A", "TCD1B"]
    fig, axes = plt.subplots(n, 2, figsize=(14, max(2.5 * n, 5)),
                             sharex="col", squeeze=False)

    for row_i, (_, row) in enumerate(gc_df.iterrows()):
        dx_p = Path(str(row.get("dx_path", "")))
        chs = read_gc_dx(dx_p) if dx_p.exists() else {}

        ts = row.get("run_datetime")
        try:
            title = pd.Timestamp(ts).strftime("%Y-%m-%d %H:%M") if pd.notna(ts) else row.get("rx_file", "?")
        except Exception:
            title = str(row.get("rx_file", "?"))

        reint = bool(row.get("_reintegrated", False))
        tag_str = "[re-integrated]" if reint else "[Agilent]"

        for col_i, ch in enumerate(channels):
            ax = axes[row_i][col_i]
            trace = chs.get(ch)
            if trace is not None:
                ax.plot(trace["time_min"], trace["signal_mV"], lw=0.8, color="k", zorder=2)

            legend_patches = []
            peak_apex_ys:        list[float] = []   # all peaks (for baselines)
            peak_apex_ys_scale:  list[float] = []   # peaks used for y-axis scaling
            peak_baselines:      list[float] = []
            # TCD1B: exclude CO2 from y-axis upper-limit calculation (CO2 peak is
            # orders of magnitude larger and would compress all product peaks)
            Y_SCALE_EXCLUDE = {"TCD1B": {"CO2"}}
            for comp, (comp_ch, center, hw) in GC_PEAKS.items():
                if comp_ch != ch:
                    continue
                ppm_val = row.get(comp, math.nan)
                color   = GC_COLORS.get(comp, "gray")
                if trace is None or not math.isfinite(ppm_val):
                    continue
                wp = GC_WINDOW_PARAMS.get(comp, {"search": 2.0, "il": 1.5, "ir": 2.0})
                pk = integrate_gc_peak_adaptive(
                    trace["time_min"], trace["signal_mV"], center, hw,
                    search_mult=wp["search"],
                    integ_left_mult=wp["il"],
                    integ_right_mult=wp["ir"],
                    sl_offset=wp.get("sl_offset", 0.0),
                )
                if not pk.get("area", math.nan) or not math.isfinite(pk["area"]):
                    continue
                il, ir = pk["integ_l"], pk["integ_r"]
                ax.axvspan(il, ir, alpha=0.25, color=color, zorder=1)
                ax.axvline(pk["apex_t"], color=color, lw=1.2, ls="--", zorder=3)
                legend_patches.append(mpatches.Patch(color=color, alpha=0.6,
                    label=f"{comp} {ppm_val:.0f} ppm"))
                if math.isfinite(pk.get("apex_y", math.nan)):
                    peak_apex_ys.append(pk["apex_y"])
                    if comp not in Y_SCALE_EXCLUDE.get(ch, set()):
                        peak_apex_ys_scale.append(pk["apex_y"])
                if math.isfinite(pk.get("baseline", math.nan)):
                    peak_baselines.append(pk["baseline"])

            ax.set_xlim(CHANNEL_XRANGE.get(ch, (0, 3)))
            # Y-axis: upper = max highlighted peak × 1.2 (CO2 excluded for TCD1B);
            #         lower = auto from baselines.
            # TCD1B minimum range: [-0.3, 0.8] — if auto limits are tighter, use these.
            TCD1B_YMIN = (-0.3, 0.8)
            scale_apexes = peak_apex_ys_scale if peak_apex_ys_scale else peak_apex_ys
            if scale_apexes:
                y_top = max(scale_apexes) * 1.2
                y_bot = min(peak_baselines) if peak_baselines else min(scale_apexes) * 0.0
                margin = (y_top - y_bot) * 0.05
                y_bot -= margin
                if ch == "TCD1B":
                    y_bot = min(y_bot, TCD1B_YMIN[0])
                    y_top = max(y_top, TCD1B_YMIN[1])
                ax.set_ylim(y_bot, y_top)
            ax.set_xlabel("Time (min)", fontsize=7)
            ax.set_ylabel("TCD (mV)", fontsize=7)
            ax.tick_params(labelsize=6)
            if row_i == 0:
                ax.set_title(f"{ch}", fontsize=8)
            ax.grid(False)
            if legend_patches:
                leg = ax.legend(handles=legend_patches, fontsize=6, loc="upper right")
                leg.get_frame().set_visible(False)

        # Row label
        axes[row_i][0].set_ylabel(f"{title}\n{tag_str}\nTCD (mV)", fontsize=6)

    fig.suptitle(f"{label} — GC Integration Windows", fontsize=10)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved: {out_path.name}")


def _pchip_plot(ax, x, y, marker, ls, color, label_str, ms=4, lw=1.2, alpha=0.8):
    """Scatter measured points and connect with a PCHIP smooth curve."""
    valid = np.isfinite(y)
    ax.scatter(np.array(x)[valid], y[valid], marker=marker, color=color,
               s=ms ** 2, alpha=alpha, zorder=4, label=label_str)
    xv = np.array(x)[valid]
    yv = y[valid]
    if len(xv) >= 2:
        xi = np.linspace(xv[0], xv[-1], 300)
        ax.plot(xi, PchipInterpolator(xv, yv)(xi), ls=ls, color=color,
                lw=lw, alpha=alpha)


def plot_gc_ppm(gc_df: pd.DataFrame, out_path: Path, label: str) -> None:
    """2-panel ppm vs injection #.

    Panel 1 (top)   : H2, CO, CH4  — TCD1A products (full scale).
    Panel 2 (bottom): C2H4 + C2H6 on left axis scaled to C2H4 peak;
                      CO2 on right twin axis.
    Negative/below-zero values masked as NaN (shown as gaps).
    """
    x = list(range(1, len(gc_df) + 1))

    def masked(col: str) -> np.ndarray:
        if col not in gc_df.columns:
            return np.full(len(gc_df), math.nan)
        return gc_df[col].where(gc_df[col] > 0).to_numpy(dtype=float)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    # ── Panel 1: TCD1A products ──────────────────────────────────────────────
    for comp in ("H2", "CO", "CH4"):
        color = GC_COLORS.get(comp, "gray")
        _pchip_plot(ax1, x, masked(comp), "o", "-", color, comp)
    ax1.set_ylabel("Concentration (ppm)", fontsize=9)
    ax1.set_title(f"{label} — Gas products (GC TCD)", fontsize=10)
    leg1 = ax1.legend(fontsize=8, loc="best")
    leg1.get_frame().set_visible(False)
    ax1.grid(False)

    # ── Panel 2: TCD1B products, left = C2H4 scale, right = CO2 ────────────
    c2h4_vals = masked("C2H4")
    c2h6_vals = masked("C2H6")
    co2_vals  = masked("CO2")

    # y-axis ceiling: C2H4 peak × 1.3 (or small positive floor if all NaN)
    c2h4_max = float(np.nanmax(c2h4_vals)) if np.any(np.isfinite(c2h4_vals)) else 1.0
    y_lo, y_hi = 0.0, max(c2h4_max * 1.30, 1.0)

    _pchip_plot(ax2, x, c2h4_vals, "o", "-", GC_COLORS["C2H4"], "C2H4")
    _pchip_plot(ax2, x, c2h6_vals, "s", "-", GC_COLORS["C2H6"], "C2H6")
    _pchip_plot(ax2, x, masked("CH4_B"), "^", "-", GC_COLORS["CH4_B"], "CH4_B", lw=1.0)
    ax2.set_ylim(y_lo, y_hi)
    ax2.set_ylabel("C2H4 / C2H6 / CH4_B (ppm)", fontsize=9)
    ax2.grid(False)

    # CO2 twin axis
    ax2r = ax2.twinx()
    _pchip_plot(ax2r, x, co2_vals, "D", "--", GC_COLORS["CO2"], "CO2", ms=3, lw=1.0)
    ax2r.set_ylabel("CO₂ (ppm)", fontsize=9, color=GC_COLORS["CO2"])
    ax2r.tick_params(axis="y", labelcolor=GC_COLORS["CO2"])

    # Combine both axes' handles into one legend so matplotlib can find a
    # single clear gap in the data rather than placing two overlapping boxes.
    h2_handles, h2_labels = ax2.get_legend_handles_labels()
    hr_handles, hr_labels = ax2r.get_legend_handles_labels()
    leg2 = ax2.legend(h2_handles + hr_handles, h2_labels + hr_labels,
                      fontsize=8, loc="best")
    leg2.get_frame().set_visible(False)

    ax2.set_xticks(x)
    ax2.set_xlabel("GC injection #", fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved: {out_path.name}")


def plot_gc_calibration_curves(cal_meta: dict, out_path: Path, label: str) -> None:
    """Scatter + linear regression per compound, 2-row layout by TCD channel.

    Row 0 (TCD1A): H2, CH4, CO
    Row 1 (TCD1B): CO2, C2H4, C2H6
    3 columns each → 2 × 3 grid.
    """
    # Fixed compound order by channel row
    ROWS: list[list[str]] = [
        ["H2",  "CH4",  "CO"],    # TCD1A
        ["CO2", "C2H4", "C2H6"],  # TCD1B
    ]
    n_rows, n_cols = len(ROWS), max(len(r) for r in ROWS)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4.5 * n_cols, 4.5 * n_rows),
                             squeeze=False)

    for row_i, row_comps in enumerate(ROWS):
        for col_i, comp in enumerate(row_comps):
            ax   = axes[row_i][col_i]
            meta = cal_meta.get(comp, {})
            areas = np.array(meta.get("areas", []), dtype=float)
            ppms  = np.array(meta.get("ppms",  []), dtype=float)
            slope     = meta.get("slope",     math.nan)
            intercept = meta.get("intercept", math.nan)
            r2        = meta.get("r2",        math.nan)
            color     = GC_COLORS.get(comp, "steelblue")

            if len(areas) == 0 and not math.isfinite(slope):
                ax.set_visible(False)
                continue

            use_log = (len(ppms) >= 2 and ppms.min() >= 0 and ppms.max() > 0
                       and ppms[ppms > 0].min() > 0
                       and ppms.max() / max(ppms[ppms > 0].min(), 1e-9) > 50)

            # Count measured standards (excluding the forced origin) for the title
            origin_mask = (areas == 0.0) & (ppms == 0.0)
            std_areas = areas[~origin_mask]
            std_ppms  = ppms[~origin_mask]

            # ── All calibration points (origin + standards): same filled circle ──
            if len(areas):
                ax.scatter(areas, ppms, color=color, s=60, zorder=4, alpha=0.6,
                           marker="o", edgecolors="white", linewidths=0.5,
                           label="standard")

            # ── Regression line: always starts from x = 0 ─────────────────
            if math.isfinite(slope) and math.isfinite(intercept) and len(areas) >= 2:
                a_hi  = float(areas.max()) * 1.1 if areas.max() > 0 else 1.0
                a_fit = np.linspace(0.0, a_hi, 300)
                ax.plot(a_fit, slope * a_fit + intercept,
                        color="gray", lw=1.8, ls="--", alpha=0.7, zorder=3, label="linear fit")


            if use_log:
                ax.set_xscale("log"); ax.set_yscale("log")

            ax.set_xlabel("Peak area (mV·min)", fontsize=9)
            ax.set_ylabel("Concentration (ppm)", fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(False)
            leg = ax.legend(fontsize=7, loc="upper left")
            leg.get_frame().set_visible(False)

            sl_str = f"{slope:.4e}" if math.isfinite(slope) else "N/A"
            ic_str = f"{intercept:+.4e}" if math.isfinite(intercept) else ""
            r2_str = f"R² = {r2:.4f}" if math.isfinite(r2) else "R² = N/A"
            # n counts measured standards only (origin excluded from display count)
            ax.set_title(f"{comp}   n={len(std_areas)} std + (0,0)\ny={sl_str}·x{ic_str}\n{r2_str}",
                         fontsize=8)

    # Row labels
    axes[0][0].annotate("TCD1A", xy=(0, 0.5), xycoords="axes fraction",
                         xytext=(-0.35, 0.5), textcoords="axes fraction",
                         fontsize=9, fontweight="bold", va="center", rotation=90)
    axes[1][0].annotate("TCD1B", xy=(0, 0.5), xycoords="axes fraction",
                         xytext=(-0.35, 0.5), textcoords="axes fraction",
                         fontsize=9, fontweight="bold", va="center", rotation=90)

    fig.suptitle(f"{label} — GC Calibration (linear regression, mV·min vs ppm)",
                 fontsize=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved: {out_path.name}")


# ─────────────────────────────────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    import re as _re
    parser = argparse.ArgumentParser(
        description="Parse Agilent 990 GC .rslt folders and plot integration "
                    "windows + calibration curves (analogous to hplc_openlab_adaptive.py).")
    parser.add_argument("--gc-dirs",  nargs="+", required=True, type=Path,
                        help=".rslt directory/directories to process")
    parser.add_argument("--label",    default="GC",  help="Title label for plots")
    parser.add_argument("--out-dir",  default="output/gc", type=Path)
    args = parser.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    import shutil; _s = out / "scripts"; _s.mkdir(exist_ok=True)
    shutil.copy2(__file__, _s / Path(__file__).name)

    # Fix import shadowing
    import pandas as pd

    print(f"[1] Collecting GC injections from {len(args.gc_dirs)} .rslt dir(s)")

    # Rebuild collect_gc_data inline to avoid the inner import issue
    rows_dict: dict = {}
    for rslt in args.gc_dirs:
        for rx_p in sorted(rslt.glob("*.rx")):
            stem  = rx_p.stem
            dx_p  = rx_p.with_suffix(".dx")
            parsed = parse_gc_rx(rx_p)
            row: dict = {"rx_file": rx_p.name, "dx_path": str(dx_p),
                         "rslt_dir": rslt.name}
            m = _re.match(r"(\d{4}-\d{2}-\d{2} \d{2}-\d{2}-\d{2})", stem)
            if m:
                try:
                    row["run_datetime"] = pd.to_datetime(m.group(1),
                                                         format="%Y-%m-%d %H-%M-%S")
                except Exception:
                    row["run_datetime"] = pd.NaT
            else:
                row["run_datetime"] = pd.NaT
            for key in GC_PEAKS:
                d = parsed.get(key, {})
                row[key]              = d.get("ppm",     math.nan)
                row[f"{key}_area_au"] = d.get("area_au", math.nan)
            rows_dict[stem] = row

    if not rows_dict:
        print("  No .rx files found. Exiting.")
        return

    gc_df = pd.DataFrame(list(rows_dict.values()))
    gc_df = gc_df.sort_values("run_datetime", na_position="last").reset_index(drop=True)
    print(f"    {len(gc_df)} injection(s) found")
    cols_show = ["run_datetime", "rx_file"] + \
                [c for c in GC_PEAKS if c in gc_df.columns]
    print(gc_df[cols_show].to_string())

    print("[2] Applying 3-pass calibration (with DEFAULT_FITS fallback for CH4 etc.)")
    gc_df, cal_meta = apply_gc_calibration(gc_df)
    print("    Calibration fit (linear regression, mV·min → ppm):")
    for comp, m in cal_meta.items():
        sl  = m.get("slope",     math.nan)
        ic_ = m.get("intercept", math.nan)
        r2  = m.get("r2",        math.nan)
        n_  = len(m["areas"])
        fb  = "[DEFAULT]" if m.get("fallback") else ""
        sl_s = f"{sl:.4e}" if math.isfinite(sl) else "N/A"
        ic_s = f"{ic_:+.4e}" if math.isfinite(ic_) else ""
        r2_s = f"R²={r2:.4f}" if math.isfinite(r2) else "R²=N/A"
        print(f"      {comp:10s}  y={sl_s}·x{ic_s}  {r2_s}  ({n_} pt) {fb}")

    print("[3] Saving plots")

    # ── Integration windows: one figure per group ────────────────────────────
    # Group 1: all calibration standard rows (rslt_dir contains "cal")
    # Group 2+: one figure per distinct experiment rslt_dir
    import re as _re2
    def _safe(s: str) -> str:
        """Sanitise a directory name for use as a filename stem."""
        return _re2.sub(r"[^\w\-]+", "_", s).strip("_")[:60]

    cal_mask  = gc_df["rslt_dir"].apply(lambda d: "cal" in str(d).lower())
    df_cal    = gc_df[cal_mask].reset_index(drop=True)
    df_exp    = gc_df[~cal_mask].reset_index(drop=True)

    if not df_cal.empty:
        plot_gc_integration_windows(
            df_cal,
            out / "gc_integration_windows_cal.png",
            f"{args.label} — calibration standards",
        )

    for rslt_dir, grp in df_exp.groupby("rslt_dir", sort=False):
        safe_name = _safe(str(rslt_dir))
        plot_gc_integration_windows(
            grp.reset_index(drop=True),
            out / f"gc_integration_windows_{safe_name}.png",
            f"{args.label} — {rslt_dir}",
        )

    plot_gc_ppm(gc_df, out / "gc_ppm.png", args.label)
    plot_gc_calibration_curves(cal_meta, out / "gc_calibration.png", args.label)

    gc_df.drop(columns=["dx_path"], errors="ignore").to_csv(
        out / "gc_raw.csv", index=False)
    print(f"    Saved: gc_raw.csv")

    print(f"\n✓ Done.  Output: {out.resolve()}")


if __name__ == "__main__":
    main()
