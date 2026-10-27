"""
hplc_openlab_adaptive.py — Standalone HPLC parser + adaptive peak integrator
for Agilent OpenLAB .rslt folders (products: HCOOH, AcOH, EtOH, Propionate).

The bundled peak windows and compound map are illustrative method defaults;
review them against standards and blanks for every new method.

Usage (standalone):
    uv run python hplc_openlab_adaptive.py \
        --hplc-dirs "/path/to/sequence.rslt" \
        --label "Example sequence" \
        --out-dir output/hplc

Usage (as library):
    from hplc_openlab_adaptive import (
        collect_hplc_data,
        apply_hplc_calibration_factors,
        best_hplc_concentrations,
        integrate_hplc_peak_adaptive,
        plot_hplc_integration_windows,
        plot_hplc_calibration_curves,
    )
    hplc_df = collect_hplc_data(hplc_dirs)
    hplc_df, cal_meta = apply_hplc_calibration_factors(hplc_df, hplc_dirs)
    best = best_hplc_concentrations(hplc_df)
"""

from __future__ import annotations

import argparse
import math
import sys
import warnings
import zipfile
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

# ─────────────────────────────────────────────────────────────────────────────
#  Constants
# ─────────────────────────────────────────────────────────────────────────────

# Compound name mapping from Agilent OpenLAB XML names → short keys
HPLC_COMPOUND_MAP: dict[str, str] = {
    "Formic acid": "HCOOH",
    "Acetic acid": "AcOH",
    "Ethanol":     "EtOH",
    "propionate":  "Propionate",
}

# Colours for per-compound integration window highlights
HPLC_COLORS: dict[str, str] = {
    "HCOOH":      "teal",
    "AcOH":       "coral",
    "EtOH":       "purple",
    "Propionate": "orange",
}

# Half-window (minutes) used in the initial adaptive integration search
HPLC_HALF_WINDOW_MIN: float = 0.5

# Per-compound window multipliers.
#   search : Pass-1 window half-width = search × hw  (symmetric)
#   il, ir : Pass-2 integration window = [apex - il×hw, apex + ir×hw]
HPLC_WINDOW_PARAMS: dict[str, dict[str, float]] = {
    "HCOOH":      {"search": 2.0, "il": 1.0, "ir": 1.2, "left_extra": 0.0, "min_start_min": 25.0},
    "AcOH":       {"search": 3.2, "il": 1.0, "ir": 1.2, "left_extra": 1.5, "min_start_min": 27.5, "max_end_min": 30.0},
    "EtOH":       {"search": 8.0, "il": 1.0, "ir": 1.2, "left_extra": 0.0},
    "Propionate": {"search": 1.6, "il": 1.2, "ir": 1.44, "left_extra": 0.0, "max_end_min": 35.0},
}

# Fallback reference end-times (min) used when Agilent reports no peak data.
# Derived from manually confirmed injections on this HPLC RID method.
# center_est = end_time_min - 1.0 × hw_min
HPLC_REF_END_MIN: dict[str, float] = {
    "HCOOH":      26.9,   # centre ≈ 26.4 min
    "AcOH":       29.3,   # centre ≈ 28.8 min
    "EtOH":       43.2,   # centre ≈ 42.7 min
    "Propionate": 34.9,   # centre ≈ 34.4 min
}

# ─────────────────────────────────────────────────────────────────────────────
#  File I/O
# ─────────────────────────────────────────────────────────────────────────────

def read_hplc_rid_chromatogram(dx_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Read the RID1A chromatogram from an Agilent OpenLAB .dx file.

    .dx files are ZIP archives.  The signal parameters are in
    ``injection.acmd`` (XML); the binary trace is in ``<TraceId>.CH``.

    The binary trace is stored as float64 (``<f8``) values with a header
    that occupies everything before the last ``n * 8`` bytes.

    Returns
    -------
    t_min : np.ndarray
        Time axis in minutes, linearly spaced between TimeStart and TimeEnd.
    signal_nRIU : np.ndarray
        RID1A signal in nRIU (Nano Refractive Index Units), scaled by Slope.
    """
    with zipfile.ZipFile(dx_path) as zf:
        root = ET.parse(zf.open("injection.acmd")).getroot()
    ns = root.tag.split("}")[0].lstrip("{") if "}" in root.tag else ""
    T = lambda t: f"{{{ns}}}{t}" if ns else t
    info = root.find(T("InjectionInfo"))
    for sig in info.findall(f".//{T('Signal')}"):
        ch = sig.findtext(T("ChannelName")) or ""
        if ch != "RID1A":
            continue
        tid = sig.findtext(T("TraceId")) or ""
        n   = int(float(sig.findtext(T("NumberOfValues")) or 0))
        t0  = float(sig.findtext(T("TimeStart")) or 0)
        t1  = float(sig.findtext(T("TimeEnd"))   or 0)
        slo = float(sig.findtext(T("Slope"))     or 1.0)
        with zipfile.ZipFile(dx_path) as zf:
            raw = zf.open(f"{tid}.CH").read()
        data = np.frombuffer(raw[len(raw) - n * 8:], dtype="<f8")[:n] * slo
        t_min = np.linspace(t0, t1, n) / 60_000.0   # ms → min
        return t_min, data
    raise ValueError(f"RID1A channel not found in {dx_path}")


def parse_hplc_rx(rx_path: Path) -> dict:
    """
    Parse one Agilent OpenLAB .rx result file.

    .rx files are ZIP archives.  The result XML is at
    ``Base/InjectionACAML``.

    Timestamp is extracted from the **filename**, not the XML, because the
    XML metadata timestamp is unreliable across firmware versions.  Expected
    filename format: ``YYYY-MM-DD HH-MM-SS+HH-MM-NN.rx``
    (e.g. ``2026-05-14 11-52-54+08-00-01.rx``).

    Returns a dict with:
      - ``"run_datetime"`` : ISO-format string parsed from filename
      - one key per compound (short name): dict with
          ``conc_mM``, ``area_au``, ``end_min``

    Notes
    -----
    Agilent firmware quirk: ``<EndTime unit='s'>`` values are numerically
    in **minutes**, not seconds.  We treat 's'-labelled values as minutes
    and prefer 'min'-labelled values if present.
    """
    stem = rx_path.stem
    try:
        date_part = stem[:10]
        time_part = stem[11:19].replace("-", ":")
        tz_part   = stem[19:25].replace("-", ":")
        run_dt_str = f"{date_part}T{time_part}{tz_part}"
    except Exception:
        run_dt_str = ""

    with zipfile.ZipFile(rx_path) as zf:
        root = ET.parse(zf.open("Base/InjectionACAML")).getroot()
    ns = root.tag.split("}")[0].lstrip("{") if "}" in root.tag else ""
    T = lambda t: f"{{{ns}}}{t}" if ns else t

    result: dict = {"run_datetime": run_dt_str}
    for ic in root.findall(f".//{T('InjectionCompound')}"):
        raw_name   = (ic.findtext(T("CompoundName")) or "").strip()
        short_name = HPLC_COMPOUND_MAP.get(raw_name, raw_name)
        if not short_name:
            continue

        # Concentration (mM) — try Amount then Concentration attribute
        conc = math.nan
        for el in [ic.find(f".//{T('Amount')}"), ic.find(f".//{T('Concentration')}")]:
            if el is None:
                continue
            val = (el.get("val") or "").strip()
            if val and val.lower() != "nan":
                try:
                    conc = float(val); break
                except ValueError:
                    pass

        # Peak area (Agilent internal units)
        area = math.nan
        area_el = ic.find(f".//{T('Area')}")
        if area_el is not None:
            try:
                area = float(area_el.get("val") or "nan")
            except ValueError:
                pass

        # Peak end time — unit='s' values are numerically in minutes
        end_min = math.nan
        end_s_fallback = math.nan
        for e in ic.findall(f".//{T('EndTime')}"):
            try:
                v = float(e.get("val", "nan"))
            except ValueError:
                continue
            if e.get("unit") == "min" and math.isfinite(v):
                end_min = v
                break
            elif e.get("unit") == "s" and math.isfinite(v):
                end_s_fallback = v
        if not math.isfinite(end_min) and math.isfinite(end_s_fallback):
            end_min = end_s_fallback   # use 's'-labelled value as minutes

        prev = result.get(short_name, {})
        prev_conc = prev.get("conc_mM", math.nan) if isinstance(prev, dict) else prev
        if math.isnan(prev_conc) or (math.isfinite(conc) and conc > prev_conc):
            result[short_name] = {
                "conc_mM": conc,
                "area_au": area,
                "end_min": end_min,
            }
    return result


def collect_hplc_data(hplc_dirs: list[Path]) -> pd.DataFrame:
    """
    Scan one or more .rslt directories and return one DataFrame row per
    injection (sorted by run_datetime).

    Rerun-analysis files (``-r001``, ``-r002``, …) are included because
    for short runs they may carry the only usable data.

    Injections where Agilent found no peaks (``InjectionCompound`` absent)
    still get a row, populated with ``NaN`` concentrations and the fallback
    reference end-times from ``HPLC_REF_END_MIN`` so raw chromatograms can
    still be drawn and re-integrated.

    Columns per compound (e.g. ``AcOH``):
      - ``AcOH``           : Agilent-reported concentration (mM) or NaN
      - ``AcOH_area``      : Agilent-reported area (native units) or NaN
      - ``AcOH_end_min``   : peak end-time used as integration anchor (min)
    """
    rows = []
    for rslt in hplc_dirs:
        for rx in sorted(rslt.glob("*.rx")):
            if "-r0" in rx.name:
                continue
            try:
                data = parse_hplc_rx(rx)
                flat = {
                    "run_datetime": data.pop("run_datetime"),
                    "rx_file":  rx.name,
                    "rslt":     rslt.name,
                    "dx_path":  str(rslt / rx.stem) + ".dx",
                }
                compounds_found = False
                for cname, cdata in data.items():
                    if isinstance(cdata, dict):
                        flat[cname]              = cdata.get("conc_mM", math.nan)
                        flat[f"{cname}_area"]    = cdata.get("area_au", math.nan)
                        flat[f"{cname}_end_min"] = cdata.get("end_min", math.nan)
                        compounds_found = True
                    else:
                        flat[cname] = cdata
                if not compounds_found:
                    for cname, ref_end in HPLC_REF_END_MIN.items():
                        flat.setdefault(cname, math.nan)
                        flat.setdefault(f"{cname}_area", math.nan)
                        flat[f"{cname}_end_min"] = ref_end
                    flat["_no_agilent_data"] = True
                rows.append(flat)
            except Exception as e:
                warnings.warn(f"HPLC parse failed {rx}: {e}")
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["run_datetime"] = pd.to_datetime(df["run_datetime"], utc=True, errors="coerce")
    df = df.sort_values("run_datetime").reset_index(drop=True)
    return df


# ─────────────────────────────────────────────────────────────────────────────
#  Adaptive peak integration
# ─────────────────────────────────────────────────────────────────────────────

def integrate_hplc_peak_adaptive(
    t_min: np.ndarray,
    sig: np.ndarray,
    end_time_min: float,
    hw_min: float = HPLC_HALF_WINDOW_MIN,
    search_mult: float = 2.0,
    integ_left_mult: float = 1.0,
    integ_right_mult: float = 1.2,
    search_left_extra: float = 0.0,
    min_start_min: float | None = None,
    max_end_min: float | None = None,
) -> dict:
    """
    Two-pass adaptive integration of one HPLC peak.

    The Agilent-reported ``EndTime`` is used as the right anchor for the
    initial search window.  Centre estimate = ``end_time_min - 1.0 × hw``.

    Pass 1 — apex search
        Window: ``[centre - search_mult×hw,  centre + search_mult×hw]``
        Find the highest signal point → apex_t

    Pass 2 — integration
        Window: ``[apex_t - integ_left_mult×hw,  apex_t + integ_right_mult×hw]``
        Sloped baseline: straight line connecting sig[left_edge] → sig[right_edge]
        (valley-to-valley method; matches Agilent's standard baseline)
        Area = trapezoid integral of (signal − baseline), clipped ≥ 0

    Per-compound defaults are stored in ``HPLC_WINDOW_PARAMS``:
      AcOH           : search_mult = 3.2 (×0.8); min_start_min = 27.5; max_end_min = 30.0
      HCOOH          : min_start_min = 25.0
      Propionate     : search_mult = 1.6 (×0.8); max_end_min = 35.0

    Parameters
    ----------
    t_min : array of time values in minutes
    sig : array of RID signal values (nRIU)
    end_time_min : Agilent EndTime used as right anchor (minutes)
    hw_min : half-window in minutes (default 0.5 min)
    search_mult : symmetric search half-width multiplier (default 2.0)
    integ_left_mult : left integration half-width multiplier (default 1.0)
    integ_right_mult : right integration half-width multiplier (default 1.2)
    min_start_min : if set, force area = 0 when left integration edge is below this (min)
    max_end_min : if set, force area = 0 when right integration edge exceeds this (min)

    Returns
    -------
    dict with keys:
      ``area``       — integrated area in nRIU·min (NaN if integration failed)
      ``apex_min``   — time of peak apex (min)
      ``start_min``  — left edge of integration window
      ``end_min``    — right edge of integration window
      ``t_window``   — time array of integration window
      ``sig_window`` — signal array of integration window
      ``baseline``   — sloped baseline array (same length as t_window)
    """
    if not math.isfinite(end_time_min):
        return {"area": math.nan}

    search_left  = search_mult * hw_min + search_left_extra
    search_right = search_mult * hw_min
    integ_left   = integ_left_mult * hw_min
    integ_right  = integ_right_mult * hw_min
    center_est   = end_time_min - 1.0 * hw_min

    mask_s = (t_min >= center_est - search_left) & (t_min <= center_est + search_right)
    if mask_s.sum() < 2:
        return {"area": math.nan}
    apex_t = t_min[mask_s][np.argmax(sig[mask_s])]

    mask_i = (t_min >= apex_t - integ_left) & (t_min <= apex_t + integ_right)
    tw = t_min[mask_i]
    sw = sig[mask_i]
    if len(tw) < 2:
        return {"area": math.nan, "apex_min": float(apex_t)}

    baseline = np.linspace(float(sw[0]), float(sw[-1]), len(tw))
    start_min = float(tw[0])
    end_min   = float(tw[-1])
    reject = (
        (min_start_min is not None and start_min < min_start_min)
        or (max_end_min is not None and end_min > max_end_min)
    )
    area = 0.0 if reject else max(float(np.trapezoid(sw - baseline, tw)), 0.0)

    return {
        "area":       area,
        "apex_min":   float(apex_t),
        "start_min":  start_min,
        "end_min":    end_min,
        "t_window":   tw,
        "sig_window": sw,
        "baseline":   baseline,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  Unit-consistent calibration (3-pass)
# ─────────────────────────────────────────────────────────────────────────────

def apply_hplc_calibration_factors(
    hplc_df: pd.DataFrame,
    hplc_dirs: list[Path],
) -> tuple[pd.DataFrame, dict]:
    """
    Build a unit-consistent calibration (mM / nRIU·min) from Agilent-
    quantified injections and apply it to rows where Agilent returned NaN.

    Why this is necessary
    ---------------------
    Agilent's internal ``Area`` (``area_au``) uses proprietary scaling that
    is **not** in nRIU·min.  If you mix Agilent areas with manually
    re-integrated areas you get a unit mismatch and wildly wrong
    concentrations.  This function avoids the mismatch by:

    Pass 1 — re-integrate **every** injection from its raw .dx trace so all
             areas are in consistent ``nRIU·min`` units.

    Pass 2 — for injections where Agilent reported a valid concentration,
             compute the calibration factor:
             ``factor [mM / nRIU·min] = mean(Agilent_conc / our_area)``

    Pass 3 — for injections where Agilent returned NaN, fill the
             concentration using:
             ``conc_mM = our_area × factor``

    Parameters
    ----------
    hplc_df : DataFrame returned by ``collect_hplc_data``
    hplc_dirs : same list passed to ``collect_hplc_data`` (used for dx lookup)

    Returns
    -------
    (modified_df, hplc_cal_meta)

    hplc_cal_meta[compound] has keys:
      ``factor``       — mM / nRIU·min  (NaN if no calibration data)
      ``areas``        — list of our nRIU·min areas from calibration rows
      ``concs``        — list of Agilent concentrations from calibration rows
      ``source_files`` — list of .rx filenames used for calibration
    """
    compounds = list(HPLC_REF_END_MIN.keys())

    # ── Pass 1: re-integrate every injection ─────────────────────────────────
    our_areas: dict[int, dict[str, float]] = {}
    for idx, row in hplc_df.iterrows():
        dx_p = Path(str(row.get("dx_path", "")))
        if not dx_p.exists():
            continue
        try:
            t_min, sig = read_hplc_rid_chromatogram(dx_p)
        except Exception:
            continue
        our_areas[idx] = {}
        for c in compounds:
            end_t = row.get(f"{c}_end_min", math.nan)
            if not math.isfinite(end_t):
                end_t = HPLC_REF_END_MIN.get(c, math.nan)
            wp = HPLC_WINDOW_PARAMS.get(c, {"search": 2.0, "il": 1.0, "ir": 1.2, "left_extra": 0.0})
            pk = integrate_hplc_peak_adaptive(t_min, sig, end_t,
                                              search_mult=wp["search"],
                                              integ_left_mult=wp["il"],
                                              integ_right_mult=wp["ir"],
                                              search_left_extra=wp.get("left_extra", 0.0),
                                              min_start_min=wp.get("min_start_min"),
                                              max_end_min=wp.get("max_end_min"))
            our_areas[idx][c] = pk.get("area", math.nan)

    # ── Pass 2: derive calibration from STANDARD rows only ───────────────────
    # A row is a calibration standard when its rslt folder name contains "cal"
    # (e.g. "0_calibration_...", "Cal1.rslt").  Experiment rows are excluded
    # from the fit; their concentrations are computed in Pass 3 from the fit.
    #
    # max_agilent is computed from standards only so that experiment-sample
    # concentrations never inflate the noise threshold.
    max_agilent: dict[str, float] = {}
    for c in compounds:
        if c in hplc_df.columns:
            std_mask = hplc_df["rslt"].apply(lambda d: "cal" in str(d).lower())
            pos = hplc_df.loc[std_mask, c].dropna() if c in hplc_df.columns else pd.Series(dtype=float)
            pos = pos[pos > 0]
            max_agilent[c] = float(pos.max()) if not pos.empty else 0.0
        else:
            max_agilent[c] = 0.0

    std_data:  dict[str, dict] = {c: {"areas": [], "concs": [], "files": []} for c in compounds}

    for idx, row in hplc_df.iterrows():
        rslt = str(row.get("rslt", ""))
        if "cal" not in rslt.lower():
            continue  # skip experiment rows in Pass 2
        for c in compounds:
            agilent_conc = row.get(c, math.nan)
            our_area     = our_areas.get(idx, {}).get(c, math.nan)
            min_thresh   = max(0.5, max_agilent.get(c, 0.0) * 0.05)
            if (math.isfinite(agilent_conc) and agilent_conc >= min_thresh
                    and math.isfinite(our_area) and our_area > 0):
                std_data[c]["areas"].append(float(our_area))
                std_data[c]["concs"].append(float(agilent_conc))
                std_data[c]["files"].append(str(row.get("rx_file", "unknown")))

    fits: dict[str, dict] = {}
    for c in compounds:
        areas = std_data[c]["areas"]
        concs = std_data[c]["concs"]
        if len(areas) >= 2:
            coeffs = np.polyfit(areas, concs, 1)
            slope, intercept = float(coeffs[0]), float(coeffs[1])
            concs_pred = slope * np.array(areas) + intercept
            ss_res = float(np.sum((np.array(concs) - concs_pred) ** 2))
            ss_tot = float(np.sum((np.array(concs) - np.mean(concs)) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
        elif len(areas) == 1:
            slope = float(concs[0] / areas[0]) if areas[0] != 0 else math.nan
            intercept = 0.0
            r2 = math.nan
        else:
            slope = intercept = r2 = math.nan
        fits[c] = {"slope": slope, "intercept": intercept, "r2": r2}

    # ── Pass 3: calculate concentrations for ALL non-standard rows ────────────
    # Experiment rows: always derive concentration from the calibration fit so
    # that reported mM is consistent with the calibration line.
    # Collect (our_area, fit_conc) for the calibration plot — these points lie
    # exactly on the line by construction.
    samp_plot: dict[str, dict] = {c: {"areas": [], "concs": [], "files": []} for c in compounds}

    for idx, row in hplc_df.iterrows():
        rslt = str(row.get("rslt", ""))
        if "cal" in rslt.lower():
            continue  # standards keep their Agilent values
        for c in compounds:
            our_area = our_areas.get(idx, {}).get(c, math.nan)
            fit = fits.get(c, {})
            sl = fit.get("slope", math.nan)
            ic = fit.get("intercept", math.nan)
            if not (math.isfinite(our_area) and our_area > 0
                    and math.isfinite(sl) and math.isfinite(ic)):
                continue
            fit_conc = max(sl * our_area + ic, 0.0)
            hplc_df.at[idx, c] = fit_conc
            hplc_df.at[idx, f"{c}_area"] = our_area
            hplc_df.at[idx, "_reintegrated"] = True
            samp_plot[c]["areas"].append(float(our_area))
            samp_plot[c]["concs"].append(fit_conc)
            samp_plot[c]["files"].append(str(row.get("rx_file", "unknown")))

    hplc_cal_meta: dict[str, dict] = {}
    for c in compounds:
        fit = fits[c]
        hplc_cal_meta[c] = {
            "slope":        fit["slope"],
            "intercept":    fit["intercept"],
            "r2":           fit["r2"],
            "areas":        std_data[c]["areas"],
            "concs":        std_data[c]["concs"],
            "source_files": std_data[c]["files"],
            # Experiment sample points — on the line by construction
            "samp_areas":        samp_plot[c]["areas"],
            "samp_concs":        samp_plot[c]["concs"],
            "samp_source_files": samp_plot[c]["files"],
        }

    return hplc_df, hplc_cal_meta


def best_hplc_concentrations(hplc_df: pd.DataFrame) -> dict[str, float]:
    """Return the maximum (most-accumulated) concentration for each compound."""
    best = {}
    for c in list(HPLC_REF_END_MIN.keys()):
        if c not in hplc_df.columns:
            best[c] = math.nan
            continue
        valid = hplc_df[c].dropna()
        valid = valid[valid > 0]
        best[c] = float(valid.max()) if not valid.empty else math.nan
    return best


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting
# ─────────────────────────────────────────────────────────────────────────────

def plot_hplc_integration_windows(
    hplc_df: pd.DataFrame,
    out_path: Path,
    label: str,
    compounds: list[str] | None = None,
    extra_ref_end_min: dict[str, float] | None = None,
    extra_window_params: dict[str, dict] | None = None,
    extra_colors: dict[str, str] | None = None,
) -> None:
    """
    One row per injection: RID1A chromatogram with coloured integration
    windows and annotated concentrations.

    Layout rules
    ------------
    - All rows share the same x-axis range (global min/max across all
      integration windows + padding), enabling direct visual comparison.
    - y-axis auto-scales to the 2nd–98th percentile of the visible window,
      expanded to include all integration window peaks.
    - Each compound window is highlighted with ``axvspan`` (alpha=0.25)
      in its assigned colour; the apex is marked with a dashed vertical line.
    - Annotation shows compound name and concentration above the apex.
    - Source tag in subplot title: ``[Agilent]`` = Agilent-quantified,
      ``[re-integrated]`` = our calibration applied, ``[ref windows only]``
      = no Agilent data at all.
    """
    _ref_end  = {**HPLC_REF_END_MIN,  **(extra_ref_end_min  or {})}
    _win_par  = {**HPLC_WINDOW_PARAMS, **(extra_window_params or {})}
    _colors   = {**HPLC_COLORS,        **(extra_colors       or {})}
    compounds = compounds if compounds is not None else list(HPLC_REF_END_MIN.keys())
    valid_rows = [row for _, row in hplc_df.iterrows()
                  if Path(str(row.get("dx_path", ""))).exists()]
    if not valid_rows:
        print("    ⚠ No readable HPLC .dx files — skipping integration window plot.")
        return

    n_inj = len(valid_rows)
    fig, axes = plt.subplots(n_inj, 1, figsize=(14, 3.5 * n_inj), squeeze=False)

    # Pass 1: global x range (use per-compound search_mult for accurate bounds)
    all_x0s, all_x1s = [], []
    for row in valid_rows:
        for comp in compounds:
            end_t = row.get(f"{comp}_end_min", math.nan)
            if not math.isfinite(end_t):
                end_t = _ref_end.get(comp, math.nan)
            if not math.isfinite(end_t):
                continue
            hw = HPLC_HALF_WINDOW_MIN
            wp = _win_par.get(comp, {"search": 2.0, "left_extra": 0.0})
            sm = wp["search"]
            le = wp.get("left_extra", 0.0)
            all_x0s.append(end_t - 1.0 * hw - sm * hw - le - 1.5)
            all_x1s.append(end_t + sm * hw + 1.5)

    global_x_lo = max(0.0, min(all_x0s)) if all_x0s else 15.0
    global_x_hi = max(all_x1s) if all_x1s else 55.0

    # Pass 2: draw
    for row_idx, row in enumerate(valid_rows):
        ax = axes[row_idx][0]
        dx_p = Path(str(row["dx_path"]))
        try:
            t_min, sig = read_hplc_rid_chromatogram(dx_p)
        except Exception as e:
            ax.text(0.5, 0.5, f"read error: {e}", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9)
            continue

        run_ts = row.get("run_datetime", "")
        try:
            run_ts = run_ts.strftime("%Y-%m-%d %H:%M")
        except Exception:
            run_ts = str(row.get("rx_file", ""))
        no_agilent    = bool(row.get("_no_agilent_data", False))
        reintegrated  = bool(row.get("_reintegrated", False))
        source_tag = (" [ref windows only]" if no_agilent
                      else (" [re-integrated]" if reintegrated else " [Agilent]"))
        ax.set_title(f"HPLC {run_ts} — integration windows{source_tag}", fontsize=9)

        legend_patches = []
        window_spans   = []
        for comp in compounds:
            end_col = f"{comp}_end_min"
            end_t = float(row[end_col]) if (end_col in row and math.isfinite(row.get(end_col, math.nan))) \
                    else _ref_end.get(comp, math.nan)
            if not math.isfinite(end_t):
                continue
            wp = _win_par.get(comp, {"search": 2.0, "il": 1.0, "ir": 1.2, "left_extra": 0.0})
            pk = integrate_hplc_peak_adaptive(t_min, sig, end_t,
                                              search_mult=wp["search"],
                                              integ_left_mult=wp["il"],
                                              integ_right_mult=wp["ir"],
                                              search_left_extra=wp.get("left_extra", 0.0),
                                              min_start_min=wp.get("min_start_min"),
                                              max_end_min=wp.get("max_end_min"))
            if not math.isfinite(pk.get("area", math.nan)):
                continue
            color    = _colors.get(comp, "gray")
            conc     = row.get(comp, math.nan)
            conc_txt = f"{conc:.3f} mM" if math.isfinite(conc) else "N/D"
            window_spans.append(
                (pk["start_min"], pk["end_min"], pk["apex_min"], comp, conc_txt, color, pk))
            legend_patches.append(
                mpatches.Patch(color=color, alpha=0.55, label=f"{comp}: {conc_txt}"))

        ax.set_xlim(global_x_lo, global_x_hi)
        vis = (t_min >= global_x_lo) & (t_min <= global_x_hi)
        if vis.sum() > 1:
            y_lo = float(np.nanpercentile(sig[vis], 2))
            y_hi = float(np.nanpercentile(sig[vis], 98))
            for x0_w, x1_w, *_ in window_spans:
                wm = (t_min >= x0_w) & (t_min <= x1_w)
                if wm.sum() > 0:
                    y_lo = min(y_lo, float(np.nanmin(sig[wm])))
                    y_hi = max(y_hi, float(np.nanmax(sig[wm])))
            pad = max((y_hi - y_lo) * 0.25, 5.0)
            ax.set_ylim(y_lo - pad * 0.2, y_hi + pad)

        ax.plot(t_min, sig, color="black", lw=0.6, zorder=2)
        y_ann = ax.get_ylim()[1]
        for x0, x1, apex_t, comp, conc_txt, color, _pk in window_spans:
            ax.axvspan(x0, x1, alpha=0.25, color=color, zorder=1)
            ax.axvline(apex_t, color=color, lw=1.2, ls="--", zorder=3)
            ax.text(apex_t, y_ann * 0.97,
                    f" {comp}\n {conc_txt}", fontsize=7, color=color,
                    va="top", clip_on=True)

        ax.set_xlabel("Time (min)", fontsize=8)
        ax.set_ylabel("RID (nRIU)", fontsize=8)
        ax.grid(False)
        if legend_patches:
            leg = ax.legend(handles=legend_patches, fontsize=8, loc="best")
            leg.get_frame().set_visible(False)
        ax.tick_params(labelsize=7)

    fig.suptitle(f"{label} — HPLC Integration Windows (RID1A)", fontsize=10)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved: {out_path.name}")


def _pchip_plot(ax, x, y, marker, ls, color, label_str, ms=4, lw=1.2, alpha=0.8):
    """Scatter measured points and connect with a PCHIP smooth curve."""
    xa = np.array(list(x))
    ya = np.asarray(y, dtype=float)
    valid = np.isfinite(ya)
    ax.scatter(xa[valid], ya[valid], marker=marker, color=color,
               s=ms ** 2, alpha=alpha, zorder=4, label=label_str)
    xv, yv = xa[valid], ya[valid]
    if len(xv) >= 2:
        xi = np.linspace(xv[0], xv[-1], 300)
        ax.plot(xi, PchipInterpolator(xv, yv)(xi), ls=ls, color=color,
                lw=lw, alpha=alpha)


def plot_hplc_liquid(
    hplc_df: pd.DataFrame,
    out_path: Path,
    label: str,
) -> None:
    """
    Line plot of concentration (mM) vs. injection number (1-based).
    NaN values are rendered as gaps (not filled with zero).
    """
    compounds = [c for c in HPLC_REF_END_MIN if c in hplc_df.columns]
    fig, ax = plt.subplots(figsize=(8, 4))
    x = range(1, len(hplc_df) + 1)
    for comp in compounds:
        color = HPLC_COLORS.get(comp, "gray")
        vals = hplc_df[comp].where(hplc_df[comp] > 0).to_numpy(dtype=float)
        _pchip_plot(ax, x, vals, "o", "-", color, comp)
    ax.set_xticks(list(x))
    ax.set_xlabel("HPLC injection #", fontsize=9)
    ax.set_ylabel("Concentration (mM)", fontsize=9)
    ax.set_title(f"{label} — Liquid products (HPLC RID)", fontsize=10)
    ax.grid(False)
    leg = ax.legend(fontsize=8, loc="best")
    leg.get_frame().set_visible(False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved: {out_path.name}")


def plot_hplc_calibration_curves(
    hplc_cal_meta: dict,
    out_path: Path,
    label: str,
) -> None:
    """
    One subplot per compound: calibration scatter + proper linear regression
    (with intercept, via np.polyfit).

    Layout
    ------
    - x-axis: our re-integrated area (nRIU·min)
    - y-axis: Agilent-reported concentration (mM)
    - Blue scatter: calibration points, labelled with source filename
    - Red dashed line: best-fit line extended to zero
    - Title: compound, equation  y = slope·x + intercept, R², n points
    - Log scale applied on both axes when concentration range > 50×
    """
    compounds = [c for c in HPLC_REF_END_MIN
                 if hplc_cal_meta.get(c, {}).get("areas")]
    if not compounds:
        return
    n = len(compounds)
    # Layout rule: rows = floor(sqrt(n)), cols = floor(sqrt(n)) + 1
    import math as _math
    n_rows = max(1, int(_math.sqrt(n)))
    n_cols = n_rows if n_rows * n_rows == n else n_rows + 1
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5.0 * n_cols, 4.5 * n_rows), squeeze=False)

    for i, comp in enumerate(compounds):
        r, c = divmod(i, n_cols)
        ax   = axes[r][c]
        meta = hplc_cal_meta[comp]
        areas = np.array(meta["areas"], dtype=float)
        concs = np.array(meta["concs"],  dtype=float)
        slope     = meta.get("slope",     math.nan)
        intercept = meta.get("intercept", math.nan)
        r2        = meta.get("r2",        math.nan)
        color     = HPLC_COLORS.get(comp, "steelblue")

        use_log = (len(concs) >= 2
                   and concs.min() > 0
                   and concs.max() / concs.min() > 50)

        # ── Standards: filled circles — used to build the fit ──────────────
        if len(areas):
            ax.scatter(areas, concs, color=color, zorder=4, s=60, alpha=0.6,
                       marker="o", edgecolors="white", linewidths=0.5,
                       label="standard")

        # ── Regression line (standards range only) ──────────────────────
        if math.isfinite(slope) and math.isfinite(intercept) and len(areas) >= 1:
            a_span = float(areas.max() - areas.min()) if len(areas) > 1 else float(areas[0])
            a_lo = max(0.0, float(areas.min()) - a_span * 0.1)
            a_hi = float(areas.max()) + a_span * 0.1
            a_fit = np.linspace(a_lo, a_hi, 300)
            ax.plot(a_fit, slope * a_fit + intercept, color="gray", lw=1.8,
                    ls="--", alpha=0.7, zorder=3, label="linear fit")
            ax.plot([0, a_lo], [intercept, slope * a_lo + intercept],
                    color="gray", lw=0.8, ls=":", alpha=0.35, zorder=2)

        if use_log:
            ax.set_xscale("log")
            ax.set_yscale("log")

        ax.set_xlabel("Peak area (nRIU·min)", fontsize=9)
        ax.set_ylabel("Concentration (mM)", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(False)
        leg = ax.legend(fontsize=7, loc="best")
        leg.get_frame().set_visible(False)

        sl_str = f"{slope:.4e}" if math.isfinite(slope) else "N/A"
        ic_str = (f"{intercept:+.4e}" if math.isfinite(intercept) else "")
        r2_str = f"R² = {r2:.4f}" if math.isfinite(r2) else "R² = N/A"
        ax.set_title(
            f"{comp}   n={len(areas)} std\n"
            f"y={sl_str}·x {ic_str}\n{r2_str}",
            fontsize=8)

    # Hide unused axes
    for idx in range(len(compounds), n_rows * n_cols):
        r, c = divmod(idx, n_cols)
        axes[r][c].set_visible(False)

    fig.suptitle(
        f"{label} — HPLC Calibration (linear regression, nRIU·min vs. mM)",
        fontsize=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved: {out_path.name}")


# ─────────────────────────────────────────────────────────────────────────────
#  CLI entry-point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse Agilent OpenLAB HPLC .rslt folders and plot "
                    "integration windows + calibration curves.")
    parser.add_argument("--hplc-dirs", nargs="+", required=True, type=Path,
                        help=".rslt directory/directories to process")
    parser.add_argument("--label", default="HPLC", help="Title label for plots")
    parser.add_argument("--out-dir", default="output/hplc", type=Path,
                        help="Output directory")
    args = parser.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    print(f"[1] Collecting HPLC injections from {len(args.hplc_dirs)} .rslt dir(s)")
    hplc_df = collect_hplc_data(args.hplc_dirs)
    if hplc_df.empty:
        print("    No injections found.  Exiting.")
        sys.exit(1)
    print(f"    {len(hplc_df)} injection(s) found")
    print(hplc_df[["run_datetime", "rx_file"] +
                  [c for c in HPLC_REF_END_MIN if c in hplc_df.columns]].to_string())

    print("[2] Applying unit-consistent calibration (3-pass)")
    hplc_df, cal_meta = apply_hplc_calibration_factors(hplc_df, args.hplc_dirs)
    print("    Calibration fit (linear regression, nRIU·min → mM):")
    for c, m in cal_meta.items():
        sl = m.get("slope", math.nan)
        ic = m.get("intercept", math.nan)
        r2 = m.get("r2", math.nan)
        n  = len(m["areas"])
        sl_s = f"{sl:.4e}" if math.isfinite(sl) else "N/A"
        ic_s = f"{ic:+.4e}" if math.isfinite(ic) else ""
        r2_s = f"R²={r2:.4f}" if math.isfinite(r2) else "R²=N/A"
        print(f"      {c:12s}  y={sl_s}·x{ic_s}  {r2_s}  ({n} pt)")

    print("[3] Best concentrations:")
    best = best_hplc_concentrations(hplc_df)
    for c, v in best.items():
        print(f"      {c:12s}  {v:.4f} mM" if math.isfinite(v) else f"      {c:12s}  N/A")

    import re as _re2, shutil as _sh
    _s = out / "scripts"; _s.mkdir(exist_ok=True)
    _sh.copy2(__file__, _s / Path(__file__).name)

    print("[4] Saving plots")

    # ── Integration windows: one figure per group ────────────────────────────
    # Group 1: all calibration standard rows (rslt folder contains "cal")
    # Group 2+: one figure per distinct experiment rslt folder
    def _safe(s: str) -> str:
        return _re2.sub(r"[^\w\-]+", "_", s).strip("_")[:60]

    cal_mask = hplc_df["rslt"].apply(lambda d: "cal" in str(d).lower())
    df_cal   = hplc_df[cal_mask].reset_index(drop=True)
    df_exp   = hplc_df[~cal_mask].reset_index(drop=True)

    if not df_cal.empty:
        plot_hplc_integration_windows(
            df_cal,
            out / "hplc_integration_windows_cal.png",
            f"{args.label} — calibration standards",
        )
    for rslt_dir, grp in df_exp.groupby("rslt", sort=False):
        safe_name = _safe(str(rslt_dir))
        plot_hplc_integration_windows(
            grp.reset_index(drop=True),
            out / f"hplc_integration_windows_{safe_name}.png",
            f"{args.label} — {rslt_dir}",
        )

    plot_hplc_liquid(hplc_df, out / "hplc_liquid.png", args.label)
    plot_hplc_calibration_curves(cal_meta, out / "hplc_calibration.png", args.label)

    hplc_df.to_csv(out / "hplc_raw.csv", index=False)
    print(f"    Saved: hplc_raw.csv")
    print(f"\n✓ Done.  Output: {out.resolve()}")


if __name__ == "__main__":
    main()
