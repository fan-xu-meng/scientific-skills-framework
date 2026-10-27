"""
bruker_raw4_analyzer.py  –  Parse Bruker D8 binary .raw (RAW4.00) XRD files
and produce stacked / overlay plots. No external XRD library required.

Usage (CLI):
    python bruker_raw4_analyzer.py --data-dir <folder> --output-dir <out>

Usage (import):
    from bruker_raw4_analyzer import parse_raw4, plot_stacked, plot_overlay

Bruker RAW v4 (RAW4.00) binary format (reverse-engineered from sample files):
  Magic bytes : 'RAW4.00\\0'
  start_2theta: float64 LE found by byte-scan in 0x02A0–0x02F5 (value 3–15°)
  step_size   : float64 LE at start_offset + 8 bytes (value 0.001–0.15°)
  Count data  : longest consecutive float32 LE run (30–500 000 cts)
                starting near 0x47E, at whichever 4-byte alignment fits.
  Typical range: 5–90° 2θ, Cu Kα, step 0.0237° (3588 pts) or 0.0395° (2153 pts)
"""

import struct
import re
import glob
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── RAW4 Parser ─────────────────────────────────────────────────────────────

def _find_scan_params(data: bytes):
    """Find (start_2theta, step_size) by scanning for adjacent float64 pair.
    The pair is 8 bytes apart; scan every byte to be robust to header size variation."""
    for i in range(0x02A0, 0x02F5, 1):
        if i + 16 > len(data):
            break
        v1 = struct.unpack_from("<d", data, i)[0]
        v2 = struct.unpack_from("<d", data, i + 8)[0]
        if 3.0 < v1 < 15.0 and 0.001 < v2 < 0.15:
            return v1, v2
    return None, None


def _find_count_data(data: bytes):
    """
    Find the float32 count array: longest consecutive run of values in (30, 500000),
    trying all 4 byte alignments starting near offset 0x47E.
    Returns (offset, counts_array).
    """
    best_off, best_n = 0, 0
    for align in range(4):
        i = 0x47E + align
        curr_n = 0
        curr_start = i
        while i < len(data) - 4:
            v = struct.unpack_from("<f", data, i)[0]
            if 30 < v < 500_000 and v == v:
                if curr_n == 0:
                    curr_start = i
                curr_n += 1
                if curr_n > best_n:
                    best_n = curr_n
                    best_off = curr_start
            else:
                curr_n = 0
            i += 4
    if best_n < 100:
        return None, None
    counts = np.array(
        [struct.unpack_from("<f", data, best_off + k * 4)[0] for k in range(best_n)],
        dtype=float,
    )
    return best_off, counts


def parse_raw4(filepath: str):
    """
    Parse a Bruker RAW v4 file.
    Returns dict with keys: name, two_theta (ndarray), counts (ndarray).
    """
    data = Path(filepath).read_bytes()
    if not data.startswith(b"RAW4"):
        raise ValueError(f"Not a RAW4 file: {filepath}")

    start, step = _find_scan_params(data)
    _, counts = _find_count_data(data)

    if counts is None:
        raise ValueError(f"Could not find count data in {filepath}")
    if start is None:
        raise ValueError(f"Could not find scan parameters in {filepath}")

    two_theta = start + np.arange(len(counts)) * step
    name = Path(filepath).stem
    return {"name": name, "path": filepath, "two_theta": two_theta, "counts": counts}


# ── Plotting ─────────────────────────────────────────────────────────────────

# Known XRD reference peaks (Cu Kα, JCPDS / literature).
# Each entry: list of (2θ°, relative_intensity 0–100).
REFERENCES = {
    "Pt (FCC)": {
        "peaks": [(39.76, 100), (46.24, 52), (67.46, 32),
                  (81.27, 36), (85.73, 12)],
        "color": "#1f77b4",
    },
    "Ni(OH)₂": {
        "peaks": [(19.3, 90), (33.1, 100), (38.6, 55),
                  (52.1, 30), (59.0, 60), (63.1, 20), (70.2, 15)],
        "color": "#d62728",
    },
    "Co(OH)₂": {
        "peaks": [(19.1, 100), (32.5, 60), (37.9, 40),
                  (51.5, 20), (58.0, 30)],
        "color": "#2ca02c",
    },
    "Carbon": {
        "peaks": [(26.3, 100)],
        "color": "#7F8C8D",
    },
}

_PALETTE = [
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
    "#9467bd", "#8c564b", "#e377c2", "#17becf",
]


def _pub_style():
    plt.rcParams.update({
        "font.family":      "DejaVu Sans",
        "font.size":        16,
        "axes.labelsize":   18,
        "axes.titlesize":   17,
        "xtick.labelsize":  15,
        "ytick.labelsize":  15,
        "axes.linewidth":   1.4,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.direction":  "in",
        "ytick.direction":  "in",
        "figure.dpi":       150,
    })


def _relevant_refs(sample_name: str) -> dict:
    """Return only the REFERENCES relevant to this sample based on its name."""
    name = sample_name.lower()
    keep = {}
    if "pt" in name:
        keep["Pt (FCC)"] = REFERENCES["Pt (FCC)"]
    if "ni(oh)" in name or "nioh" in name:
        keep["Ni(OH)₂"] = REFERENCES["Ni(OH)₂"]
    if "co(oh)" in name or "cooh" in name:
        keep["Co(OH)₂"] = REFERENCES["Co(OH)₂"]
    # Fallback: show all if nothing matched
    return keep if keep else REFERENCES


def _estimate_background(counts, window_frac=0.08):
    """Morphological opening (min→max filter) background estimate for XRD."""
    from scipy.ndimage import minimum_filter1d, maximum_filter1d
    from scipy.signal import savgol_filter
    w = max(int(len(counts) * window_frac), 11)
    bg = minimum_filter1d(counts.astype(float), size=w)
    bg = maximum_filter1d(bg, size=w)
    wsg = w | 1  # must be odd
    bg = savgol_filter(bg, window_length=min(wsg, len(counts) // 3 * 2 + 1),
                       polyorder=3)
    return np.clip(bg, 0, None)


def _draw_panel(ax, sp, color, ref_lines, is_last, subtitle=""):
    """Draw one XRD panel: normalised curve + reference sticks below y=0."""
    from matplotlib.lines import Line2D

    norm = sp["counts"] / sp["counts"].max()
    ax.plot(sp["two_theta"], norm, color=color, linewidth=1.0, zorder=3)
    # No fill/shading on the curve

    ref_bot  = -0.30   # bottom of reference stick zone (full height = 100%)
    text_y   = ref_bot + 0.02
    x_lo, x_hi = 10, 88

    # Reference sticks — only phases relevant to this sample
    if ref_lines:
        refs = _relevant_refs(sp["name"])
        for phase, info in refs.items():
            for pk, intensity in info["peaks"]:
                if x_lo <= pk <= x_hi:
                    # Scale stick height by relative intensity
                    stick_top = ref_bot * (1.0 - intensity / 100.0)
                    ax.vlines(pk, ref_bot, stick_top,
                              color=info["color"], linewidth=0.9,
                              alpha=0.80, zorder=2)

        # Text labels at bottom-left and bottom-right (no legend entries)
        phases = list(refs.items())
        positions = [x_lo + 0.5, x_hi - 0.5]    # left anchor, right anchor
        haligns   = ["left", "right"]
        for i, (phase, info) in enumerate(phases):
            px = positions[i % 2]
            ha = haligns[i % 2]
            ax.text(px, text_y, phase,
                    color=info["color"], fontsize=11,
                    ha=ha, va="bottom", zorder=5)

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(ref_bot - 0.02, 1.18)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)

    label = sp["name"] + (f"  ({subtitle})" if subtitle else "")
    ax.set_title(label, fontsize=13, loc="left", pad=4)

    if not is_last:
        ax.set_xticks([])
    else:
        ax.set_xlabel("2θ (°, Cu Kα)", fontsize=14)


def plot_stacked(spectra, output_path):
    """One subplot per sample: pattern + reference sticks inside each panel."""
    _pub_style()
    n = len(spectra)
    fig, axes = plt.subplots(nrows=n, ncols=1,
                              figsize=(10, 3.2 * n),
                              constrained_layout=True)
    if n == 1:
        axes = [axes]

    for idx, sp in enumerate(spectra):
        color = _PALETTE[idx % len(_PALETTE)]
        _draw_panel(axes[idx], sp, color,
                    ref_lines=True, is_last=(idx == n - 1))

    title = Path(output_path).stem.replace("_", " ").title()
    fig.suptitle(f"XRD Patterns – {title}", fontsize=15)
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_background_subtracted(spectra, output_path):
    """One subplot per sample: background-subtracted curve + reference sticks."""
    _pub_style()
    n = len(spectra)
    fig, axes = plt.subplots(nrows=n, ncols=1,
                              figsize=(10, 3.2 * n),
                              constrained_layout=True)
    if n == 1:
        axes = [axes]

    for idx, sp in enumerate(spectra):
        color  = _PALETTE[idx % len(_PALETTE)]
        bg     = _estimate_background(sp["counts"])
        net    = np.clip(sp["counts"] - bg, 0, None)
        sp_net = {**sp, "counts": net}
        _draw_panel(axes[idx], sp_net, color,
                    ref_lines=True, is_last=(idx == n - 1),
                    subtitle="bg-subtracted")

    title = Path(output_path).stem.replace("_", " ").title()
    fig.suptitle(f"XRD Patterns – {title} (Background Subtracted)", fontsize=15)
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_overlay(spectra, output_path):
    """
    Simple overlay without offset – useful for comparing peak positions.
    """
    _pub_style()
    fig, ax = plt.subplots(figsize=(10.4, 5.85))

    for idx, sp in enumerate(spectra):
        color = _PALETTE[idx % len(_PALETTE)]
        norm  = sp["counts"] / sp["counts"].max()
        ax.plot(sp["two_theta"], norm, color=color,
                linewidth=1.2, label=sp["name"])

    ax.set_xlim(5, 92)
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlabel("2θ (°, Cu Kα)")
    ax.set_ylabel("Normalised Intensity")
    title = Path(output_path).stem.replace("_", " ").title()
    ax.set_title(f"XRD Patterns – {title}")
    ax.legend(fontsize=10, loc="upper right", framealpha=0.8,
              bbox_to_anchor=(1.28, 1), borderaxespad=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Plot Bruker RAW4 XRD files")
    parser.add_argument("--data-dir",   required=True,
                        help="Folder containing .raw files (searched recursively)")
    parser.add_argument("--output-dir", default="output/xrd_raw4",
                        help="Where to write xrd_stacked.png and xrd_overlay.png")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_files = sorted(glob.glob(f"{args.data_dir}/**/*.raw", recursive=True))
    print(f"Found {len(raw_files)} .raw files\n")

    # De-duplicate by (n_pts, name stem)
    seen = set()
    spectra = []
    for f in raw_files:
        name = Path(f).stem
        data = open(f, "rb").read()
        key  = (len(data), name)
        if key in seen:
            continue
        seen.add(key)
        try:
            sp = parse_raw4(f)
            spectra.append(sp)
            print(f"  OK  {name:35s}  n={len(sp['counts']):5d}  "
                  f"2θ {sp['two_theta'][0]:.2f}–{sp['two_theta'][-1]:.2f}°  "
                  f"max={sp['counts'].max():.0f} cts")
        except Exception as e:
            print(f"  SKIP {name}: {e}")

    if not spectra:
        print("No spectra loaded.")
        return

    # Sort by sample name for consistent order
    spectra.sort(key=lambda s: s["name"])

    print(f"\nPlotting {len(spectra)} spectra …")
    plot_stacked(spectra, out_dir / "xrd_stacked.png")
    plot_background_subtracted(spectra, out_dir / "xrd_bg_subtracted.png")
    print("Done.")


if __name__ == "__main__":
    main()
