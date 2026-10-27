"""
read_bsw_uv.py  –  Parse Agilent/Varian Cary BSW binary files and plot UV-Vis spectra.

Binary format discovered by reverse-engineering:
  - File starts with length-prefixed string "Varian UV-VIS-NIR"
  - Followed by 9 × TContinuumStore(9/7/5) blocks, each 4405–4409 bytes
  - Each block holds 8-byte records: [absorbance (LE f32), wavelength_nm (LE f32)]
  - Spectral data begins at byte offset 0x413 (1043) within each block
  - A 'nacSd' marker at offset +0x47 in each block stores scan parameters
    (600 nm/min, 2 nm SBW, 800–500 nm range, 1 nm step → 301 points per spectrum)
  - The last portion of the file contains TGraphStore, TReportStore (RTF), and
    TDatabaseStore sections with additional metadata.

Usage:
    python read_bsw_uv.py <file.BSW> [--output-dir DIR]

Outputs:
    <output_dir>/<stem>_uv.png    – stacked multi-spectrum plot
    <output_dir>/<stem>_cal.png   – calibration curve (peak absorbance vs. label)
"""

import struct
import re
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ── Format constants ────────────────────────────────────────────────────────
DATA_OFFSET_IN_BLOCK = 0x413   # bytes from block-start to first 8-byte record
RECORD_SIZE          = 8       # bytes per (absorbance, wavelength) pair
MARKER               = b"TContinuumStor"


# ── Binary parser ───────────────────────────────────────────────────────────
def _find_blocks(data: bytes):
    """Return list of (offset, version_char) for every TContinuumStor block."""
    blocks = []
    for m in re.finditer(re.escape(MARKER), data):
        offset = m.start()
        version_byte = data[offset + len(MARKER) : offset + len(MARKER) + 1]
        blocks.append((offset, version_byte.decode("ascii", "replace")))
    return blocks


def _read_spectrum(data: bytes, block_start: int, n_max: int = 400):
    """
    Extract wavelength and absorbance arrays from one TContinuumStor block.
    Returns (wavelengths_nm, absorbance) numpy arrays (both ascending wl order).
    """
    base = block_start + DATA_OFFSET_IN_BLOCK
    wls, abs_ = [], []
    for i in range(n_max):
        off = base + i * RECORD_SIZE
        if off + 8 > len(data):
            break
        a  = struct.unpack_from("<f", data, off)[0]
        wl = struct.unpack_from("<f", data, off + 4)[0]
        if wl < 100 or wl > 1100:   # out-of-range → end of data
            break
        wls.append(wl)
        abs_.append(a)

    wls  = np.array(wls,  dtype=float)
    abs_ = np.array(abs_, dtype=float)

    # Sort by ascending wavelength
    idx  = np.argsort(wls)
    return wls[idx], abs_[idx]


def _get_labels(data: bytes, blocks):
    """
    Best-effort: pair each block with a sample label from the file's metadata.
    Falls back to "Scan N" for unlabelled blocks.
    """
    # Find all "NNN.NNmmol" or similar labels and the collection-time strings
    conc_re  = re.compile(rb"[\d.]+mmol")
    name_re  = re.compile(rb"([\d.\-]+(?:mmol|%))")
    time_re  = re.compile(rb"Collection Time: ([^\r\n]+)")

    # Gather label/time positions
    conc_hits  = [(m.start(), m.group().decode("ascii","replace").rstrip("%"))
                  for m in re.finditer(conc_re, data)]
    times      = [(m.start(), m.group(1).decode("ascii","replace"))
                  for m in re.finditer(time_re, data)]

    # Determine block boundaries
    block_starts = [b[0] for b in blocks]
    block_ends   = block_starts[1:] + [len(data)]

    labels = []
    for i, (bstart, _) in enumerate(blocks):
        bend = block_ends[i]
        # First concentration label inside this block
        label = None
        for pos, lbl in conc_hits:
            if bstart <= pos < bend:
                label = lbl
                break
        if label is None:
            label = f"Scan {i+1}"
        labels.append(label)
    return labels


def parse_bsw(filepath: str):
    """
    Parse a Cary BSW file.

    Returns
    -------
    spectra : list of dict
        Each dict has keys "label" (str), "wl" (ndarray), "abs" (ndarray).
    """
    data = Path(filepath).read_bytes()
    blocks = _find_blocks(data)
    if not blocks:
        raise ValueError("No TContinuumStor blocks found – not a Cary BSW file?")

    labels = _get_labels(data, blocks)
    spectra = []
    for (bstart, _), label in zip(blocks, labels):
        wl, ab = _read_spectrum(data, bstart)
        if len(wl) < 10:
            continue
        spectra.append({"label": label, "wl": wl, "abs": ab})
    return spectra


# ── Plotting ─────────────────────────────────────────────────────────────────
def _publication_style():
    plt.rcParams.update({
        "font.family":       "DejaVu Sans",
        "font.size":         10,
        "axes.linewidth":    1.2,
        "xtick.major.size":  4,
        "ytick.major.size":  4,
        "xtick.direction":   "in",
        "ytick.direction":   "in",
        "figure.dpi":        150,
    })


def plot_spectra(spectra, output_path, title="UV-Vis Absorbance Spectra"):
    """
    Multi-line overlay of all spectra with a colour ramp.
    """
    _publication_style()
    fig, ax = plt.subplots(figsize=(7, 4.5))

    colors = cm.viridis(np.linspace(0.1, 0.9, len(spectra)))
    for sp, color in zip(spectra, colors):
        ax.plot(sp["wl"], sp["abs"], color=color, linewidth=1.5,
                label=sp["label"])

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Absorbance (AU)")
    ax.set_title(title)
    ax.set_xlim(500, 800)
    ax.legend(framealpha=0.7, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.01, 1), borderaxespad=0)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_calibration(spectra, output_path, cal_wl=700):
    """
    Calibration curve: absorbance at cal_wl vs. concentration (if labels are numeric).
    """
    _publication_style()

    conc, abs_at_wl = [], []
    for sp in spectra:
        # Try to extract numeric concentration in mmol
        m = re.match(r"([\d.]+)mmol", sp["label"])
        if not m:
            continue
        c = float(m.group(1))
        # interpolate absorbance at cal_wl
        a = float(np.interp(cal_wl, sp["wl"], sp["abs"]))
        conc.append(c)
        abs_at_wl.append(a)

    if len(conc) < 2:
        return  # not enough points for a calibration curve

    conc = np.array(conc)
    abs_at_wl = np.array(abs_at_wl)

    # Linear fit
    coeffs = np.polyfit(conc, abs_at_wl, 1)
    poly   = np.poly1d(coeffs)
    x_fit  = np.linspace(conc.min() * 0.9, conc.max() * 1.05, 100)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(conc, abs_at_wl, s=60, zorder=5, color="#2C3E50",
               label="Data")
    ax.plot(x_fit, poly(x_fit), "--", color="#E74C3C", linewidth=1.5,
            label=f"y = {coeffs[0]:.3f}x + {coeffs[1]:.4f}")

    r2 = np.corrcoef(conc, abs_at_wl)[0, 1] ** 2
    ax.text(0.05, 0.92, f"$R^2$ = {r2:.4f}", transform=ax.transAxes,
            fontsize=9)

    ax.set_xlabel("Concentration (mmol)")
    ax.set_ylabel(f"Absorbance at {cal_wl} nm (AU)")
    ax.set_title("Calibration Curve")
    ax.legend(fontsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Plot Cary BSW UV-Vis files")
    parser.add_argument("bsw_file", help="Path to .BSW file")
    parser.add_argument("--output-dir", default=None,
                        help="Directory for output PNGs (default: same as input)")
    parser.add_argument("--cal-wl", type=float, default=700,
                        help="Wavelength (nm) used for calibration curve (default 700)")
    args = parser.parse_args()

    bsw_path = Path(args.bsw_file)
    out_dir  = Path(args.output_dir) if args.output_dir else bsw_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading: {bsw_path}")
    spectra = parse_bsw(str(bsw_path))
    print(f"  Found {len(spectra)} spectra:")
    for sp in spectra:
        print(f"    {sp['label']:15s}  wl {sp['wl'].min():.0f}–{sp['wl'].max():.0f} nm  "
              f"abs {sp['abs'].min():.3f}–{sp['abs'].max():.3f}")

    stem = bsw_path.stem
    plot_spectra(spectra, out_dir / f"{stem}_uv.png",
                 title=f"UV-Vis Spectra – {stem}")
    plot_calibration(spectra, out_dir / f"{stem}_cal.png",
                     cal_wl=args.cal_wl)


if __name__ == "__main__":
    main()
