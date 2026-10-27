#!/usr/bin/env python3
"""Read cleaned TEM DM3/DM4 images and create quantitative QC summaries."""

from pathlib import Path
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ncempy.io import dm

DATA_ROOT = Path(os.environ["TEM_DATA_DIR"]).resolve()
OUT_ROOT = Path(os.environ["TEM_OUTPUT_DIR"]).resolve() / "overview"


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", text).strip("_")[:150]


def scalar_size(value: object) -> object:
    if isinstance(value, (list, tuple)):
        return "x".join(str(v) for v in value)
    try:
        return int(value)
    except Exception:
        return str(value)


def normalize_image(data: np.ndarray) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    if arr.ndim > 2:
        arr = arr.squeeze()
    if arr.ndim != 2:
        return np.empty((0, 0))
    low, high = np.nanpercentile(arr, [1, 99])
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        low, high = np.nanmin(arr), np.nanmax(arr)
    return np.clip((arr - low) / (high - low), 0, 1) if high > low else np.zeros_like(arr)


def fft_image(data: np.ndarray) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    arr = arr - np.nanmean(arr)
    ft = np.fft.fftshift(np.fft.fft2(arr))
    return np.log1p(np.abs(ft))


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    image_dir = OUT_ROOT / "images"
    fft_dir = OUT_ROOT / "fft"
    image_dir.mkdir(exist_ok=True)
    fft_dir.mkdir(exist_ok=True)
    rows = []
    files = sorted(list(DATA_ROOT.rglob("*.dm3")) + list(DATA_ROOT.rglob("*.dm4")))
    for path in files:
        rel = path.relative_to(DATA_ROOT)
        record = {"file": str(rel), "file_name": path.name, "status": "ok"}
        try:
            obj = dm.fileDM(str(path))
            record.update({"format": path.suffix.lower().lstrip("."), "file_size_bytes": path.stat().st_size,
                           "dataset_count": int(obj.numObjects), "x_size": scalar_size(obj.xSize),
                           "y_size": scalar_size(obj.ySize), "z_size": scalar_size(obj.zSize)})
            ds = obj.getDataset(0)
            data = np.asarray(ds["data"])
            if data.ndim > 2:
                data = data.squeeze()
            if data.ndim != 2:
                raise ValueError(f"dataset is not 2D after squeeze: {data.shape}")
            pix = ds.get("pixelSize", np.nan)
            unit = ds.get("pixelUnit", "")
            try:
                pix_values = np.asarray(pix, dtype=float).ravel()
                pix = float(np.nanmean(pix_values)) if len(pix_values) else np.nan
            except Exception:
                pix = np.nan
            record.update({"data_shape": str(tuple(data.shape)), "dtype": str(data.dtype),
                           "pixel_size": pix, "pixel_unit": str(unit),
                           "intensity_min": float(np.nanmin(data)), "intensity_max": float(np.nanmax(data)),
                           "intensity_mean": float(np.nanmean(data)), "intensity_std": float(np.nanstd(data)),
                           "intensity_p01": float(np.nanpercentile(data, 1)), "intensity_p99": float(np.nanpercentile(data, 99))})
            label = safe_name(str(rel.with_suffix("")))
            img = normalize_image(data)
            plt.imsave(image_dir / f"{label}.png", img, cmap="gray", vmin=0, vmax=1)
            ft = fft_image(data)
            plt.imsave(fft_dir / f"{label}_fft.png", ft, cmap="magma")
            record["analysis_note"] = "2D image; FFT saved for qualitative periodicity/lattice review"
        except Exception as exc:
            record.update({"status": "error", "error": str(exc)})
        rows.append(record)

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_ROOT / "tem_dm_summary.csv", index=False)
    ok = summary[summary["status"].eq("ok")].copy() if not summary.empty else summary
    if not ok.empty:
        # A compact contact sheet makes morphology differences easier to inspect.
        thumbs = []
        for _, row in ok.iterrows():
            path = image_dir / f"{safe_name(str(Path(row['file']).with_suffix('')))}.png"
            if path.exists():
                thumbs.append((row["file_name"], plt.imread(path)))
        cols = 4
        rows_n = max(1, int(np.ceil(len(thumbs) / cols)))
        fig, axes = plt.subplots(rows_n, cols, figsize=(14, 3.2 * rows_n))
        axes = np.atleast_1d(axes).ravel()
        for ax in axes:
            ax.axis("off")
        for ax, (name, image) in zip(axes, thumbs, strict=False):
            ax.imshow(image, cmap="gray")
            ax.set_title(name, fontsize=8)
            ax.axis("off")
        fig.suptitle("TEM image contact sheet")
        fig.tight_layout()
        fig.savefig(OUT_ROOT / "tem_dm_contact_sheet.png", dpi=170)
        plt.close(fig)
    (OUT_ROOT / "tem_dm_analysis_method.txt").write_text(
        "TEM DM3/DM4 analysis method\n\n"
        "All DM3 and DM4 files under the selected input folder were read with an Electron Microscopy DM reader. "
        "Dataset 0 was treated as the primary image. Image dimensions, dtype, intensity statistics, and pixel "
        "calibration metadata were recorded. Each image was exported with 1-99 percentile contrast for visual QC. "
        "A log-magnitude 2D FFT was also generated for qualitative periodicity or lattice review. FFT peaks were "
        "not converted to d-spacing unless a valid pixel calibration and diffraction geometry were available.\n",
        encoding="utf-8",
    )
    print(f"Files found: {len(files)}; readable: {int((summary['status'] == 'ok').sum())}; output: {OUT_ROOT}")


if __name__ == "__main__":
    main()
