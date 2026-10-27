#!/usr/bin/env python3
"""Pt/C particle and lattice candidates from calibrated DM3/DM4 images."""

from __future__ import annotations

from pathlib import Path
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.signal import find_peaks
from skimage import filters, measure, morphology, segmentation
from ncempy.io import dm

DATA_ROOT = Path(os.environ["TEM_DATA_DIR"]).resolve()
OUT_ROOT = Path(os.environ["TEM_OUTPUT_DIR"]).resolve() / "pt-candidates"
PT_FCC_D = {"111": 0.226, "200": 0.196, "220": 0.139}
PARTICLE_PIXEL_SIZE_MIN_NM = 0.02
PARTICLE_PIXEL_SIZE_MAX_NM = 0.08


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", text).strip("_")[:160]


def load_primary(path: Path) -> tuple[np.ndarray, float, str]:
    obj = dm.fileDM(str(path))
    ds = obj.getDataset(0)
    arr = np.asarray(ds["data"], dtype=float).squeeze()
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2D image, got shape {arr.shape}")
    pix = np.asarray(ds.get("pixelSize", [np.nan]), dtype=float).ravel()
    pixel_size = float(np.nanmean(pix)) if len(pix) else np.nan
    units = np.asarray(ds.get("pixelUnit", "")).ravel()
    unique_units = {str(value).strip().lower() for value in units}
    if len(unique_units) != 1:
        raise ValueError(f"Inconsistent pixel units: {sorted(unique_units)}")
    unit = unique_units.pop()
    factors = {"nm": 1.0, "nanometer": 1.0, "nanometre": 1.0,
               "å": 0.1, "a": 0.1, "angstrom": 0.1, "µm": 1000.0,
               "μm": 1000.0, "um": 1000.0, "m": 1e9}
    if unit not in factors or not np.isfinite(pixel_size) or pixel_size <= 0:
        raise ValueError(f"Unusable pixel calibration: {pixel_size} {unit!r}")
    return arr, pixel_size * factors[unit], unit


def norm(arr: np.ndarray) -> np.ndarray:
    lo, hi = np.nanpercentile(arr, [1, 99])
    return np.clip((arr - lo) / (hi - lo), 0, 1) if hi > lo else np.zeros_like(arr)


def remove_objects_smaller_than(mask: np.ndarray, minimum_size: int) -> np.ndarray:
    try:
        return morphology.remove_small_objects(mask, max_size=minimum_size - 1)
    except TypeError:
        return morphology.remove_small_objects(mask, min_size=minimum_size)


def segment_particles(arr: np.ndarray, pixel_size: float) -> tuple[np.ndarray, list[measure._regionprops.RegionProperties]]:
    image = norm(arr)
    # Strongly dark pixels seed particles; a lower threshold grows each seed to
    # the full locally dark boundary. This avoids measuring only the dark core.
    denoise_sigma = max(0.8, min(2.0, 0.06 / max(pixel_size, 1e-4)))
    background_sigma = max(12.0, min(48.0, 1.2 / max(pixel_size, 1e-4)))
    denoised = filters.gaussian(image, sigma=denoise_sigma)
    background = filters.gaussian(denoised, sigma=background_sigma)
    dark_score = background - denoised
    finite_score = dark_score[np.isfinite(dark_score)]
    if finite_score.size < 100:
        return np.zeros_like(image, dtype=int), []
    high_threshold = max(0.040, float(np.nanpercentile(finite_score, 97.5)))
    low_threshold = max(0.012, 0.45 * high_threshold)
    strong_seeds = dark_score >= high_threshold
    min_area_px = max(30, int(np.pi * (0.40 / max(pixel_size, 1e-4)) ** 2))
    strong_seeds = morphology.opening(strong_seeds, morphology.disk(2))
    strong_seeds = morphology.closing(strong_seeds, morphology.disk(3))
    strong_seeds = remove_objects_smaller_than(strong_seeds, min_area_px)
    if not np.any(strong_seeds):
        return np.zeros_like(image, dtype=int), []
    dark = filters.apply_hysteresis_threshold(dark_score, low_threshold, high_threshold)
    max_growth_px = max(8, int(round(3.0 / max(pixel_size, 1e-4))))
    dark &= ndi.distance_transform_edt(~strong_seeds) <= max_growth_px
    close_radius = max(2, int(round(0.12 / max(pixel_size, 1e-4))))
    dark = morphology.closing(dark, morphology.disk(close_radius))
    dark = morphology.opening(dark, morphology.disk(max(1, close_radius // 2)))
    dark = remove_objects_smaller_than(dark, min_area_px)
    dark = morphology.remove_small_holes(dark, max_size=max(100, 2 * min_area_px))

    distance = ndi.distance_transform_edt(dark)
    markers = measure.label(strong_seeds & dark)
    labels = segmentation.watershed(-distance, markers, mask=dark)

    max_area_px = np.pi * (6.0 / max(pixel_size, 1e-4)) ** 2
    kept_labels = np.zeros_like(labels, dtype=int)
    kept_count = 0
    for p in measure.regionprops(labels, intensity_image=image):
        if p.area < min_area_px or p.area > max_area_px:
            continue
        if p.bbox[0] == 0 or p.bbox[1] == 0 or p.bbox[2] == image.shape[0] or p.bbox[3] == image.shape[1]:
            continue
        kept_count += 1
        kept_labels[labels == p.label] = kept_count
    return kept_labels, list(measure.regionprops(kept_labels, intensity_image=image))


def fft_candidates(arr: np.ndarray, pixel_size: float) -> list[dict[str, object]]:
    image = norm(arr)
    h, w = image.shape
    crop = image[max(0, h // 2 - min(h, w) // 4): h // 2 + min(h, w) // 4,
                 max(0, w // 2 - min(h, w) // 4): w // 2 + min(h, w) // 4]
    crop = crop - filters.gaussian(crop, sigma=8)
    ft = np.abs(np.fft.fftshift(np.fft.fft2(crop)))
    cy, cx = np.array(ft.shape) // 2
    yy, xx = np.indices(ft.shape)
    radius = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    radial = []
    for r in range(3, int(min(ft.shape) * 0.45)):
        band = ft[(radius >= r - 0.5) & (radius < r + 0.5)]
        radial.append(np.nanmedian(band) if len(band) else np.nan)
    radial = np.asarray(radial)
    if len(radial) < 5 or not np.isfinite(pixel_size):
        return []
    smooth = filters.gaussian(radial, sigma=1)
    peaks, props = find_peaks(smooth, prominence=max(np.nanstd(smooth) * 0.5, 1e-9), distance=4)
    order = peaks[np.argsort(smooth[peaks])[::-1]][:8]
    rows = []
    n = min(crop.shape)
    for peak in order:
        freq = float((peak + 3) / n)
        d_nm = float(pixel_size / freq) if freq > 0 else np.nan
        nearest = min(PT_FCC_D, key=lambda k: abs(PT_FCC_D[k] - d_nm))
        delta = abs(PT_FCC_D[nearest] - d_nm)
        rows.append({"fft_radius_px": int(peak + 3), "d_spacing_nm": d_nm,
                     "fft_relative_strength": float(smooth[peak] / np.nanmax(smooth)),
                     "nearest_pt_fcc_reflection": nearest if delta <= 0.015 else "unknown",
                     "reference_d_spacing_nm": PT_FCC_D[nearest], "assignment_delta_nm": delta,
                     "assignment_confidence": "candidate_only"})
    return rows


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    particle_dir = OUT_ROOT / "particle_qc_overlays"
    fft_dir = OUT_ROOT / "fft_qc"
    particle_dir.mkdir(exist_ok=True); fft_dir.mkdir(exist_ok=True)
    particle_rows, fft_rows, image_rows, manual_rows = [], [], [], []
    files = sorted(DATA_ROOT.rglob("*.dm3")) + sorted(DATA_ROOT.rglob("*.dm4"))
    for path in files:
        try:
            arr, pixel_size, pixel_unit = load_primary(path)
            rel = safe_name(str(path.relative_to(DATA_ROOT))); label = safe_name(str(path.relative_to(DATA_ROOT).with_suffix("")))
            if pixel_size < PARTICLE_PIXEL_SIZE_MIN_NM:
                analysis_group = "high_resolution_lattice_candidate"
            elif pixel_size <= PARTICLE_PIXEL_SIZE_MAX_NM:
                analysis_group = "low_medium_particle_size_analysis"
            else:
                analysis_group = "coarse_overview_not_used_for_size"
            image_rows.append({"file": rel, "file_name": path.name, "pixel_size_nm": pixel_size, "source_pixel_unit": pixel_unit,
                               "image_shape": str(arr.shape), "analysis_group": analysis_group})
            labels, props = (segment_particles(arr, pixel_size) if analysis_group == "low_medium_particle_size_analysis"
                             else (np.zeros_like(arr, dtype=int), []))
            for idx, p in enumerate(props, 1):
                area_px = float(p.area); eq_d_px = 2 * np.sqrt(area_px / np.pi)
                feret_px = float(getattr(p, "feret_diameter_max", np.nan))
                minor_axis_px = float(getattr(p, "axis_minor_length", np.nan))
                diameter_nm = eq_d_px * pixel_size
                radius_nm = 0.5 * diameter_nm
                roundness = float(4 * np.pi * area_px / (p.perimeter ** 2)) if p.perimeter else np.nan
                particle_rows.append({"file": rel, "file_name": path.name, "particle_id": idx,
                    "area_px2": area_px, "area_nm2": area_px * pixel_size ** 2,
                    "equivalent_radius_nm": radius_nm, "equivalent_diameter_nm": diameter_nm,
                    "feret_max_nm": feret_px * pixel_size,
                    "feret_min_nm": minor_axis_px * pixel_size,
                    "aspect_ratio": feret_px / minor_axis_px if minor_axis_px else np.nan,
                    "roundness": roundness,
                    "centroid_x_px": float(p.centroid[1]), "centroid_y_px": float(p.centroid[0]),
                    "segmentation_method": "dark_core_hysteresis_growth_and_watershed",
                    "automated_qc_status": "recommended_review" if 0.8 <= diameter_nm <= 10.0 and roundness >= 0.25 else "likely_false_or_merged",
                    "manual_qc_pass": "pending"})
            if analysis_group == "low_medium_particle_size_analysis":
                fig, ax = plt.subplots(figsize=(8, 8)); ax.imshow(norm(arr), cmap="gray")
                ax.contour(labels > 0, levels=[0.5], colors="#ff4d4d", linewidths=0.45)
                ax.set_title(f"Particle candidates | {path.name} | n={len(props)}")
                ax.axis("off"); fig.tight_layout(); fig.savefig(particle_dir / f"{label}_particles.png", dpi=170); plt.close(fig)
            candidates = fft_candidates(arr, pixel_size)
            for rank, row in enumerate(candidates, 1):
                fft_rows.append({"file": rel, "file_name": path.name, "candidate_rank": rank, **row})
            if pixel_size < 0.02:
                ft = np.log1p(np.abs(np.fft.fftshift(np.fft.fft2(norm(arr) - filters.gaussian(norm(arr), sigma=8)))))
                plt.imsave(fft_dir / f"{label}_fft.png", ft, cmap="magma")
                manual_rows.append({"file": rel, "file_name": path.name, "pixel_size_nm": pixel_size,
                    "zone_axis": "manual_review_required", "facet_assignment": "manual_review_required",
                    "evidence": "lattice-resolved image and FFT candidate table", "confidence": "pending"})
        except Exception as exc:
            image_rows.append({"file": safe_name(str(path.relative_to(DATA_ROOT))), "file_name": path.name, "status": "error", "error": str(exc)})
    if not any(row.get("analysis_group") for row in image_rows):
        raise ValueError("No calibrated two-dimensional DM images could be analyzed; inspect input files and pixel units")
    pd.DataFrame(image_rows).to_csv(OUT_ROOT / "tem_analysis_image_groups.csv", index=False)
    pd.DataFrame(particle_rows).to_csv(OUT_ROOT / "particle_measurements_pending_qc.csv", index=False)
    particle_df = pd.DataFrame(particle_rows)
    if not particle_df.empty:
        accepted = particle_df[particle_df["automated_qc_status"].eq("recommended_review")]
        accepted.to_csv(OUT_ROOT / "particle_measurements_recommended_review.csv", index=False)
        overall_summary = pd.DataFrame(
            [{
                "candidate_count": len(particle_df),
                "recommended_review_count": len(accepted),
                "mean_equivalent_radius_nm": accepted["equivalent_radius_nm"].mean(),
                "median_equivalent_radius_nm": accepted["equivalent_radius_nm"].median(),
                "mean_equivalent_diameter_nm": accepted["equivalent_diameter_nm"].mean(),
                "median_equivalent_diameter_nm": accepted["equivalent_diameter_nm"].median(),
                "std_equivalent_diameter_nm": accepted["equivalent_diameter_nm"].std(),
                "minimum_equivalent_diameter_nm": accepted["equivalent_diameter_nm"].min(),
                "maximum_equivalent_diameter_nm": accepted["equivalent_diameter_nm"].max(),
            }]
        )
        overall_summary.to_csv(OUT_ROOT / "particle_size_summary.csv", index=False)
        accepted.groupby("file_name").agg(
            particle_count=("particle_id", "count"),
            mean_equivalent_radius_nm=("equivalent_radius_nm", "mean"),
            mean_equivalent_diameter_nm=("equivalent_diameter_nm", "mean"),
            median_equivalent_diameter_nm=("equivalent_diameter_nm", "median"),
            std_equivalent_diameter_nm=("equivalent_diameter_nm", "std"),
        ).reset_index().to_csv(OUT_ROOT / "particle_size_summary_by_image.csv", index=False)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(accepted["equivalent_diameter_nm"], bins=30, color="#245b8a", alpha=0.85)
        ax.set_xlabel("Equivalent circular diameter (nm)"); ax.set_ylabel("Candidate count")
        ax.set_title(f"Dark-particle candidates recommended for manual QC (n={len(accepted)})")
        ax.grid(axis="y", alpha=0.2); fig.tight_layout(); fig.savefig(OUT_ROOT / "particle_size_distribution_pending_qc.png", dpi=180); plt.close(fig)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(accepted["equivalent_radius_nm"], bins=30, color="#2a7f62", alpha=0.85)
        ax.set_xlabel("Equivalent circular radius (nm)"); ax.set_ylabel("Candidate count")
        ax.set_title(f"Dark-particle radius candidates recommended for manual QC (n={len(accepted)})")
        ax.grid(axis="y", alpha=0.2); fig.tight_layout(); fig.savefig(OUT_ROOT / "particle_radius_distribution_pending_qc.png", dpi=180); plt.close(fig)
    fft_columns = ["file", "file_name", "candidate_rank", "fft_radius_px", "d_spacing_nm",
                   "fft_relative_strength", "nearest_pt_fcc_reflection", "reference_d_spacing_nm",
                   "assignment_delta_nm", "assignment_confidence"]
    pd.DataFrame(fft_rows, columns=fft_columns).to_csv(
        OUT_ROOT / "lattice_spacing_candidates_pending_indexing.csv", index=False)
    pd.DataFrame(manual_rows).to_csv(OUT_ROOT / "facet_statistics_manual_review_template.csv", index=False)
    (OUT_ROOT / "tem_particle_facet_method.txt").write_text(
        "Pt/C TEM particle and reflection-candidate workflow\n\n"
        "Particle sizes were measured only from images calibrated between 0.02 and 0.08 nm/pixel; higher-resolution "
        "images were reserved for lattice review and coarser overview images were excluded from size statistics. Strong "
        "dark cores seeded a lower-threshold hysteresis growth mask so the measurement follows the complete locally dark "
        "particle boundary rather than only its darkest core. Distance-transform watershed separated touching particles. "
        "Measurements include area, equivalent circular radius and diameter, Feret maximum/minimum, aspect ratio, "
        "roundness, and centroid. All particle labels are "
        "marked manual_qc_pass=pending because automated masks must be reviewed against carbon-support contrast. The "
        "recommended-review subset uses 0.8-10.0 nm equivalent diameter and roundness >= 0.25 as a triage rule only.\n\n"
        "High-resolution images were assigned to FFT/lattice review when pixel size was below 0.02 nm/pixel. Radial FFT "
        "peaks were converted to candidate d-spacings using the DM calibration. Pt fcc reference checks were 0.226 nm "
        "for {111}, 0.196 nm for {200}, and 0.139 nm for {220}. These are candidate assignments only; they do not prove "
        "an exposed facet or zone axis. Facet statistics require manual boundary/facet tracing and are provided as a "
        "review template.\n",
        encoding="utf-8",
    )
    print(f"Files: {len(files)}; particle candidates: {len(particle_rows)}; FFT candidates: {len(fft_rows)}")


if __name__ == "__main__":
    main()
