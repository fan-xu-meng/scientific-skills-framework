---
name: tem-analysis
description: Analyze calibrated TEM and HRTEM images, including DM3/DM4 previews, nanoparticle-size candidates, FFT lattice-spacing candidates, and cautious reflection or facet review. Use for image QC, particle masks, size distributions, or crystallographic checks.
metadata:
  version: "1.1.0"
---

# TEM analysis

This single skill covers image inventory and, when appropriate, dark-particle and Pt/C-specific FFT candidate analysis. Keep raw images unchanged. Record the image ID, calibration, pixel units, selection rules, and every rejected or ambiguous field.

## Choose a mode

- For any DM3/DM4 dataset, make metadata, image, FFT, and contact-sheet previews with the [runner](scripts/run_tem_analysis.py):

  ```bash
  python scripts/run_tem_analysis.py --data-dir <dm-folder> --out-dir <output-folder>
  ```

- For Pt/C images where dark particles are the target, add `--pt-candidates`. This also creates particle masks, candidate size tables, FFT d-spacing tables, Pt fcc reflection matches, and manual-review templates. These numerical thresholds are a starting profile, not universal material properties. Read [the particle and FFT method](references/particle-and-fft-qc.md) before using it.

The runner uses `ncempy`, `numpy`, `pandas`, `scipy`, `scikit-image`, and `matplotlib`, installing them through `uv` if needed. It takes explicit input/output paths and has no built-in study directory.

## Interpret and report

1. Verify scale and pixel unit from DM metadata against any visible scale bar. Reject uncalibrated images from nanometre measurements.
2. Define image and particle inclusion criteria before summarizing size. Inspect contour overlays; keep all automated values marked `manual_qc_pass=pending` until a human checks the original images.
3. State whether size is radius, equivalent circular diameter, or Feret diameter; report particle count, image count, uncertainty, and between-image variation. Do not pool unlike magnifications without justification.
4. Compare FFT or fringe spacings with an identified reference and stated tolerance. A d-spacing match is a *reflection candidate*, not proof of an exposed surface facet, zone axis, or 3D surface-area fraction.
5. Assign exposed facets only after crystallographic orientation and particle-boundary tracing are reviewed. Report projected boundary-length fractions separately from 3D surface fractions.

The [synthetic QC example](examples/synthetic-qc.md) illustrates the distinction between an automated size candidate and a confirmed facet claim.

