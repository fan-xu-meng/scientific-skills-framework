---
name: tem-analysis
description: Analyze TEM and HRTEM images for nanoparticle size and candidate lattice spacings. Use when calibrated DM3/DM4 images, FFTs, particle masks, or possible crystal-plane assignments need review.
metadata:
  version: "1.0.0"
---

# TEM particle and lattice analysis

1. Read image calibration and pixel units from the source metadata. Verify scale bars against that calibration before measuring.
2. Define particle-selection and overlap rules. Segment or annotate particles and inspect the mask on representative images; retain rejected and ambiguous candidates.
3. Report diameter or radius with the stated estimator, sample count, image count, and distribution. Keep uncertainty from segmentation separate from between-image variation.
4. For HRTEM, calculate FFT or real-space fringe spacing from calibrated pixels. Match only against appropriate reference spacings and report tolerance.
5. Treat a spacing match as a candidate reflection. Confirm zone axis and image orientation before making exposed-facet claims.

Provide source image IDs, measurement overlays, size tables, and candidate reflection tables so a reviewer can reproduce the assignments.

