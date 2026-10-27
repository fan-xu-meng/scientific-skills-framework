---
name: nmr-hcook-analysis
description: Quantify formate or potassium formate and related liquid products from 1H NMR spectra. Use with Bruker raw or processed spectra, or exported peak tables, when calibration and spectral quality control matter.
metadata:
  version: "1.0.0"
---

# Quantitative 1H NMR

1. Inventory standards, blanks, internal references, and samples. Preserve Bruker acquisition and processing metadata, including pulse program, number of scans, receiver gain, and relaxation delay.
2. Check phase, baseline, chemical-shift reference, and peak overlap before integration. Use one consistent integration method across standards and samples.
3. Build a matrix-matched concentration calibration from known standards. Document any internal-standard concentration and response factor.
4. Integrate the formate signal and other requested products within justified windows. Apply dilution and sample-volume corrections, then report concentration with calibration uncertainty.
5. Save an overlay, integration-window plot, calibration curve, sample table, and exclusions. Flag saturated, poorly phased, or overlapped spectra.

For Faradaic efficiency, pair the quantified amount with the corresponding electrochemical charge interval and explicitly state the electron count.

