---
name: uv-vis-analysis
description: Parse Cary UV-Vis BSW or CSV spectra, build Beer-Lambert calibration, and plot spectra and concentrations. Use when standards and samples must be compared with wavelength and baseline checks.
metadata:
  version: "1.0.0"
---

# UV-Vis calibration

1. Identify wavelength, absorbance, blank, and sample columns; check export units and wavelength ordering.
2. Subtract or account for the blank according to the acquisition protocol. Inspect spectra for outliers, scattering, baseline drift, and saturation.
3. Choose an analytical wavelength or fit region before calculating sample concentration. Use standards with documented preparation and report the calibration range, model, residuals, and dilution.
4. Plot overlaid spectra and the calibration with each standard visible. Flag samples outside the calibrated range.
5. Export a tidy concentration table with source file, selected wavelength, absorbance, calibration ID, and QC status.

For Cary BSW, use the bundled [binary reader](scripts/read_bsw_uv.py). For CSV exports, use the [CSV analyzer](scripts/uv_cary_csv_adaptive.py) after checking that its expected column layout matches the file.

For the tested BSW layout and known limits, read the [BSW parsing reference](references/bsw-workflow.md). An invented [calibration example](examples/synthetic-calibration.md) illustrates range checks.

