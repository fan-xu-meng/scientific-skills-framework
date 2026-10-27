---
name: xps-analysis
description: Review XPS survey and high-resolution regions, calibrate binding energy, fit constrained components, and quantify elemental signals. Use for VAMAS or exported XPS data when background and chemical-state assignments require audit.
metadata:
  version: "1.0.0"
---

# XPS fitting and quantification

1. Inspect survey and high-resolution spectra, acquisition settings, charging, and contamination. Preserve the original binding-energy and intensity arrays.
2. Choose a binding-energy reference appropriate to the sample and record the shift. Adventitious C 1s may be useful but is not universal, especially on carbon-rich or differentially charged materials.
3. Select and justify the background for each region. Fit the smallest chemically defensible set of components; constrain spin-orbit splitting, area ratios, and widths where physics requires them.
4. Inspect residuals, parameter bounds, and sensitivity to background or fit-window choices. Avoid assigning a chemical state from peak position alone.
5. For atomic fractions, use the appropriate sensitivity factors and report excluded regions, uncertainty, and any elements outside detection.

Export raw points, background, components, total fit, residuals, parameter table, and source provenance for each region.

