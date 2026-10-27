---
name: xas-analysis
description: Process XANES and EXAFS datasets with traceable energy alignment, normalization, background removal, and constrained structural fitting. Use when synchrotron XAS or Athena project files need quality controlled comparison.
metadata:
  version: "1.0.0"
---

# XAS and EXAFS

1. Inventory scans, standards, calibration foils, edge energy, detector mode, and exclusions. Keep raw and processed spectra separately.
2. Align energy using a documented reference, normalize XANES with recorded pre- and post-edge ranges, and inspect scan agreement before averaging.
3. For EXAFS, record background-removal parameters, k range, k weighting, transform window, and R range. Check sensitivity to reasonable parameter changes.
4. Construct structural paths from an identified model and constrain shared physical parameters consistently. Keep the number of independent fit parameters below the information content of the fitted range.
5. Report coordination, distance, disorder, energy shift, amplitude assumptions, uncertainties, correlations, and residuals. Compare fit and data in k space and both magnitude and real/imaginary R-space components.
6. State which conclusions are model dependent and retain alternative fits when they remain plausible.

Do not treat a low residual alone as evidence for a unique coordination structure.

