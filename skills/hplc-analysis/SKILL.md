---
name: hplc-analysis
description: Integrate Agilent OpenLAB HPLC RID chromatograms and quantify liquid products from standards. Use for organic-acid and alcohol calibration, peak review, and concentration reporting.
metadata:
  version: "1.0.0"
---

# HPLC RID analysis

1. Preserve raw chromatograms and the standard-preparation record. Confirm detector channel, injection volume, retention time, and baseline units.
2. Select compound windows from standards and inspect coelution. Plot integration boundaries for every calibration run and representative samples.
3. Fit each compound's area-to-concentration calibration within its measured range. Report fit type, residuals, blank, limit of quantification, and outliers.
4. Apply dilution and sample-volume corrections. Do not extrapolate beyond the calibration range without flagging the result.
5. Export peak-area, calibration, and concentration tables with source filenames and method parameters.

If OpenLAB provides a stored concentration, verify that it uses the same standards and units as the actual preparation sheet before using it.

For the result-container structure, read the [original OpenLAB parsing note](references/agilent-openlab-rslt.md). The [OpenLAB RID parser](scripts/hplc_openlab_adaptive.py) has method-specific compound names and windows; verify them with standards before use. The [synthetic calibration example](examples/synthetic-calibration.md) shows the reporting fields without publishing experiment data.

