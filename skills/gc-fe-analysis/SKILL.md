---
name: gc-fe-analysis
description: Integrate gas-chromatography peaks, calibrate detector response, align injections with electrochemical intervals, and calculate gas-product Faradaic efficiency. Use for GC-supported electrochemical product analysis.
metadata:
  version: "1.0.0"
---

# GC calibration and gas FE

1. Inspect raw chromatograms and detector/channel metadata. Choose integration windows from standards, blanks, and sample retention times; save window diagnostics.
2. Build a detector-specific area-to-composition calibration from measured standards. Report regression form, range, residuals, and excluded standards. Do not substitute an instrument's stored concentration table without checking its standard preparation.
3. Convert calibrated composition to product molar flow using the measured outlet flow, pressure and temperature convention, and sampling interval.
4. Align injection time to electrochemical charge. Apply the electron count for each gas product, then compute `FE_i = n_e,i × F × n_i / |Q|`.
5. Flag coelution, saturation, retention-time drift, integration boundaries, missing injections, and FE closure failures.

Keep the raw peak areas, calibration table, aligned charge table, and final FE table as separate reproducible outputs.

The [OpenLAB GC parser](scripts/gc_openlab_adaptive.py) supports `.rslt` sequences. Its compound names and retention windows are method defaults that must be reviewed against standards. See the [original parser notes](references/openlab-qc.md) and [synthetic example](examples/synthetic-peak-table.csv).

