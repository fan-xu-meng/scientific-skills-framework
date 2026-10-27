---
name: xrd-analysis
description: Parse Bruker RAW4 or text XRD patterns, compare peaks with reference phases, and prepare annotated diffraction plots. Use for phase checks, temperature series, or powder-diffraction reports.
metadata:
  version: "1.0.0"
---

# XRD pattern analysis

1. Preserve acquisition metadata: radiation, wavelength, scan range, step size, and instrument geometry. Verify 2θ and intensity units.
2. Inspect the raw pattern before smoothing, background subtraction, or peak picking. Record preprocessing settings and retain the raw trace.
3. Compare candidate phases using reference peak positions and relative intensities. Account for peak shifts, preferred orientation, overlap, amorphous background, and finite crystallite size.
4. Annotate candidate reflections and distinguish an unambiguous match from a tentative phase assignment. Do not claim a phase from one overlapping peak.
5. Export a peak table, phase-reference IDs, matching tolerances, plots, and unresolved features.

For Bruker RAW4 input, the bundled [RAW4 parser](scripts/bruker_raw4_analyzer.py) accepts explicit data and output folders:

```bash
python scripts/bruker_raw4_analyzer.py --data-dir <raw-files> --output-dir <output>
```

