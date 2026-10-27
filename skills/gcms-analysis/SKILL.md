---
name: gcms-analysis
description: Parse GC-MS chromatograms and mass spectra, inspect TIC peaks, and compare isotope channels. Use with Agilent MassHunter exports or similar GC-MS datasets when peak identity and isotope evidence need review.
metadata:
  version: "1.0.0"
---

# GC-MS analysis

1. Inventory runs, blanks, standards, and failed acquisitions. Preserve each source filename and acquisition time.
2. Parse TIC and mass-spectral channels using the exported schema. Check scan counts, time units, ion masses, and detector saturation.
3. Select peak windows from chromatograms before reading isotope ratios. Use standards and blanks to assess background and overlapping fragments.
4. For isotope comparisons, report the exact ion channels, integration windows, natural-abundance or fragmentation corrections, and uncertainty. Keep raw and corrected ratios separate.
5. Plot TIC overlays and mass spectra at stated retention times. Mark excluded runs and unresolved peak assignments.

Do not infer a molecule's isotope origin from a single ion channel without checking fragmentation and coelution.

The [isotope-ratio audit](references/isotope-ratio-audit.md) and [synthetic ion table](examples/synthetic-ion-table.csv) describe an exported-data check. The [scan-ratio helper](scripts/scan_ion_ratios.py) computes raw channel ratios only; it does not replace raw-scan parsing or compound identification.

