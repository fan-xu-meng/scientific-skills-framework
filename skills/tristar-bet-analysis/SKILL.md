---
name: tristar-bet-analysis
description: Extract and audit adsorption/desorption isotherms and BET surface-area results from Micromeritics TriStar II SMP and REP files. Use when a reproducible table and fit plot are needed from native instrument files.
metadata:
  version: "1.0.0"
---

# TriStar II BET analysis

1. Pair each sample's SMP and REP files. Keep the instrument report and sample identifier with the extracted values.
2. Parse relative pressure and adsorbed quantity. Check that arrays align and identify adsorption and desorption branches.
3. Reproduce the BET linear fit over the chosen pressure range. Report that range, slope, intercept, fit diagnostics, and the surface area from both fit and instrument report.
4. Preserve Langmuir, t-plot, and BJH values as reported, with units and the method labels; do not present them as interchangeable with BET area.
5. Save the isotherm, fit plot, numeric tables, and any parse/QC holds.

The bundled [TriStar processor](scripts/process_all_bet.py) accepts explicit input and output folders:

```bash
python scripts/process_all_bet.py --data-root <SMP-and-REP-folder> --out-dir <output>
```

Read the [fit-range checklist](references/bet-fit-qc.md) before accepting an automated result. The [synthetic isotherm rows](examples/synthetic-isotherm.csv) illustrate the output schema, not a valid BET fit.

