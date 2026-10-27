---
name: ic-hcoo-co3-analysis
description: Calibrate ion-chromatography peaks and quantify formate and carbonate species in standards and samples. Use for Metrohm or Chromeleon IC exports with dilution, blank, and time-point pairing.
metadata:
  version: "1.0.0"
---

# IC formate and carbonate analysis

1. Identify analytes and detector channels from standards. Check retention-time stability and integration-window overlap.
2. Derive calibration concentrations from preparation records. Fit peak area against concentration using the stated range and model; keep raw areas.
3. Apply blank and dilution corrections. Pair time points only when sample identifiers and collection records support the match.
4. Report concentration, uncertainty, calibration residuals, and any below-range or unresolved peaks separately.
5. For carbonate, state whether the reported quantity is carbonate, bicarbonate, or total inorganic carbon and the chemistry used to convert between them.

The workflow should remain reproducible from chromatograms, standard recipes, and a saved parameter table.

Use the [IC calibration contract](references/calibration-and-speciation.md) and [synthetic standards](examples/synthetic-standards.csv) to check units and carbonate terminology before quantification. The [standard-table fit helper](scripts/fit_standard_table.py) is only an arithmetic check after peak integration.

