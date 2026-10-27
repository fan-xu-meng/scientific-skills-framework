---
name: corr-fe-analysis
description: Calculate gas and liquid Faradaic efficiency for CO reduction experiments from charge, calibrated GC data, and calibrated HPLC or NMR products. Use when the feed is CO rather than CO2.
metadata:
  version: "1.0.0"
---

# CO reduction Faradaic efficiency

Use CO as the feed and do not count unreacted CO as a product. Keep the electron-count and carbon-accounting assumptions visible for each reported species.

1. Map every GC injection and liquid sample to an electrochemical interval, documenting flow, delay, collection volume, and dilution.
2. Convert detector response to amount using run-specific calibration and blank correction. Keep gas and liquid uncertainties separate.
3. Calculate `FE_i = n_e,i × F × n_i / |Q|` for each product. For cumulative liquids, use cumulative charge or a justified interval difference consistently.
4. Report gas FE, liquid FE, total measured FE, and a separate carbon balance. Flag missing products and apparent totals above 100% for review.
5. Save a traceable table with source filenames, timestamps, calibration IDs, charge windows, electron counts, and final values.

If CO2 is detected, distinguish feed impurity or crossover from products of CO conversion before assigning a pathway.

The [FE calculator](scripts/calculate_fe.py) checks arithmetic after calibration and charge-window alignment. Its [synthetic input](examples/synthetic-fe-input.csv) and [input contract](references/fe-input-contract.md) illustrate the format; unreacted CO must not appear as a product row.

