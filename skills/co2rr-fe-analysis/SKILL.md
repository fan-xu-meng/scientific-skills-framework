---
name: co2rr-fe-analysis
description: Calculate and audit product Faradaic efficiency and carbon balance for CO2 electroreduction. Use when combining electrochemical charge with calibrated gas and liquid product measurements.
metadata:
  version: "1.0.0"
---

# CO2 reduction Faradaic efficiency

Work from measured charge, calibrated product amount, and an explicit reaction stoichiometry. Keep the raw measurement, calibration result, and final FE in separate tables.

1. Align each gas injection or liquid sample to the charge interval it represents. State any residence-time or sampling-delay correction.
2. Calculate product moles using the appropriate gas flow and calibration, or liquid concentration, volume, dilution, and sampling history.
3. For product `i`, compute `FE_i = n_e,i × F × n_i / |Q|`. State the electron count and whether `Q` is interval or cumulative charge.
4. Sum measured product efficiencies. Compute carbon balance separately, using each product's carbon count and a stated carbon-input basis.
5. Report calibration range, detection limits, blank corrections, uncertainty, unassigned products, and any FE or carbon-balance closure gap.

Do not use stored calibration coefficients from another instrument or date without a documented cross-check. Distinguish CO2 feed, CO2 carryover, and carbon in products when defining utilization.

