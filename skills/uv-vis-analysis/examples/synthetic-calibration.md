# Synthetic absorbance calibration

Invented values for documentation only.

| standard concentration (mM) | blank-corrected absorbance at 500 nm |
| ---: | ---: |
| 0 | 0.00 |
| 1 | 0.10 |
| 2 | 0.20 |
| 3 | 0.30 |

The illustrative relation is `A = 0.10 × concentration_mM`. A sample with `A=0.15` would be `1.5 mM` before dilution. A sample with `A=0.45` is outside the calibrated range; do not report `4.5 mM` without flagging extrapolation or rerunning a diluted sample. Real analysis also needs baseline, wavelength, residual, and saturation checks.
