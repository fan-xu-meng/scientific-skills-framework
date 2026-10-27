# Synthetic RID calibration

Invented values for documentation only; they are not instrument data.

| standard | concentration (mM) | blank-corrected area (nRIU·min) |
| --- | ---: | ---: |
| blank | 0 | 0 |
| S1 | 1 | 10 |
| S2 | 2 | 20 |
| S3 | 4 | 40 |

If a sample has area `25 nRIU·min`, the illustrative slope `10 nRIU·min/mM` gives `2.5 mM` before dilution correction. Save the actual fitted slope, intercept, residuals, calibration ID, chromatogram window, blank, and sample dilution. A sample above `40 nRIU·min` is outside this example range and requires dilution or a justified extrapolation flag.
