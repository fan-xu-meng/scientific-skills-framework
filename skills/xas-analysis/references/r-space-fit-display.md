# R-space Fit Display and Leading-quarter Diagnostics

Use this reference whenever an EXAFS fit is performed in R space, especially
when the leading quarter of a raw k-space plot looks poor.

## Contents

1. Interpretation rule
2. Required transform procedure
3. Leading-quarter metric
4. Figure and table outputs
5. FEFF reference consistency
6. Acceptance checklist

## 1. Interpretation rule

An R-space fit minimizes residuals only over its declared R interval. AUTOBK
components below `Rbkg`, Fourier leakage, and signal outside `Rmin:Rmax` are not
structural FEFF targets. Therefore:

- require `Rmin >= Rbkg` for new Athena/AUTOBK structural fits;
- judge the structural model from `|chi(R)|`, Re/Im `chi(R)`, and the
  back-transform of the fitted R window;
- retain raw `k^w chi(k)` and full R-space as a separate processing diagnostic;
- never present an R-filtered curve as raw data or omit the full diagnostic;
- do not add FEFF paths merely to reproduce `R < Rmin` background intensity.

If both raw and R-filtered leading quarters are poor, investigate processing,
E0, path phase, and the chemical model. If only the raw leading quarter is
poor, first inspect AUTOBK/Rbkg and low-R background rather than expanding the
structural model.

## 2. Required transform procedure

Recompute a common display transform for data and model. This avoids a Larch
cache pitfall: after a multi-kweight fit, `dset.data.chir_*` may correspond to a
different k-weight from the one assumed by the plot.

```python
import larch
import numpy as np
from larch.xafs import xftf, xftr

def fitted_window_views(data_group, dset, session, display_kw=2):
    views = []
    for chi in (data_group.chi, dset.model.chi):
        view = larch.Group(
            k=np.asarray(data_group.k, dtype=float).copy(),
            chi=np.asarray(chi, dtype=float).copy(),
        )
        xftf(
            view,
            kmin=dset.transform.kmin,
            kmax=dset.transform.kmax,
            dk=dset.transform.dk,
            dk2=dset.transform.dk2,
            window=dset.transform.window,
            kweight=display_kw,
            _larch=session,
        )
        xftr(
            view,
            rmin=dset.transform.rmin,
            rmax=dset.transform.rmax,
            dr=dset.transform.dr,
            dr2=dset.transform.dr2,
            window=dset.transform.rwindow,
            _larch=session,
        )
        views.append(view)
    return tuple(views)  # data_view, model_view
```

Important: `xftf(..., kweight=2)` transforms `k^2*chi(k)`. Its `xftr` result
`chiq_re` is already the R-filtered `k^2*chi(k)`. Do not multiply `chiq_re` by
`q**2` again. Divide by `q**2` only when an unweighted filtered `chi(k)` column
is explicitly needed.

Use the exact reverse-window settings stored in `dset.transform`. Do not choose
a visually convenient `dr` after fitting, because that changes the displayed
information relative to the fit definition.

## 3. Leading-quarter metric

Define the leading quarter from the fitted k interval, not from the physical
file length:

```python
quarter_max = kmin + 0.25 * (kmax-kmin)

def relative_rmse(observed, predicted):
    observed = np.asarray(observed)
    predicted = np.asarray(predicted)
    scale = max(float(np.sqrt(np.mean(observed**2))), 1e-12)
    return float(np.sqrt(np.mean((observed-predicted)**2)) / scale)
```

Calculate and report:

- raw leading-quarter relative RMSE from raw `k^2*chi(k)`;
- fitted-R-window leading-quarter relative RMSE from `chiq_re`;
- fitted-R-window full-k relative RMSE;
- exact fitted k range and R filter range.

Use these metrics comparatively, not as universal pass/fail thresholds. Inspect
phase and residual shape. A low scalar RMSE cannot validate an unphysical path
model.

## 4. Figure and table outputs

Write two per-sample figures.

Primary fitted-window figure:

1. R-filtered `k^2*chi(k)` data, fit, and residual over fitted k;
2. `|chi(R)|` data and fit only over fitted R;
3. Re/Im `chi(R)` data and fit only over fitted R.

Full diagnostic figure:

1. raw `k^2*chi(k)` data, model, and residual;
2. full R-space magnitude;
3. full R-space Re/Im;
4. gray shading and labels for `R < Rmin` and `R > Rmax`.

Keep the word `R-filtered` in the primary k-axis label and `raw` in the full
diagnostic label. Use marker-only data and line-only fits.

Write:

```text
chi_k.csv                         # raw chi and raw model
chi_k_fit_window_filtered.csv     # R-filtered data/model/residual
ft_r.csv                          # freshly recomputed common-kw transforms
front_quarter_diagnostics.csv     # raw and filtered metrics
exafs_fit_<sample>.png            # primary fitted-window figure
exafs_fit_diagnostic_full_<sample>.png
```

## 5. FEFF reference consistency

When one dataset combines paths from independent FEFF runs, check their path
headers and potential treatment before adding E0 variables.

For Pt L3 Pt-O + FCC Pt models:

1. generate the Pt-O path with an SCF card such as `SCF 4.2 0 30 0.2`;
2. keep the metallic FCC calculation and the Pt-O calculation self-contained;
3. map paths by `nleg`, degeneracy, and `reff`, not filename alone;
4. start with one shared E0 bounded near the Pt L3 range (typically 5-16 eV);
5. compare non-SCF/shared-E0, SCF/shared-E0, and any justified separate-E0
   alternative with BIC, phase agreement, bounds, correlations, and window
   robustness;
6. accept separate E0 values only for documented distinct phases/reference
   offsets, not simply because they lower R-factor.

In the validated 2026-08-08 Pt/C case, SCF Pt-O plus shared E0 reduced R-factor
from about 0.0265 to 0.0136 and reduced the fitted-window leading-quarter
relative RMSE to about 0.083 without adding a fit variable.

## 6. Acceptance checklist

- Confirm source data remain unchanged and copy/hash provenance when requested.
- Preserve Athena normalization and all available AUTOBK spline/clamp settings.
- Require `Rmin >= Rbkg` for new fits.
- Recompute common-kw data/model transforms after robustness trials.
- Use exact fit reverse-window settings for the back-transform.
- Verify that filtered arrays were not double k-weighted.
- Inspect primary and full diagnostic figures visually.
- Record raw and filtered leading-quarter metrics.
- Run k/R-window and Rbkg sensitivity tests.
- Keep `Nvariables < Nidp`, inspect uncertainties/bounds/correlations, and
  retain QC class C wording when strong correlations remain.
