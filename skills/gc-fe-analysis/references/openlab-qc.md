# OpenLAB GC parser checks

`scripts/gc_openlab_adaptive.py` reads Agilent `.rslt` sequences with `.rx` integration results and `.dx` raw channel traces. Before processing a new method, inspect its TCD channel names, standard concentrations, retention-time drift, and compound mapping. The script's built-in windows are starting values from one method, not portable calibration constants.

Retain chromatogram-window plots, standard areas, regression residuals, blank response, and excluded injections. Pair calibrated gas composition with measured flow and the matching charge interval only after these parser checks. A concentration reported by the instrument method is not independent confirmation of the manually prepared standard.
