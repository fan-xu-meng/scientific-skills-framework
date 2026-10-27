# BSW UV-Vis Workflow Reference

## 1. Instrument & File Format

**Instrument**: Agilent (Varian) Cary 5000 UV-VIS-NIR spectrophotometer
**Software**: Cary WinUV v6.1
**File format**: `.BSW` — Batch Scan Workspace, proprietary binary

The BSW format is not publicly documented. The structure below was determined by
binary reverse-engineering (hex inspection + float scanning).

---

## 2. BSW File Layout

```
[File header]
  0x00  0x11 "Varian UV-VIS-NIR"  (Pascal-style length-prefixed string)
  0x12  padding zeros

[N × TContinuumStor blocks]   ← spectral data lives here
  Each block: ~4405–4409 bytes
  Marker: b"TContinuumStor" + version_byte (b"9", b"7", or b"5")

[TGraphStore section]         ← graph layout metadata (not spectral data)
[TReportStore section]        ← RTF-formatted report text
[TDatabaseStore section]      ← method parameters, scan settings
```

### TContinuumStor Block Internal Structure

```
+0x00  "TContinuumStor" + version_byte   (15 bytes)
+0x0F  4 bytes block header
+0x47  "nacSd" marker                    (5 bytes)
         followed by scan parameters as LE f32:
           offset +20: 0.0997 (SBW-related)
           offset +24: 2.0    (SBW, nm)
           offset +28: -1.0
           offset +32: 600.0  (scan rate, nm/min)
           offset +36: 220.5
           offset +40: 2.0    (SBW, nm — confirmed)
           offset +44: 500.0  (X stop, nm)
           offset +48: 800.0  (X start, nm)

+0x413 SPECTRAL DATA START               ← key offset
         8 bytes/point: [absorbance f32 LE, wavelength f32 LE]
         Terminates when wavelength < 100 or > 1100

[remainder] sample name, collection time (metadata strings)
```

### Locating Sample Labels

Sample names (e.g. `0.02mmol`) appear as ASCII strings after the spectral data
within each block. The parser searches each block's byte range with:
```
regex: b"[\d.]+mmol"
```
If not found, the block is labelled `Scan N`.

Collection times appear as:
```
b"Collection Time: DD/M/YYYY H:MM:SS PM"
```

---

## 3. Data Extraction (Python)

```python
import struct, re
import numpy as np

MARKER           = b"TContinuumStor"
DATA_OFFSET      = 0x413   # bytes from block start

def find_blocks(data):
    return [m.start() for m in re.finditer(re.escape(MARKER), data)]

def read_spectrum(data, block_start, n_max=400):
    base = block_start + DATA_OFFSET
    wls, abs_ = [], []
    for i in range(n_max):
        off = base + i * 8
        if off + 8 > len(data): break
        a  = struct.unpack_from("<f", data, off)[0]
        wl = struct.unpack_from("<f", data, off + 4)[0]
        if wl < 100 or wl > 1100: break   # end sentinel
        wls.append(wl); abs_.append(a)
    wls, abs_ = np.array(wls), np.array(abs_)
    idx = np.argsort(wls)
    return wls[idx], abs_[idx]
```

---

## 4. Companion CSV

Cary WinUV can export `.csv` files alongside the `.BSW`. These CSVs carry
richer sample names (e.g. `-0.4-0`, `-0.5`) that are not embedded in the binary.

CSV structure:
```
Row 0: sample names (every other column, blank between = wavelength/abs pair)
Row 1: "Wavelength (nm)", "Abs", "Wavelength (nm)", "Abs", ...
Row 2+: numeric data
```

Parsing:
```python
import pandas as pd

raw   = pd.read_csv("data.csv", header=None)
names = raw.iloc[0].fillna("")
# Pairs: columns (0,1), (2,3), (4,5)... → (wavelength, absorbance) per sample
```

---

## 5. Plotting

### Spectra Overlay

- Colour ramp: `matplotlib.cm.viridis`, 0.1–0.9
- Axes: wavelength (nm) on x; absorbance (AU) on y
- Legend outside-right (bbox_to_anchor=(1.01, 1))
- Spines: top + right hidden; tick direction inward

### Calibration Curve

- X axis: concentration (mmol) — parsed from `NNN.NNmmol` labels
- Y axis: absorbance interpolated at `--cal-wl` (default 700 nm)
- Fit: `numpy.polyfit(conc, abs_, 1)` → linear; display equation + R²
- Markers: scatter (dark), dashed fit line (red)

---

## 6. Known Limitations

| Issue | Detail |
|-------|--------|
| Label extraction | Only `NNN.NNmmol` pattern is auto-parsed; other names need CSV |
| Block count | Assumes all `TContinuumStor` blocks contain spectra; some may be reference/baseline |
| Scan range | Hardcoded sentinel `wl < 100 or > 1100`; works for any visible/NIR range |
| Version bytes | Blocks with version `9`, `7`, `5` tested; other versions untested |
| Calibration | Requires ≥2 blocks with `mmol`-format labels to generate the curve |
