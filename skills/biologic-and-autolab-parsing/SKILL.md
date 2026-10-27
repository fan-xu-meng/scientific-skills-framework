---
name: biologic-and-autolab-parsing
description: Parse BioLogic EC-Lab .mpr/.mpt and Autolab NOVA .nox electrochemistry files into time, current, potential, and charge traces. Use when preparing CA, CV, or CP data for plotting or product-efficiency calculations.
metadata:
  version: "1.0.0"
---

# BioLogic and Autolab parsing

Use the original instrument export as the numeric source. Record the filename, instrument format, time and current units, and any missing or truncated segments before combining runs.

1. Identify the experiment type and select the corresponding data columns. Convert time to seconds and current to amperes; preserve the unconverted columns for audit.
2. For charge, use an instrument charge channel when trustworthy and compare it with numerical integration of current over time. State the sign convention.
3. Check that time is monotonic, sample counts match the file metadata, and potential/current traces have plausible ranges. Do not silently fill gaps.
4. Export a tidy table with `time_s`, `current_A`, `potential_V` when available, and `charge_C`, plus a source-file identifier.
5. Plot traces with units and report excluded files or intervals.

The bundled [parser](scripts/ec_biologic_autolab_parser.py) supports a folder of CA sessions:

```bash
python scripts/ec_biologic_autolab_parser.py --data-dir <sessions> --output-dir <output>
```

Review the parser's format assumptions against one representative file from each instrument before processing a batch.

Use the [column and charge audit](references/column-and-charge-audit.md) when selecting channels; the [synthetic trace](examples/synthetic-ca.csv) is only an input-format illustration.

