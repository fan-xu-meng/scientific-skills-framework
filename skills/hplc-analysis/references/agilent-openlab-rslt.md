# Agilent OpenLAB HPLC Result Files

Use this reference when parsing Agilent OpenLAB `.rslt` folders containing RID runs. Confirm the schema against the instrument software version before relying on it.

## Files

| File | Role |
|---|---|
| `.rslt/` | OpenLAB result folder containing sequence, methods, injections, and raw data |
| `.acaml` | XML result-level metadata; maps raw data files to sample names, vial positions, sample type, replicate, and order |
| `.sqx` | zipped sequence table; `SampleListPart/SampleListPart` contains planned samples and vials |
| `.rx` | zipped injection result; `Base/InjectionACAML` contains integrated peaks and compound assignments |
| `.dx` | zipped injection raw data; contains `injection.acmd` and binary channel traces such as RID1A `.CH` |
| `hplc.xlsx` | exported OpenLAB report; useful for audit but may be incomplete or difficult for openpyxl to parse |

## `.acaml` Metadata

Parse `InjectionMetaData` nodes. Important attributes:

- `RawDataFileName`: an acquisition timestamp ending in `.dx`
- `SampleName`
- `SampleType`: `Blank`, `Calibration`, or `Sample`
- `VialNumber`
- `InjectionAcqDateTime`

Important child nodes:

- `ReplicateNumber` with `val`
- `SampleOrderNumber` with `val`
- `SampleInjectionsCount` with `val`

Use `RawDataFileName.replace(".dx", ".rx")` to find the paired injection result.

## `.sqx` Sequence

Open as zip and read `SampleListPart/SampleListPart`. The XML namespace is verbose; parse by local tag name.

For each `anyType` entry, useful fields include:

- `SampleName`
- `SampleType`
- `Vial`
- `InjPerVial`
- `Volume`
- `AcquisitionMethod`
- `AnalysisMethod`

Use `.sqx` to understand planned sequence layout. Prefer `.acaml` for actual injections, because it maps directly to raw files.

## `.rx` Peak Data

Open as zip and parse `Base/InjectionACAML`.

Peak nodes:

- local name `Peak`
- attribute `id`
- children such as `RetentionTime`, `Area`, `Height`, `BeginTime`, `EndTime`, `BaselineStart`, `BaselineEnd`, `WidthBase`
- values may be in child text or `val` attribute

Compound mapping:

- local name `InjectionCompound`
- child `CompoundName`
- nested `Peak_ID` with peak `id`

Join `InjectionCompound` to `Peak` by `Peak_ID`. If no `Peak_ID` exists, the compound was not integrated/detected in that injection.

## `.dx` RID Raw Trace

Open as zip.

1. Parse `injection.acmd`.
2. Find a `Signal` where:
   - `ChannelName == RID1A`
   - `IsIntegrable == true`
3. Read:
   - `TraceId`
   - `TimeStart`
   - `TimeEnd`
   - `NumberOfValues`
4. Read binary channel file `<TraceId>.CH`.
5. Inspect the `.CH` header and numeric encoding for the instrument version; a 6144-byte header followed by little-endian float64 values is one observed layout, not a universal specification.
6. Convert Agilent time values from ms to min:
   - `time_start_min = TimeStart / 60000`
   - `time_end_min = TimeEnd / 60000`
7. Generate the time axis with `numpy.linspace(time_start_min, time_end_min, NumberOfValues)`.

## Workbook Parsing Fallback

Some exported `.xlsx` files contain styles that make `openpyxl` fail. In that case:

1. Treat `.xlsx` as zip.
2. Read `xl/sharedStrings.xml`.
3. Read `xl/worksheets/sheet.xml` or `sheet1.xml`.
4. Resolve cells with `t="s"` through shared strings.

This is enough to recover OpenLAB summary tables with columns such as name, vial, RT, amount, concentration, and area.

## Plot Conventions

- Full raw trace: plot all RID1A data and mark expected retention windows from calibration points.
- Product-window checks: for each calibrated product, plot a zoomed region around the expected calibration integration window.
- Detected peak plot: show actual integration begin/end, baseline, filled integrated area, RT, area, and calculated concentration.
