# DynQual Synthetic ISIMIP Pressure Demo

## Purpose

This demo script (`examples/dynqual_synthetic_isimip_pressure_demo.jl`) illustrates how `TwoTimescaleResilience` can process gridded water-quality data outputs. It is targeted at water-quality modellers who are familiar with DynQual or ISIMIP outputs.

The demo reads real DynQual spatial-temporal patterns over Europe and derives synthetic ISIMIP-style pressure archetypes, such as organic oxygen demand, pathogen exposure, and ionic stress. It routes these pressures through DEB-informed axes using a threshold-free response model, clustering the responses into relative vulnerability regimes.

**Important Note:** The raw NetCDF variables (BOD, pathogen/FC, TDS load, and BOD load) are NOT treated as measured individual chemicals. Instead, they serve as realistic spatial-temporal patterns used to construct synthetic, transparent pressure proxies (e.g., `organic_oxygen_demand_proxy`). The vulnerability regimes produced are threshold-free relative response regimes, not strict safe/unsafe regulatory classes.

## Required Data

The script requires real DynQual NetCDF files, which are not distributed with this package. You must provide the paths to your local files using environment variables:

- `TTR_DYNQUAL_BOD_FILE`: Path to `organic_monthlyAvg_1980_2019.nc`
- `TTR_DYNQUAL_FC_FILE`: Path to `pathogen_monthlyAvg_1980_2019.nc`
- `TTR_DYNQUAL_TDS_FILE`: Path to `TDSload_monthlyAvg_1980_2019.nc`
- `TTR_DYNQUAL_BODLOAD_FILE`: Path to `BODload_monthlyAvg_1980_2019.nc`

If these environment variables are not set, the script will attempt to use a default Windows local path as a convenience fallback. If neither is found, the script will stop with an error.

## Running the Demo

The demo also supports optional plotting via `CairoMakie`. You can control whether plots are generated using the `TTR_DYNQUAL_MAKE_PLOTS` environment variable (defaults to `true`). You can also override the grid subsetting bounding box.

**Using Bash (Linux / macOS):**

```bash
TTR_DYNQUAL_BOD_FILE="/path/to/organic_monthlyAvg_1980_2019.nc" \
TTR_DYNQUAL_FC_FILE="/path/to/pathogen_monthlyAvg_1980_2019.nc" \
TTR_DYNQUAL_TDS_FILE="/path/to/TDSload_monthlyAvg_1980_2019.nc" \
TTR_DYNQUAL_BODLOAD_FILE="/path/to/BODload_monthlyAvg_1980_2019.nc" \
julia --project=. examples/dynqual_synthetic_isimip_pressure_demo.jl
```

**Using PowerShell (Windows):**

```powershell
$env:TTR_DYNQUAL_BOD_FILE="C:\Users\peete074\Downloads\organic_monthlyAvg_1980_2019.nc"
$env:TTR_DYNQUAL_FC_FILE="C:\Users\peete074\Downloads\pathogen_monthlyAvg_1980_2019.nc"
$env:TTR_DYNQUAL_TDS_FILE="C:\Users\peete074\Downloads\TDSload_monthlyAvg_1980_2019.nc"
$env:TTR_DYNQUAL_BODLOAD_FILE="C:\Users\peete074\Downloads\BODload_monthlyAvg_1980_2019.nc"
julia --project=. examples/dynqual_synthetic_isimip_pressure_demo.jl
```

## Generated Outputs

All outputs are generated inside `output/dynqual_synthetic_isimip_pressure_demo/`.

- **CSVs:** Metadata and tables for scaling logic, derived pressure weighting mapping, species selection, tranche definitions, and cluster vulnerability summaries.
- **JSON:** A simulation metadata file (`dynqual_synthetic_isimip_metadata.json`) documenting parameters and constraints.
- **Figures:** If plotting is enabled (and `CairoMakie` is installed), four key figures are generated:
  1. `dynqual_raw_climatology_maps.png`: Real spatial patterns of the DynQual model outputs over Europe.
  2. `dynqual_derived_pressure_layers.png`: Processed proxy layers used in the model chain.
  3. `dynqual_vulnerability_regime_maps.png`: Baseline and recent-decade vulnerability regime maps built via fixed-reference clustering.
  4. `dynqual_regime_explanation_heatmap.png`: Standardised feature means per cluster to explain how vulnerability regimes are differentiated.

## Limitations

- The generated figures provide diagnostic communication value. They do not constitute formal risk assessments.
- This is a demonstration script and does not represent a generalised raster ingestion pipeline for `TwoTimescaleResilience`.
