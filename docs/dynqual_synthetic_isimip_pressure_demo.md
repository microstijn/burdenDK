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

## Memory-conscious execution

This demo has been refactored for a memory-conscious execution due to the potentially massive size of the underlying NetCDF datasets.

- The script streams monthly NetCDF slices using `NCDatasets`.
- Robust scaling quantiles are estimated by deterministic sampling over time and space slices.
- Full monthly raw, scaled, and response (Q/F/E) 3D arrays are not materialized simultaneously.
- Pressure memory (`B_state`) is carried by current 2D state arrays only.
- An optional spatial stride can reduce memory usage for demonstration purposes.
- NetCDF outputs are disabled by default for memory safety.

### Examples

PowerShell low-memory first run:

```powershell
$env:TTR_DYNQUAL_SPATIAL_STRIDE="4" 
$env:TTR_DYNQUAL_WRITE_NETCDF="false" 
julia --project=. examples/dynqual_synthetic_isimip_pressure_demo.jl 
```

Fuller-resolution run if memory allows:

```powershell
$env:TTR_DYNQUAL_SPATIAL_STRIDE="1" 
julia --project=. examples/dynqual_synthetic_isimip_pressure_demo.jl 
```

## Comparing baseline and recent tranches

Comparing baseline and recent cluster maps side-by-side can be difficult. The lightweight plotting script (`examples/plot_dynqual_synthetic_isimip_pressure_demo.jl`) generates an explicit comparison panel (`dynqual_baseline_recent_comparison.png`) that highlights changes.

This multi-panel figure adds:
- A changed-regime map showing exactly where transitions occurred.
- A cluster transition heatmap quantifying movement between regimes.
- Feature delta maps highlighting continuous response changes (e.g. amplification $F$ and required adaptive margin $A$).

Additionally, the script generates a clusterwise feature-delta heatmap to summarise what changed in each resulting regime.

These plots are derived entirely from the cache (`dynqual_demo_cache.nc`); no original DynQual files are needed for replotting.

## Replotting without rerunning the DynQual analysis

The analysis script writes a reusable output cache (`dynqual_demo_cache.nc`) and feature metadata. The separate plotting script reads this cache to regenerate figures without rereading the original DynQual NetCDF files. This is extremely useful for changing colormaps, modifying titles, selecting tranches, or adjusting heatmap feature selections without duplicating the heavy pressure analysis and clustering logic.

First, run the heavy analysis to compute vulnerability features and write the cache:

**Using Bash (Linux / macOS):**
```bash
TTR_DYNQUAL_WRITE_CACHE="true" julia --project=. examples/dynqual_synthetic_isimip_pressure_demo.jl
```

**Using PowerShell (Windows):**
```powershell
$env:TTR_DYNQUAL_WRITE_CACHE="true"
julia --project=. examples/dynqual_synthetic_isimip_pressure_demo.jl
```

Then, you can freely regenerate plots only:

**Using Bash (Linux / macOS):**
```bash
julia --project=. examples/plot_dynqual_synthetic_isimip_pressure_demo.jl
```

**Using PowerShell (Windows):**
```powershell
julia --project=. examples/plot_dynqual_synthetic_isimip_pressure_demo.jl
```

If you saved your outputs to an optional custom output directory, specify it:

**Using PowerShell (Windows):**
```powershell
$env:TTR_DYNQUAL_DEMO_OUTPUT_DIR="output/dynqual_synthetic_isimip_pressure_demo"
julia --project=. examples/plot_dynqual_synthetic_isimip_pressure_demo.jl
```

The plotting script exclusively uses the computed arrays in `dynqual_demo_cache.nc` and metadata from `dynqual_feature_metadata.csv`. It does not require access to the original DynQual files.
