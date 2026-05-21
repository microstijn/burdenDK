# ==============================================================================
# nc_monthly_longterm_isimip_moa_deb_inspection.jl
#
# Manual inspection script for the ISIMIP Water Quality -> Mode-of-Action ->
# reduced-DEB response-operator pipeline.
#
# This is a manual inspection script.
# It is intentionally not included in automated tests because it can require
# local NetCDF files. If those files are not available, the script can generate
# a deterministic synthetic ISIMIP-style estuary dataset for all seven canonical
# variables so the full pipeline remains inspectable.
#
# The synthetic estuary is designed to produce meaningful spatial and temporal
# variation in amplification F:
#   - river-to-sea gradient;
#   - salinity intrusion;
#   - urban wastewater hotspot;
#   - upstream nutrient source;
#   - port/industrial hotspot;
#   - shoreline plastic accumulation;
#   - bloom/eutrophication zone;
#   - hotspot expansion and intensification over time.
#
# Style philosophy:
#   - user settings at the top;
#   - helper functions below settings;
#   - load one example month first;
#   - keep intermediate objects in top-level variables;
#   - then loop through selected months;
#   - store intermediate outputs in dictionaries;
#   - summary_rows remains a Vector{NamedTuple};
#   - DataFrames is optional only;
#   - plotting is optional and controlled by ENV["TTR_EXAMPLE_PLOTS"].
# ==============================================================================
#!/usr/bin/env julia

using Pkg
try
    Pkg.activate(joinpath(@__DIR__, ".."))
catch err
    @warn "Could not activate project" exception=(err, catch_backtrace())
end

using TwoTimescaleResilience
using NCDatasets
using Statistics
using CairoMakie

# ==============================================================================
# USER SETTINGS
# ==============================================================================

# Edit these paths and variable names to match your local ISIMIP Water Quality files.
# If the files do not exist and use_synthetic_if_files_missing=true, the script
# will use a deterministic synthetic dataset instead.
WT_file       = raw"C:\path\to\watertemp_monthlyAvg_1980_2019.nc"
BOD_file      = raw"C:\path\to\bod_monthlyAvg_1980_2019.nc"
TDS_file      = raw"C:\path\to\tds_monthlyAvg_1980_2019.nc"
FC_file       = raw"C:\path\to\fc_monthlyAvg_1980_2019.nc"
Nutrient_file = raw"C:\path\to\nutrient_monthlyAvg_1980_2019.nc"
Chemical_file = raw"C:\path\to\chemical_monthlyAvg_1980_2019.nc"
Plastic_file  = raw"C:\path\to\plastic_monthlyAvg_1980_2019.nc"

WT_var       = "watertemp"
BOD_var      = "bod"
TDS_var      = "tds"
FC_var       = "fc"
Nutrient_var = "nutrient"
Chemical_var = "chemical"
Plastic_var  = "plastic"

# If real files are missing, use a complete deterministic synthetic dataset.
use_synthetic_if_files_missing = true

# Synthetic options:
#   :estuary gives a larger, spatially structured estuary example with multiple hotspots.
#   :debug9 gives a small 9x9 debugging example.
synthetic_mock_type = :estuary
synthetic_n = 41
# synthetic_mock_type = :debug9
# synthetic_n = 9

# Synthetic layers are already in a 0..1 stress scale. If true, the script does
# not re-normalise synthetic layers year-by-year, preserving temporal trends.
skip_normalisation_for_synthetic = true

# If real files exist but some variables are missing, fill those missing variables
# with zero grids after a reference grid has been loaded.
fill_missing_variables_with_zero = true

start_year = 1980
end_year   = 2019
base_year  = 1980

# Inspect only selected years/months first.
selected_years = [1980, 1990, 2000, 2010, 2019]
selected_month = 7

# Optional crop. Both formats are accepted by helper functions:
#   bbox = (-12.0, 35.0, 34.0, 72.0)
#   bbox = (lon=(-12.0, 35.0), lat=(34.0, 72.0))
# For synthetic data, bbox is used to define mock lon/lat ranges.
# Estuary-like default using Wadden/Rhine-Meuse-style coordinates:
bbox = (lon = (3.0, 5.2), lat = (51.0, 53.2))
# bbox = nothing

# To process every month between start_year and end_year, set to true.
process_all_months = false

# Select base target profile. Alternatives include:
#   selected_profile_base = aquatic_invertebrate_profile()
#   selected_profile_base = fish_profile()
#   selected_profile_base = bivalve_profile()
#   selected_profile_base = amphibian_profile()
#   selected_profile_base = bird_profile()
#   selected_profile_base = small_mammal_profile()
#   selected_profile_base = human_profile()
selected_profile_base = fish_profile()

# Demo-sensitive profile makes synthetic estuary contrasts easier to inspect.
# Set false if you want to use the package default profile exactly.
use_demo_sensitive_profile = true

# Plotting is disabled by default to keep smoke tests/headless runs robust.
# To enable manually before running:
#   ENV["TTR_EXAMPLE_PLOTS"] = "true"
do_plots = get(ENV, "TTR_EXAMPLE_PLOTS", "true") == "true"

output_dir = joinpath(@__DIR__, "..", "output", "nc_monthly_longterm_isimip_moa_deb_inspection")
mkpath(output_dir)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

finite_values(x) = filter(isfinite, vec(collect(skipmissing(x))))

function robust_normalise(x; lower_q=0.02, upper_q=0.98, log_transform=false)
    fv = finite_values(x)
    if isempty(fv)
        return fill(NaN, size(x))
    end

    fv_work = log_transform ? log10.(1.0 .+ max.(fv, 0.0)) : fv
    lo = quantile(fv_work, lower_q)
    hi = quantile(fv_work, upper_q)

    out = similar(x, Float64)

    if (hi - lo == 0) || isapprox(hi, lo; atol=1e-12)
        @warn "Normalisation range is zero. Returning zeros for finite cells."
        for i in eachindex(x)
            out[i] = isfinite(x[i]) ? 0.0 : NaN
        end
        return out
    end

    for i in eachindex(x)
        if !isfinite(x[i])
            out[i] = NaN
        else
            val = Float64(x[i])
            if log_transform
                val = log10(1.0 + max(val, 0.0))
            end
            out[i] = clamp((val - lo) / (hi - lo), 0.0, 1.0)
        end
    end

    return out
end

replace_missing_with_nan(A) = Matrix{Float64}(coalesce.(A, NaN))

function bbox_ranges(bbox)
    if bbox === nothing
        return nothing
    elseif bbox isa NamedTuple
        lon_min, lon_max = bbox.lon
        lat_min, lat_max = bbox.lat
        return (lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)
    elseif bbox isa Tuple && length(bbox) == 4
        lon_min, lon_max, lat_min, lat_max = bbox
        return (lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)
    else
        error("bbox must be nothing, a 4-tuple (lon_min, lon_max, lat_min, lat_max), or a NamedTuple with lon and lat ranges")
    end
end

function load_nc_slice(path, varname; time_index, bbox=nothing)
    ds = NCDataset(path, "r")
    try
        lon_name = "longitude" in keys(ds) ? "longitude" : ("lon" in keys(ds) ? "lon" : error("No longitude/lon variable found"))
        lat_name = "latitude" in keys(ds) ? "latitude" : ("lat" in keys(ds) ? "lat" : error("No latitude/lat variable found"))

        lon = ds[lon_name][:]
        lat = ds[lat_name][:]

        println("Longitude range in file: ", minimum(lon), " to ", maximum(lon))

        if !(varname in keys(ds))
            error("Variable $varname not found in $path")
        end

        br = bbox_ranges(bbox)
        if br !== nothing
            lon_min, lon_max, lat_min, lat_max = br.lon_min, br.lon_max, br.lat_min, br.lat_max

            if minimum(lon) >= 0 && lon_min < 0
                # File is 0..360 but bbox contains negative longitudes.
                lon_min_mod = mod(lon_min, 360)
                lon_max_mod = mod(lon_max, 360)
                if lon_min_mod > lon_max_mod
                    lon_idx = findall(x -> (x >= lon_min_mod) || (x <= lon_max_mod), lon)
                else
                    lon_idx = findall(x -> lon_min_mod <= x <= lon_max_mod, lon)
                end
            else
                lon_idx = findall(x -> lon_min <= x <= lon_max, lon)
            end
            lat_idx = findall(x -> lat_min <= x <= lat_max, lat)
        else
            lon_idx = collect(1:length(lon))
            lat_idx = collect(1:length(lat))
        end

        data_slice = ds[varname][lon_idx, lat_idx, time_index]
        lon_selected = lon[lon_idx]
        lat_selected = lat[lat_idx]

        layer = replace_missing_with_nan(data_slice)
        return layer, lon_selected, lat_selected
    finally
        close(ds)
    end
end

time_index_from_year_month(year, month; base_year=1980) = (year - base_year) * 12 + month

function sort_for_plot(layer, lon, lat)
    layer_sorted = copy(layer)
    lon_sorted = copy(lon)
    lat_sorted = copy(lat)

    if length(lon) > 1 && lon[1] > lon[end]
        lon_sorted = reverse(lon)
        layer_sorted = reverse(layer_sorted, dims=1)
    end

    if length(lat) > 1 && lat[1] > lat[end]
        lat_sorted = reverse(lat)
        layer_sorted = reverse(layer_sorted, dims=2)
    end

    return layer_sorted, lon_sorted, lat_sorted
end

function plot_geo_grid(layer, lon, lat; title, filename, colormap=:viridis, colorrange=nothing)
    layer_s, lon_s, lat_s = sort_for_plot(layer, lon, lat)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], title=title, xlabel="Longitude", ylabel="Latitude")

    hm = if colorrange === nothing
        heatmap!(
            ax,
            lon_s,
            lat_s,
            layer_s;
            colormap=colormap,
            nan_color=:transparent
        )
    else
        heatmap!(
            ax,
            lon_s,
            lat_s,
            layer_s;
            colormap=colormap,
            colorrange=colorrange,
            nan_color=:transparent
        )
    end

    Colorbar(fig[1, 2], hm)

    save(filename, fig, px_per_unit=2)
    return fig
end

function monthly_summary(grid)
    fv = finite_values(grid)
    if isempty(fv)
        return (min=NaN, mean=NaN, median=NaN, max=NaN, q90=NaN, q95=NaN, q99=NaN)
    end
    return (
        min = minimum(fv),
        mean = mean(fv),
        median = median(fv),
        max = maximum(fv),
        q90 = quantile(fv, 0.90),
        q95 = quantile(fv, 0.95),
        q99 = quantile(fv, 0.99)
    )
end

function make_demo_sensitive_profile(profile)
    return SpeciesProfile(
        name = profile.name * "_demo_sensitive",
        exposure_filter = profile.exposure_filter,
        moa_mapping = profile.moa_mapping,
        moa_deb_mapping = profile.moa_deb_mapping,
        deb_params = DEBAxisParams(
            A0 = 1.0,
            alpha_axes = (0.35, 0.45, 0.25, 0.15),
            lambda_min = 0.04,
            lambda_max = 1.0,
            KA = 0.30,
            recovery_axes = (0.20, 1.20, 0.25, 0.10),
            use_axis_recovery_penalty = true,
            use_buffer_recovery_factor = false,
            beta_Z = 0.0
        ),
        buffer_params = profile.buffer_params,
        description = profile.description * " Demo-sensitive profile for synthetic estuary inspection."
    )
end

selected_profile = use_demo_sensitive_profile ? make_demo_sensitive_profile(selected_profile_base) : selected_profile_base

# ==============================================================================
# SYNTHETIC ISIMIP-STYLE MOCK DATASETS
# ==============================================================================

"""
    synthetic_isimip_9x9_layers(; year=2010, month=7, n=9)

Small deterministic mock dataset for debugging.
"""
function synthetic_isimip_9x9_layers(; year=2010, month=7, n=9)
    x = collect(range(0.0, 1.0, length=n))
    y = collect(range(0.0, 1.0, length=n))
    trend = clamp((year - 1980) / (2019 - 1980), 0.0, 1.0)
    seasonal = 0.5 + 0.5 * sin(2 * pi * (month - 1) / 12)

    WT       = zeros(Float64, n, n)
    BOD      = zeros(Float64, n, n)
    TDS      = zeros(Float64, n, n)
    FC       = zeros(Float64, n, n)
    Nutrient = zeros(Float64, n, n)
    Chemical = zeros(Float64, n, n)
    Plastic  = zeros(Float64, n, n)

    for i in 1:n
        for j in 1:n
            xx = x[i]
            yy = y[j]
            centre_hotspot   = exp(-((xx - 0.55)^2 + (yy - 0.55)^2) / 0.045)
            pathogen_hotspot = exp(-((xx - 0.30)^2 + (yy - 0.70)^2) / 0.030)
            chemical_hotspot = exp(-((xx - 0.75)^2 + (yy - 0.35)^2) / 0.040)

            WT[i, j]       = 0.20 + 0.35 * xx + 0.25 * trend + 0.15 * seasonal
            BOD[i, j]      = 0.10 + 0.55 * centre_hotspot + 0.20 * trend
            TDS[i, j]      = 0.10 + 0.45 * xx + 0.25 * (1.0 - yy)
            FC[i, j]       = 0.05 + 0.65 * pathogen_hotspot + 0.15 * trend
            Nutrient[i, j] = 0.10 + 0.35 * yy + 0.30 * centre_hotspot
            Chemical[i, j] = 0.05 + 0.55 * chemical_hotspot + 0.20 * trend
            Plastic[i, j]  = 0.05 + 0.35 * xx + 0.35 * (1.0 - yy)
        end
    end

    for L in (WT, BOD, TDS, FC, Nutrient, Chemical, Plastic)
        clamp!(L, 0.0, 1.0)
    end

    br = bbox_ranges(bbox)
    lon = br === nothing ? collect(range(100.0, 160.0, length=n)) : collect(range(br.lon_min, br.lon_max, length=n))
    lat = br === nothing ? collect(range(-45.0, 5.0, length=n)) : collect(range(br.lat_min, br.lat_max, length=n))

    return (WT=WT, BOD=BOD, TDS=TDS, FC=FC, Nutrient=Nutrient, Chemical=Chemical, Plastic=Plastic, lon=lon, lat=lat, synthetic=true)
end

"""
    synthetic_estuary_isimip_layers(; year=2010, month=7, n=41)

Create a deterministic synthetic estuary-like ISIMIP Water Quality dataset.

Design:
    - x direction = upstream river to coastal sea;
    - y direction = lateral estuary width;
    - salinity/TDS increases seaward;
    - nutrients/pathogens/BOD originate upstream and from an urban mid-estuary point source;
    - chemicals/plastics originate from port/industrial zones near the mouth;
    - stressors intensify and hotspots expand over time;
    - salt intrusion moves upstream over time;
    - interactions are encouraged by overlapping thermal, BOD, nutrient, and pathogen zones.
"""
function synthetic_estuary_isimip_layers(; year=2010, month=7, n=41)
    x = collect(range(0.0, 1.0, length=n))  # 0 = upstream, 1 = sea
    y = collect(range(0.0, 1.0, length=n))  # lateral estuary coordinate

    trend = clamp((year - 1980) / (2019 - 1980), 0.0, 1.0)
    seasonal = 0.5 + 0.5 * sin(2 * pi * (month - 1) / 12)

    WT       = zeros(Float64, n, n)
    BOD      = zeros(Float64, n, n)
    TDS      = zeros(Float64, n, n)
    FC       = zeros(Float64, n, n)
    Nutrient = zeros(Float64, n, n)
    Chemical = zeros(Float64, n, n)
    Plastic  = zeros(Float64, n, n)

    urban_sigma = 0.020 + 0.035 * trend
    port_sigma  = 0.025 + 0.030 * trend
    bloom_sigma = 0.035 + 0.040 * trend
    salt_intrusion_shift = 0.10 * trend

    for i in 1:n
        for j in 1:n
            xx = x[i]
            yy = y[j]

            channel = exp(-((yy - 0.50)^2) / 0.035)
            shallow_north = exp(-((yy - 0.82)^2) / 0.020)
            shallow_south = exp(-((yy - 0.18)^2) / 0.020)

            urban_hotspot = exp(-((xx - 0.38)^2 + (yy - 0.55)^2) / urban_sigma)
            upstream_source = exp(-((xx - 0.12)^2) / 0.035) * (0.5 + 0.5 * channel)
            port_hotspot = exp(-((xx - 0.78)^2 + (yy - 0.42)^2) / port_sigma)
            shore_accumulation = 0.5 * shallow_north + 0.5 * shallow_south
            bloom_zone = exp(-((xx - 0.52)^2 + (yy - 0.50)^2) / bloom_sigma)
            river_decay = exp(-2.2 * xx)
            salt_wedge = clamp((xx + salt_intrusion_shift - 0.30) / 0.70, 0.0, 1.0)

            WT[i, j] =
                0.18 +
                0.18 * seasonal +
                0.18 * trend +
                0.18 * (shallow_north + shallow_south) +
                0.10 * bloom_zone

            BOD[i, j] =
                0.06 +
                0.35 * urban_hotspot * (0.8 + 0.6 * trend) +
                0.20 * upstream_source * river_decay +
                0.25 * bloom_zone * (0.5 + trend)

            TDS[i, j] =
                0.05 +
                0.70 * salt_wedge +
                0.10 * port_hotspot

            FC[i, j] =
                0.03 +
                0.55 * urban_hotspot * (0.7 + 0.8 * trend) +
                0.25 * upstream_source * river_decay

            Nutrient[i, j] =
                0.06 +
                0.40 * upstream_source * (0.8 + 0.4 * trend) +
                0.25 * urban_hotspot +
                0.25 * bloom_zone * (0.5 + trend)

            Chemical[i, j] =
                0.04 +
                0.55 * port_hotspot * (0.7 + 0.8 * trend) +
                0.15 * urban_hotspot

            Plastic[i, j] =
                0.04 +
                0.35 * port_hotspot +
                0.30 * shore_accumulation * (0.6 + 0.8 * trend) +
                0.20 * salt_wedge * trend
        end
    end

    for L in (WT, BOD, TDS, FC, Nutrient, Chemical, Plastic)
        clamp!(L, 0.0, 1.0)
    end

    br = bbox_ranges(bbox)
    lon = br === nothing ? collect(range(3.0, 5.0, length=n)) : collect(range(br.lon_min, br.lon_max, length=n))
    lat = br === nothing ? collect(range(51.0, 53.0, length=n)) : collect(range(br.lat_min, br.lat_max, length=n))

    return (WT=WT, BOD=BOD, TDS=TDS, FC=FC, Nutrient=Nutrient, Chemical=Chemical, Plastic=Plastic, lon=lon, lat=lat, synthetic=true)
end

function synthetic_layers(; year, month)
    if synthetic_mock_type == :estuary
        return synthetic_estuary_isimip_layers(year=year, month=month, n=synthetic_n)
    elseif synthetic_mock_type == :debug9
        return synthetic_isimip_9x9_layers(year=year, month=month, n=synthetic_n)
    else
        error("Unknown synthetic_mock_type=$(synthetic_mock_type). Use :estuary or :debug9.")
    end
end

# ==============================================================================
# LOADING AND NORMALISATION LOGIC
# ==============================================================================

canonical_files = (WT_file, BOD_file, TDS_file, FC_file, Nutrient_file, Chemical_file, Plastic_file)
canonical_vars  = (WT_var, BOD_var, TDS_var, FC_var, Nutrient_var, Chemical_var, Plastic_var)

function real_files_available()
    return all(isfile, canonical_files)
end

function load_or_mock_all_layers(year, month)
    time_idx = time_index_from_year_month(year, month; base_year=base_year)

    if real_files_available()
        println("Loading real NetCDF layers for $year-$month")

        loaded_layers = Vector{Matrix{Float64}}(undef, 7)
        lon_ref = nothing
        lat_ref = nothing
        ref_shape = nothing

        for idx in 1:7
            file = canonical_files[idx]
            var = canonical_vars[idx]

            try
                layer, lon, lat = load_nc_slice(file, var; time_index=time_idx, bbox=bbox)
                loaded_layers[idx] = layer

                if lon_ref === nothing
                    lon_ref = lon
                    lat_ref = lat
                    ref_shape = size(layer)
                else
                    if size(layer) != ref_shape
                        error("Layer for variable $var has shape $(size(layer)), expected $ref_shape")
                    end
                end
            catch err
                if fill_missing_variables_with_zero && ref_shape !== nothing
                    @warn "Variable/file for $(var) could not be loaded. Filling with zeros." exception=(err, catch_backtrace())
                    loaded_layers[idx] = fill(0.0, ref_shape)
                else
                    rethrow(err)
                end
            end
        end

        return (
            WT = loaded_layers[1],
            BOD = loaded_layers[2],
            TDS = loaded_layers[3],
            FC = loaded_layers[4],
            Nutrient = loaded_layers[5],
            Chemical = loaded_layers[6],
            Plastic = loaded_layers[7],
            lon = lon_ref,
            lat = lat_ref,
            synthetic = false
        )
    elseif use_synthetic_if_files_missing
        println("Real NetCDF files not found. Using synthetic $(synthetic_mock_type) ISIMIP-style mock dataset for $year-$month.")
        return synthetic_layers(year=year, month=month)
    else
        files_vec = collect(canonical_files)
        missing_files = files_vec[.!isfile.(files_vec)]
        error("One or more NetCDF files are missing and use_synthetic_if_files_missing=false: $(missing_files)")
    end
end

function normalise_loaded_layers(loaded)
    if loaded.synthetic && skip_normalisation_for_synthetic
        # Synthetic layers are deliberately already scaled to 0..1. Do not
        # re-normalise per year, otherwise temporal intensification is removed.
        return [
            loaded.WT,
            loaded.BOD,
            loaded.TDS,
            loaded.FC,
            loaded.Nutrient,
            loaded.Chemical,
            loaded.Plastic
        ]
    else
        return [
            robust_normalise(loaded.WT),
            robust_normalise(loaded.BOD),
            robust_normalise(loaded.TDS),
            robust_normalise(loaded.FC),
            robust_normalise(loaded.Nutrient),
            robust_normalise(loaded.Chemical),
            robust_normalise(loaded.Plastic)
        ]
    end
end

# ==============================================================================
# MAIN WORKFLOW -- STEP 1: LOAD EXAMPLE MONTH
# ==============================================================================

println("\n--- ISIMIP MoA DEB NetCDF Inspection Script ---")
println("Profile: $(selected_profile.name)")
println("Synthetic fallback enabled: $use_synthetic_if_files_missing")
println("Synthetic mock type: $synthetic_mock_type, n=$synthetic_n")
println("Skip normalisation for synthetic layers: $skip_normalisation_for_synthetic")

example_year = selected_years[1]
example_month = selected_month
example_time_index = time_index_from_year_month(example_year, example_month; base_year=base_year)

println("\n--- Step 1: Loading Example Month (Year: $example_year, Month: $example_month, Time Index: $example_time_index) ---")

example_layers = load_or_mock_all_layers(example_year, example_month)

WT_raw       = example_layers.WT
BOD_raw      = example_layers.BOD
TDS_raw      = example_layers.TDS
FC_raw       = example_layers.FC
Nutrient_raw = example_layers.Nutrient
Chemical_raw = example_layers.Chemical
Plastic_raw  = example_layers.Plastic
lon          = example_layers.lon
lat          = example_layers.lat

println("Raster size: ", size(WT_raw))
println("Longitude range: ", minimum(lon), " to ", maximum(lon))
println("Latitude range: ", minimum(lat), " to ", maximum(lat))

# ==============================================================================
# STEP 2: NORMALISE EXAMPLE MONTH
# ==============================================================================

println("\n--- Step 2: Preparing Example Month Layers ---")
println(example_layers.synthetic && skip_normalisation_for_synthetic ? "Synthetic layers are already 0..1; skipping per-year normalisation." : "Applying robust normalisation.")

layers_norm = normalise_loaded_layers(example_layers)

WT_norm       = layers_norm[1]
BOD_norm      = layers_norm[2]
TDS_norm      = layers_norm[3]
FC_norm       = layers_norm[4]
Nutrient_norm = layers_norm[5]
Chemical_norm = layers_norm[6]
Plastic_norm  = layers_norm[7]

# ==============================================================================
# STEP 3: RUN ISIMIP -> MOA -> DEB PIPELINE
# ==============================================================================

println("\n--- Step 3: Running ISIMIP MoA DEB Pipeline ---")

result = isimip_deb_pipeline_grid(
    layers_norm,
    selected_profile.exposure_filter,
    selected_profile.moa_mapping,
    selected_profile.moa_deb_mapping,
    selected_profile.deb_params
)

effective_layers = result.effective_layers
modes = result.modes
axes = result.axes
Zgrid = result.Z
Agrid = result.A
lambdagrid = result.lambda
Fgrid = result.amplification

println("Successfully ran pipeline on example month.")

# ==============================================================================
# STEP 4: SUMMARIES FOR EXAMPLE MONTH
# ==============================================================================

println("\n--- Step 4: Summaries for Example Month ---")
println("WT_norm:        ", monthly_summary(WT_norm))
println("BOD_norm:       ", monthly_summary(BOD_norm))
println("TDS_norm:       ", monthly_summary(TDS_norm))
println("FC_norm:        ", monthly_summary(FC_norm))
println("Nutrient_norm:  ", monthly_summary(Nutrient_norm))
println("Chemical_norm:  ", monthly_summary(Chemical_norm))
println("Plastic_norm:   ", monthly_summary(Plastic_norm))
println("mode thermal:   ", monthly_summary(modes.thermal))
println("mode oxygen:    ", monthly_summary(modes.oxygen))
println("mode osmotic:   ", monthly_summary(modes.osmotic))
println("mode immune:    ", monthly_summary(modes.immune))
println("mode eutro:     ", monthly_summary(modes.eutrophication))
println("mode toxic:     ", monthly_summary(modes.toxic))
println("mode feeding:   ", monthly_summary(modes.feeding))
println("mode physical:  ", monthly_summary(modes.physical))
println("assimilation:   ", monthly_summary(axes.assimilation))
println("maintenance:    ", monthly_summary(axes.maintenance))
println("growth:         ", monthly_summary(axes.growth))
println("reproduction:   ", monthly_summary(axes.reproduction))
println("Agrid:          ", monthly_summary(Agrid))
println("lambdagrid:     ", monthly_summary(lambdagrid))
println("Fgrid:          ", monthly_summary(Fgrid))

# ==============================================================================
# STEP 5: PLOT EXAMPLE MONTH IF REQUESTED
# ==============================================================================

if do_plots
    println("\n--- Step 5: Plotting Example Month ---")
    plot_geo_grid(WT_norm, lon, lat, title="WT Norm ($example_year-$example_month)", filename=joinpath(output_dir, "WT_norm_$(example_year)_$(example_month).png"))
    plot_geo_grid(BOD_norm, lon, lat, title="BOD Norm ($example_year-$example_month)", filename=joinpath(output_dir, "BOD_norm_$(example_year)_$(example_month).png"))
    plot_geo_grid(TDS_norm, lon, lat, title="TDS Norm ($example_year-$example_month)", filename=joinpath(output_dir, "TDS_norm_$(example_year)_$(example_month).png"))
    plot_geo_grid(FC_norm, lon, lat, title="FC Norm ($example_year-$example_month)", filename=joinpath(output_dir, "FC_norm_$(example_year)_$(example_month).png"))
    plot_geo_grid(Nutrient_norm, lon, lat, title="Nutrient Norm ($example_year-$example_month)", filename=joinpath(output_dir, "Nutrient_norm_$(example_year)_$(example_month).png"))
    plot_geo_grid(Chemical_norm, lon, lat, title="Chemical Norm ($example_year-$example_month)", filename=joinpath(output_dir, "Chemical_norm_$(example_year)_$(example_month).png"))
    plot_geo_grid(Plastic_norm, lon, lat, title="Plastic Norm ($example_year-$example_month)", filename=joinpath(output_dir, "Plastic_norm_$(example_year)_$(example_month).png"))
    plot_geo_grid(modes.oxygen, lon, lat, title="Oxygen Mode ($example_year-$example_month)", filename=joinpath(output_dir, "mode_oxygen_$(example_year)_$(example_month).png"))
    plot_geo_grid(modes.eutrophication, lon, lat, title="Eutrophication Mode ($example_year-$example_month)", filename=joinpath(output_dir, "mode_eutrophication_$(example_year)_$(example_month).png"))
    plot_geo_grid(axes.assimilation, lon, lat, title="Assimilation Axis ($example_year-$example_month)", filename=joinpath(output_dir, "axis_assimilation_$(example_year)_$(example_month).png"))
    plot_geo_grid(axes.maintenance, lon, lat, title="Maintenance Axis ($example_year-$example_month)", filename=joinpath(output_dir, "axis_maintenance_$(example_year)_$(example_month).png"))
    plot_geo_grid(Agrid, lon, lat, title="Adaptive Margin ($example_year-$example_month)", filename=joinpath(output_dir, "Agrid_$(example_year)_$(example_month).png"))
    plot_geo_grid(lambdagrid, lon, lat, title="Restoring Force ($example_year-$example_month)", filename=joinpath(output_dir, "lambdagrid_$(example_year)_$(example_month).png"))
    plot_geo_grid(Fgrid, lon, lat, title="Amplification Factor ($example_year-$example_month)", filename=joinpath(output_dir, "Fgrid_$(example_year)_$(example_month).png"))
else
    println("\n--- Step 5: Plotting skipped. Set ENV[\"TTR_EXAMPLE_PLOTS\"] = \"true\" to enable plots. ---")
end

# ==============================================================================
# LONG-TERM MONTHLY ANALYSIS
# ==============================================================================

println("\n--- Long-term monthly analysis ---")

if process_all_months
    selected_times = [(y, m) for y in start_year:end_year for m in 1:12]
else
    selected_times = [(year, selected_month) for year in selected_years]
end

raw_layers_by_time = Dict{String, Any}()
norm_layers_by_time = Dict{String, Any}()
effective_layers_by_time = Dict{String, Any}()
modes_by_time = Dict{String, Any}()
axes_by_time = Dict{String, Any}()
Z_by_time = Dict{String, Any}()
A_by_time = Dict{String, Any}()
lambda_by_time = Dict{String, Any}()
amplification_by_time = Dict{String, Any}()
summaries_by_time = Dict{String, Any}()

summary_rows = NamedTuple[]

for (yr, mo) in selected_times
    t_idx = time_index_from_year_month(yr, mo; base_year=base_year)
    key = string(yr, "-", lpad(mo, 2, "0"))

    println("Processing $key (index $t_idx)")

    loaded = load_or_mock_all_layers(yr, mo)
    lyrs_norm = normalise_loaded_layers(loaded)

    res_t = isimip_deb_pipeline_grid(
        lyrs_norm,
        selected_profile.exposure_filter,
        selected_profile.moa_mapping,
        selected_profile.moa_deb_mapping,
        selected_profile.deb_params
    )

    raw_layers_by_time[key] = (
        WT = loaded.WT,
        BOD = loaded.BOD,
        TDS = loaded.TDS,
        FC = loaded.FC,
        Nutrient = loaded.Nutrient,
        Chemical = loaded.Chemical,
        Plastic = loaded.Plastic
    )

    norm_layers_by_time[key] = (
        WT = lyrs_norm[1],
        BOD = lyrs_norm[2],
        TDS = lyrs_norm[3],
        FC = lyrs_norm[4],
        Nutrient = lyrs_norm[5],
        Chemical = lyrs_norm[6],
        Plastic = lyrs_norm[7]
    )

    effective_layers_by_time[key] = res_t.effective_layers
    modes_by_time[key] = res_t.modes
    axes_by_time[key] = res_t.axes
    Z_by_time[key] = res_t.Z
    A_by_time[key] = res_t.A
    lambda_by_time[key] = res_t.lambda
    amplification_by_time[key] = res_t.amplification

    sums = (
        WT = monthly_summary(lyrs_norm[1]),
        BOD = monthly_summary(lyrs_norm[2]),
        TDS = monthly_summary(lyrs_norm[3]),
        FC = monthly_summary(lyrs_norm[4]),
        Nutrient = monthly_summary(lyrs_norm[5]),
        Chemical = monthly_summary(lyrs_norm[6]),
        Plastic = monthly_summary(lyrs_norm[7]),
        thermal = monthly_summary(res_t.modes.thermal),
        oxygen = monthly_summary(res_t.modes.oxygen),
        osmotic = monthly_summary(res_t.modes.osmotic),
        immune = monthly_summary(res_t.modes.immune),
        eutrophication = monthly_summary(res_t.modes.eutrophication),
        toxic = monthly_summary(res_t.modes.toxic),
        feeding = monthly_summary(res_t.modes.feeding),
        physical = monthly_summary(res_t.modes.physical),
        assimilation = monthly_summary(res_t.axes.assimilation),
        maintenance = monthly_summary(res_t.axes.maintenance),
        growth = monthly_summary(res_t.axes.growth),
        reproduction = monthly_summary(res_t.axes.reproduction),
        A = monthly_summary(res_t.A),
        lambda = monthly_summary(res_t.lambda),
        amplification = monthly_summary(res_t.amplification)
    )

    summaries_by_time[key] = sums

    push!(summary_rows, (
        year = yr,
        month = mo,
        key = key,
        WT_mean = sums.WT.mean,
        BOD_mean = sums.BOD.mean,
        TDS_mean = sums.TDS.mean,
        FC_mean = sums.FC.mean,
        Nutrient_mean = sums.Nutrient.mean,
        Chemical_mean = sums.Chemical.mean,
        Plastic_mean = sums.Plastic.mean,
        thermal_mean = sums.thermal.mean,
        oxygen_mean = sums.oxygen.mean,
        osmotic_mean = sums.osmotic.mean,
        immune_mean = sums.immune.mean,
        eutrophication_mean = sums.eutrophication.mean,
        toxic_mean = sums.toxic.mean,
        feeding_mean = sums.feeding.mean,
        physical_mean = sums.physical.mean,
        assimilation_mean = sums.assimilation.mean,
        maintenance_mean = sums.maintenance.mean,
        growth_mean = sums.growth.mean,
        reproduction_mean = sums.reproduction.mean,
        A_mean = sums.A.mean,
        lambda_mean = sums.lambda.mean,
        amplification_mean = sums.amplification.mean,
        amplification_median = sums.amplification.median,
        amplification_q95 = sums.amplification.q95,
        amplification_max = sums.amplification.max
    ))
end

println("\n--- Long-term Summary Table ---")
for row in summary_rows
    println(row)
end

try
    @eval using DataFrames
    global summary_df = DataFrame(summary_rows)
    println("\nConverted summary_rows to DataFrame. summary_df is available.")
catch
    @info "DataFrames.jl not available; using summary_rows only."
end

# ==============================================================================
# LONG-TERM PLOTS IF REQUESTED
# ==============================================================================

if do_plots && !isempty(summary_rows)
    println("\n--- Generating Long-term Plots ---")

    x_vals = 1:length(summary_rows)

    fig_ts = Figure(size=(900, 420))
    ax_ts = Axis(fig_ts[1, 1], title="Amplification Factor Summary over Selected Months", xlabel="Time", ylabel="Amplification Factor F")

    lines!(ax_ts, x_vals, [r.amplification_median for r in summary_rows], label="Median")
    lines!(ax_ts, x_vals, [r.amplification_q95 for r in summary_rows], label="95th pct")
    lines!(ax_ts, x_vals, [r.amplification_max for r in summary_rows], label="Max")

    ax_ts.xticks = (x_vals, [r.key for r in summary_rows])
    axislegend(ax_ts)

    save(joinpath(output_dir, "amplification_timeseries_summary.png"), fig_ts, px_per_unit=2)

    first_key = summary_rows[1].key
    last_key = summary_rows[end].key

    plot_geo_grid(amplification_by_time[first_key], lon, lat,
        title="Amplification Factor ($first_key)",
        filename=joinpath(output_dir, "amplification_factor_$(replace(first_key, "-"=>"_")).png"))

    plot_geo_grid(amplification_by_time[last_key], lon, lat,
        title="Amplification Factor ($last_key)",
        filename=joinpath(output_dir, "amplification_factor_$(replace(last_key, "-"=>"_")).png"))

    if first_key != last_key
        Fdiff = amplification_by_time[last_key] .- amplification_by_time[first_key]
        fv_diff = finite_values(Fdiff)
        max_abs_diff = isempty(fv_diff) ? 1.0 : maximum(abs.(fv_diff))

        fig_diff = Figure(size=(800, 600))
        ax_diff = Axis(fig_diff[1, 1], title="Amplification Difference ($last_key minus $first_key)", xlabel="Longitude", ylabel="Latitude")
        diff_s, lon_s, lat_s = sort_for_plot(Fdiff, lon, lat)

        hm_diff = heatmap!(ax_diff, lon_s, lat_s, diff_s,
            colormap=:RdBu_11,
            colorrange=(-max_abs_diff, max_abs_diff),
            nan_color=:transparent)
        Colorbar(fig_diff[1, 2], hm_diff)

        save(joinpath(output_dir, "amplification_difference_$(replace(last_key, "-"=>"_"))_minus_$(replace(first_key, "-"=>"_")).png"), fig_diff, px_per_unit=2)
    end
end

# ==============================================================================
# FINAL CONCEPTUAL COMMENT
# ==============================================================================

println("\n------------------------------------------------------------")
println("The amplification factor F = lambda(A0, 0, 0) / lambda(A, s, Z) is a conditional vulnerability metric.")
println("It does not mean disease, collapse, or clinical risk occurs.")
println("It means that, under the reduced-DEB background physiological stress state,")
println("the same additional perturbation would produce a larger response burden.")
println("------------------------------------------------------------\n")
