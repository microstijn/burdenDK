# This is a manual inspection script.
# It is intentionally not included in automated tests because it requires local NetCDF files.

using TwoTimescaleResilience
using NCDatasets
using CairoMakie
using Statistics

# ==============================================================================
# USER SETTINGS
# ==============================================================================
# Edit these paths and variable names to match your local ISIMIP Water Quality files.
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

fill_missing_variables_with_zero = true

bbox = (lon = (100.0, 160.0), lat = (-45.0, 5.0))
base_year = 1980

process_all_months = false
selected_years = [1980, 1990, 2000, 2010, 2019]
selected_month = 7

selected_profile = fish_profile()

do_plots = get(ENV, "TTR_EXAMPLE_PLOTS", "false") == "true"
do_plots = true # Force to true for manual inspection if needed

# Helper functions similar to existing script
finite_values(x) = filter(v -> isfinite(v) && !ismissing(v), x)

function robust_normalise(data)
    # Exclude missing and NaN values
    valid_data = finite_values(data)
    if isempty(valid_data)
        return fill(NaN, size(data))
    end

    p5 = quantile(valid_data, 0.05)
    p95 = quantile(valid_data, 0.95)

    norm_data = similar(data, Float64)
    if p95 == p5
        # Range is zero, return 0.0 for all finite cells
        for i in eachindex(data)
            norm_data[i] = isfinite(data[i]) ? 0.0 : NaN
        end
    else
        for i in eachindex(data)
            if isfinite(data[i])
                v = (data[i] - p5) / (p95 - p5)
                norm_data[i] = clamp(v, 0.0, 1.0)
            else
                norm_data[i] = NaN
            end
        end
    end
    return norm_data
end

replace_missing_with_nan(data) = map(x -> ismissing(x) ? NaN : Float64(x), data)

time_index_from_year_month(year, month; base_year=1980) = (year - base_year) * 12 + month

# ==============================================================================
# SCRIPT LOGIC
# ==============================================================================
println("--- ISIMIP MoA DEB NetCDF Inspection Script ---")
println("Profile: $(selected_profile.name)")

# This is a stub implementation to show the structure. Real loading requires NCDatasets logic.
# Because this script won't be run by automated tests, we simulate the loading process
# here so the user can easily replace the dummy loads with real load_nc_layer calls.

# For demo purposes, we will mock the "loaded" data if the files don't exist
function mock_load_nc_layer(file, var, time_idx, bbox)
    if isfile(file)
        # Real load if file exists
        return load_nc_layer(file, var, time_idx, bbox)
    else
        # Mock load for structure demonstration
        println("Warning: Mocking data for $var because file not found: $file")
        return fill(rand(), 50, 60), collect(range(bbox.lon[1], bbox.lon[2], length=60)), collect(range(bbox.lat[1], bbox.lat[2], length=50))
    end
end

function get_layer_or_zero(file, var, time_idx, bbox, ref_shape)
    if file != "" && isfile(file)
        data, lon, lat = mock_load_nc_layer(file, var, time_idx, bbox)
        return replace_missing_with_nan(data), lon, lat
    else
        if fill_missing_variables_with_zero
            println("Warning: Variable $var missing. Filling with zeros.")
            return fill(0.0, ref_shape...), nothing, nothing
        else
            error("Missing file for $var and fill_missing_variables_with_zero is false.")
        end
    end
end

# Load one example month first
t_idx_example = time_index_from_year_month(selected_years[1], selected_month; base_year=base_year)

WT_raw, lon, lat = mock_load_nc_layer(WT_file, WT_var, t_idx_example, bbox)
ref_shape = size(WT_raw)

BOD_raw, _, _ = get_layer_or_zero(BOD_file, BOD_var, t_idx_example, bbox, ref_shape)
TDS_raw, _, _ = get_layer_or_zero(TDS_file, TDS_var, t_idx_example, bbox, ref_shape)
FC_raw, _, _ = get_layer_or_zero(FC_file, FC_var, t_idx_example, bbox, ref_shape)
Nutrient_raw, _, _ = get_layer_or_zero(Nutrient_file, Nutrient_var, t_idx_example, bbox, ref_shape)
Chemical_raw, _, _ = get_layer_or_zero(Chemical_file, Chemical_var, t_idx_example, bbox, ref_shape)
Plastic_raw, _, _ = get_layer_or_zero(Plastic_file, Plastic_var, t_idx_example, bbox, ref_shape)

WT_norm = robust_normalise(replace_missing_with_nan(WT_raw))
BOD_norm = robust_normalise(BOD_raw)
TDS_norm = robust_normalise(TDS_raw)
FC_norm = robust_normalise(FC_raw)
Nutrient_norm = robust_normalise(Nutrient_raw)
Chemical_norm = robust_normalise(Chemical_raw)
Plastic_norm = robust_normalise(Plastic_raw)

layers_norm = [WT_norm, BOD_norm, TDS_norm, FC_norm, Nutrient_norm, Chemical_norm, Plastic_norm]

res = isimip_deb_pipeline_grid(layers_norm, selected_profile.exposure_filter, selected_profile.moa_mapping, selected_profile.moa_deb_mapping, selected_profile.deb_params)

effective_layers = res.effective_layers
modes = res.modes
axes = res.axes
Zgrid = res.Z
Agrid = res.A
lambdagrid = res.lambda
Fgrid = res.amplification

println("Successfully ran pipeline on example month.")

# Loop through selected months
raw_layers_by_time = []
norm_layers_by_time = []
effective_layers_by_time = []
modes_by_time = []
axes_by_time = []
Z_by_time = []
A_by_time = []
lambda_by_time = []
amplification_by_time = []
summaries_by_time = []
summary_rows = NamedTuple[]

for yr in selected_years
    t_idx = time_index_from_year_month(yr, selected_month; base_year=base_year)
    println("Processing $yr-$selected_month (index $t_idx)")

    wt, _, _ = mock_load_nc_layer(WT_file, WT_var, t_idx, bbox)
    bod, _, _ = get_layer_or_zero(BOD_file, BOD_var, t_idx, bbox, ref_shape)
    tds, _, _ = get_layer_or_zero(TDS_file, TDS_var, t_idx, bbox, ref_shape)
    fc, _, _ = get_layer_or_zero(FC_file, FC_var, t_idx, bbox, ref_shape)
    nutr, _, _ = get_layer_or_zero(Nutrient_file, Nutrient_var, t_idx, bbox, ref_shape)
    chem, _, _ = get_layer_or_zero(Chemical_file, Chemical_var, t_idx, bbox, ref_shape)
    plast, _, _ = get_layer_or_zero(Plastic_file, Plastic_var, t_idx, bbox, ref_shape)

    lyrs_norm = [
        robust_normalise(replace_missing_with_nan(wt)),
        robust_normalise(bod),
        robust_normalise(tds),
        robust_normalise(fc),
        robust_normalise(nutr),
        robust_normalise(chem),
        robust_normalise(plast)
    ]

    res_t = isimip_deb_pipeline_grid(lyrs_norm, selected_profile.exposure_filter, selected_profile.moa_mapping, selected_profile.moa_deb_mapping, selected_profile.deb_params)

    push!(raw_layers_by_time, [wt, bod, tds, fc, nutr, chem, plast])
    push!(norm_layers_by_time, lyrs_norm)
    push!(effective_layers_by_time, res_t.effective_layers)
    push!(modes_by_time, res_t.modes)
    push!(axes_by_time, res_t.axes)
    push!(Z_by_time, res_t.Z)
    push!(A_by_time, res_t.A)
    push!(lambda_by_time, res_t.lambda)
    push!(amplification_by_time, res_t.amplification)

    f_vals = finite_values(res_t.amplification)
    if isempty(f_vals)
        push!(summary_rows, (year=yr, mean_F=NaN, max_F=NaN))
    else
        push!(summary_rows, (year=yr, mean_F=mean(f_vals), max_F=maximum(f_vals)))
    end
end

println("Summary statistics across time:")
for row in summary_rows
    println(row)
end

if do_plots
    println("Plotting... (To disable, set ENV[\"TTR_EXAMPLE_PLOTS\"] = \"false\")")
    # Add plot calls here if manual inspection is desired.
    # e.g., heatmap!(ax, lon, lat, Fgrid)
end
