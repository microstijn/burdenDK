# ==============================================================================
# nc_monthly_longterm_deb_axes_inspection.jl
#
# Example script to inspect how long-term monthly background exposure to organic
# and pathogen stressors produces a larger amplification factor over time, mapped
# to DEB-like physiological axes.
#
# This script is written as a linear, inspectable workflow. All intermediate
# grids and time-series dictionaries are available in the top-level scope for
# REPL/VS Code inspection.
# ==============================================================================

using TwoTimescaleResilience
using NCDatasets
using Statistics
using CairoMakie

# ==============================================================================
# User settings section
# ==============================================================================

organic_file  = raw"C:\Users\peete074\Downloads\organic_monthlyAvg_1980_2019.nc"
pathogen_file = raw"C:\Users\peete074\Downloads\pathogen_monthlyAvg_1980_2019.nc"

organic_var = "organic"
pathogen_var = "pathogen"

start_year = 1980
end_year   = 2019

# inspect only selected years/months first
selected_years = [1980, 1990, 2000, 2010, 2019]
selected_month = 7

# optional crop
# Europe:
bbox = (-12.0, 35.0, 34.0, 72.0)

# For global:
# bbox = nothing

# To process all months (default false to keep memory/time manageable):
process_all_months = false

output_dir = joinpath(@__DIR__, "..", "output", "nc_monthly_longterm_deb_axes_inspection")
mkpath(output_dir)

# ==============================================================================
# Helper functions
# ==============================================================================

# 1. finite_values
function finite_values(x)
    return filter(isfinite, x)
end

# 2. robust_normalise
function robust_normalise(x; lower_q=0.02, upper_q=0.98, log_transform=false)
    fv = finite_values(x)
    if isempty(fv)
        return copy(x)
    end

    if log_transform
        fv = log10.(1.0 .+ max.(fv, 0.0))
    end

    lo = quantile(fv, lower_q)
    hi = quantile(fv, upper_q)

    if (hi - lo == 0) || isapprox(hi, lo; atol=1e-10)
        @warn "Normalisation range is zero. Returning zeros for finite cells."
        out = similar(x, Float64)
        for i in eachindex(x)
            if isfinite(x[i])
                out[i] = 0.0
            else
                out[i] = NaN
            end
        end
        return out
    end

    out = similar(x, Float64)
    for i in eachindex(x)
        if !isfinite(x[i])
            out[i] = NaN
        else
            val = x[i]
            if log_transform
                val = log10(1.0 + max(val, 0.0))
            end
            norm_val = (val - lo) / (hi - lo)
            out[i] = clamp(norm_val, 0.0, 1.0)
        end
    end
    return out
end

# 3. replace_missing_with_nan
function replace_missing_with_nan(A)
    return Matrix{Float64}(coalesce.(A, NaN))
end

# 4. load_nc_slice
function load_nc_slice(path, varname; time_index, bbox=nothing)
    ds = NCDataset(path, "r")
    lon = ds["longitude"][:]
    lat = ds["latitude"][:]

    println("Longitude range in file: ", minimum(lon), " to ", maximum(lon))

    if bbox !== nothing
        lon_min, lon_max, lat_min, lat_max = bbox
        if minimum(lon) >= 0 && lon_min < 0
            # file is 0..360 but bbox contains negative longitudes
            lon_min_mod = mod(lon_min, 360)
            lon_max_mod = mod(lon_max, 360)
            if lon_min_mod > lon_max_mod
                # crosses the 0-degree meridian
                lon_idx = findall(x -> (x >= lon_min_mod) || (x <= lon_max_mod), lon)
            else
                lon_idx = findall(x -> lon_min_mod <= x <= lon_max_mod, lon)
            end
        else
            # standard selection
            lon_idx = findall(x -> lon_min <= x <= lon_max, lon)
        end
        lat_idx = findall(x -> lat_min <= x <= lat_max, lat)
    else
        lon_idx = 1:length(lon)
        lat_idx = 1:length(lat)
    end

    data_slice = ds[varname][lon_idx, lat_idx, time_index]
    lon_selected = lon[lon_idx]
    lat_selected = lat[lat_idx]
    close(ds)

    layer = replace_missing_with_nan(data_slice)
    return layer, lon_selected, lat_selected
end

# 5. time_index_from_year_month
function time_index_from_year_month(year, month)
    return (year - 1980) * 12 + month
end

# 6. sort_for_plot
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

# 7. plot_geo_grid
function plot_geo_grid(layer, lon, lat; title, filename, colormap=:viridis)
    layer_s, lon_s, lat_s = sort_for_plot(layer, lon, lat)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], title=title, xlabel="Longitude", ylabel="Latitude")
    hm = heatmap!(ax, lon_s, lat_s, layer_s, colormap=colormap, nan_color=:transparent)
    Colorbar(fig[1, 2], hm)

    save(filename, fig, px_per_unit=2)
    return fig
end

# 8. monthly_summary
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

# ==============================================================================
# DEB mapping
# ==============================================================================

mapping = default_pathogen_organic_deb_mapping(
    interaction_strength = 0.25,
    clamp_unit = false
)

params = DEBAxisParams(
    A0 = 1.0,
    alpha_axes = (0.30, 0.35, 0.20, 0.15),
    lambda_min = 0.04,
    lambda_max = 1.0,
    KA = 0.30
)

# ==============================================================================
# Main workflow - Step 1: Load example month
# ==============================================================================

example_year = 2010
example_month = 7
example_time_index = time_index_from_year_month(example_year, example_month)

println("\n--- Step 1: Loading Example Month (Year: $example_year, Month: $example_month, Time Index: $example_time_index) ---")

pathogen_raw, lon, lat = load_nc_slice(pathogen_file, pathogen_var; time_index=example_time_index, bbox=bbox)
organic_raw, lon2, lat2 = load_nc_slice(organic_file, organic_var; time_index=example_time_index, bbox=bbox)

@assert lon == lon2 "Longitude coordinates do not match between datasets!"
@assert lat == lat2 "Latitude coordinates do not match between datasets!"

println("Raster sizes: ", size(pathogen_raw))
println("Longitude range: ", minimum(lon), " to ", maximum(lon))
println("Latitude range: ", minimum(lat), " to ", maximum(lat))

# ==============================================================================
# Step 2: Normalise
# ==============================================================================

println("\n--- Step 2: Normalising Example Month ---")
pathogen_norm = robust_normalise(pathogen_raw; lower_q=0.02, upper_q=0.98, log_transform=true)
organic_norm  = robust_normalise(organic_raw; lower_q=0.02, upper_q=0.98, log_transform=false)

# ==============================================================================
# Step 3: Run DEB axis pipeline
# ==============================================================================

println("\n--- Step 3: Running DEB Pipeline ---")
layers = [pathogen_norm, organic_norm]
result = deb_amplification_pipeline(layers, mapping, params)

axes = result.axes
assimilation_grid  = axes.assimilation
maintenance_grid   = axes.maintenance
growth_grid        = axes.growth
reproduction_grid  = axes.reproduction

Agrid      = result.A
lambdagrid = result.lambda
Fgrid      = result.amplification

if any(finite_values(Fgrid) .< 1.0)
    @warn "Amplification factor has values below 1! Stress should not increase restoring force under current assumptions."
end

# ==============================================================================
# Step 4: Print summaries
# ==============================================================================

println("\n--- Step 4: Summaries for Example Month ---")
println("pathogen_norm:     ", monthly_summary(pathogen_norm))
println("organic_norm:      ", monthly_summary(organic_norm))
println("assimilation_grid: ", monthly_summary(assimilation_grid))
println("maintenance_grid:  ", monthly_summary(maintenance_grid))
println("growth_grid:       ", monthly_summary(growth_grid))
println("reproduction_grid: ", monthly_summary(reproduction_grid))
println("Agrid:             ", monthly_summary(Agrid))
println("lambdagrid:        ", monthly_summary(lambdagrid))
println("Fgrid:             ", monthly_summary(Fgrid))

# ==============================================================================
# Step 5: Plot selected intermediate rasters
# ==============================================================================

println("\n--- Step 5: Plotting Example Month ---")
plot_geo_grid(pathogen_norm, lon, lat, title="Pathogen Norm ($example_year-$example_month)", filename=joinpath(output_dir, "pathogen_norm_$(example_year)_$(example_month).png"))
plot_geo_grid(organic_norm, lon, lat, title="Organic Norm ($example_year-$example_month)", filename=joinpath(output_dir, "organic_norm_$(example_year)_$(example_month).png"))
plot_geo_grid(assimilation_grid, lon, lat, title="Assimilation Cost ($example_year-$example_month)", filename=joinpath(output_dir, "assimilation_grid_$(example_year)_$(example_month).png"))
plot_geo_grid(maintenance_grid, lon, lat, title="Maintenance Cost ($example_year-$example_month)", filename=joinpath(output_dir, "maintenance_grid_$(example_year)_$(example_month).png"))
plot_geo_grid(growth_grid, lon, lat, title="Growth Cost ($example_year-$example_month)", filename=joinpath(output_dir, "growth_grid_$(example_year)_$(example_month).png"))
plot_geo_grid(reproduction_grid, lon, lat, title="Reproduction Cost ($example_year-$example_month)", filename=joinpath(output_dir, "reproduction_grid_$(example_year)_$(example_month).png"))
plot_geo_grid(Agrid, lon, lat, title="Adaptive Margin ($example_year-$example_month)", filename=joinpath(output_dir, "Agrid_$(example_year)_$(example_month).png"))
plot_geo_grid(lambdagrid, lon, lat, title="Restoring Force ($example_year-$example_month)", filename=joinpath(output_dir, "lambdagrid_$(example_year)_$(example_month).png"))
plot_geo_grid(Fgrid, lon, lat, title="Amplification Factor ($example_year-$example_month)", filename=joinpath(output_dir, "Fgrid_$(example_year)_$(example_month).png"))

# ==============================================================================
# Long-term monthly analysis
# ==============================================================================

println("\n--- Long-term monthly analysis ---")

if process_all_months
    selected_times = []
    for y in start_year:end_year
        for m in 1:12
            push!(selected_times, (y, m))
        end
    end
else
    selected_times = [(year, selected_month) for year in selected_years]
end

raw_layers_by_time = Dict()
norm_layers_by_time = Dict()
axes_by_time = Dict()
A_by_time = Dict()
lambda_by_time = Dict()
amplification_by_time = Dict()
summaries_by_time = Dict()

summary_rows = NamedTuple[]

for (y, m) in selected_times
    t_idx = time_index_from_year_month(y, m)
    key = string(y, "-", lpad(m, 2, "0"))

    println("Processing $key...")

    p_raw, plon, plat = load_nc_slice(pathogen_file, pathogen_var; time_index=t_idx, bbox=bbox)
    o_raw, olon, olat = load_nc_slice(organic_file, organic_var; time_index=t_idx, bbox=bbox)

    p_norm = robust_normalise(p_raw; lower_q=0.02, upper_q=0.98, log_transform=true)
    o_norm = robust_normalise(o_raw; lower_q=0.02, upper_q=0.98, log_transform=false)

    res = deb_amplification_pipeline([p_norm, o_norm], mapping, params)

    if any(finite_values(res.amplification) .< 1.0)
        @warn "Amplification factor for $key has values below 1!"
    end

    # Store in memory
    raw_layers_by_time[key] = (pathogen = p_raw, organic = o_raw)
    norm_layers_by_time[key] = (pathogen = p_norm, organic = o_norm)
    axes_by_time[key] = res.axes
    A_by_time[key] = res.A
    lambda_by_time[key] = res.lambda
    amplification_by_time[key] = res.amplification

    sums = (
        pathogen = monthly_summary(p_norm),
        organic = monthly_summary(o_norm),
        assimilation = monthly_summary(res.axes.assimilation),
        maintenance = monthly_summary(res.axes.maintenance),
        growth = monthly_summary(res.axes.growth),
        reproduction = monthly_summary(res.axes.reproduction),
        A = monthly_summary(res.A),
        lambda = monthly_summary(res.lambda),
        amplification = monthly_summary(res.amplification)
    )
    summaries_by_time[key] = sums

    row = (
        year = y,
        month = m,
        key = key,
        pathogen_mean = sums.pathogen.mean,
        organic_mean = sums.organic.mean,
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
    )
    push!(summary_rows, row)
end

println("\n--- Long-term Summary Table ---")
for r in summary_rows
    println(r)
end

try
    @eval using DataFrames
    global summary_df = DataFrame(summary_rows)
    println("\nConverted summary_rows to DataFrames (summary_df is available).")
catch
    @info "DataFrames.jl not available; using summary_rows only."
end

# ==============================================================================
# Long-term plots
# ==============================================================================

println("\n--- Generating Long-term Plots ---")
if length(selected_times) > 0
    # Timeseries plot of amplification summary
    fig_ts = Figure(size=(800, 400))
    ax_ts = Axis(fig_ts[1, 1], title="Amplification Factor Summary over Selected Months",
                 xlabel="Time (Index)", ylabel="Amplification Factor F")

    x_vals = 1:length(summary_rows)
    lines!(ax_ts, x_vals, [r.amplification_median for r in summary_rows], label="Median")
    lines!(ax_ts, x_vals, [r.amplification_q95 for r in summary_rows], label="95th Pctl")
    lines!(ax_ts, x_vals, [r.amplification_max for r in summary_rows], label="Max")

    ax_ts.xticks = (x_vals, [r.key for r in summary_rows])
    axislegend(ax_ts)

    save(joinpath(output_dir, "amplification_timeseries_summary.png"), fig_ts, px_per_unit=2)

    # Map for first and last month
    first_key = summary_rows[1].key
    last_key = summary_rows[end].key

    plot_geo_grid(amplification_by_time[first_key], lon, lat,
        title="Amplification Factor ($first_key)",
        filename=joinpath(output_dir, "amplification_factor_$(replace(first_key, "-"=>"_")).png"))

    plot_geo_grid(amplification_by_time[last_key], lon, lat,
        title="Amplification Factor ($last_key)",
        filename=joinpath(output_dir, "amplification_factor_$(replace(last_key, "-"=>"_")).png"))

    # Difference map (last - first)
    if first_key != last_key
        Fdiff = amplification_by_time[last_key] .- amplification_by_time[first_key]

        # Use a diverging colormap centered at 0 if possible
        # Finding the max absolute value for symmetric colorrange
        max_abs_diff = maximum(abs.(finite_values(Fdiff)))

        fig_diff = Figure(size=(800, 600))
        ax_diff = Axis(fig_diff[1, 1], title="Amplification Difference ($last_key minus $first_key)", xlabel="Longitude", ylabel="Latitude")
        diff_s, lon_s, lat_s = sort_for_plot(Fdiff, lon, lat)

        hm_diff = heatmap!(ax_diff, lon_s, lat_s, diff_s,
            colormap=:RdBu_11, colorrange=(-max_abs_diff, max_abs_diff), nan_color=:transparent)
        Colorbar(fig_diff[1, 2], hm_diff)

        save(joinpath(output_dir, "amplification_difference_$(replace(last_key, "-"=>"_"))_minus_$(replace(first_key, "-"=>"_")).png"), fig_diff, px_per_unit=2)
    end
end

# ==============================================================================
# Final conceptual comment
# ==============================================================================

println("\n------------------------------------------------------------")
println("The amplification factor F = lambda(A0) / lambda(A_DEB) is a conditional vulnerability metric.")
println("It does not mean that disease or collapse occurs.")
println("It means that, under the DEB-like background physiological stress state,")
println("the same additional perturbation would produce a larger response burden.")
println("------------------------------------------------------------\n")
