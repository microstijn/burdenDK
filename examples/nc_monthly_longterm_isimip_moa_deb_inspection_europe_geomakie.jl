# ==============================================================================
# nc_monthly_longterm_isimip_moa_deb_inspection.jl
#
# Manual inspection script for real DynQual / ISIMIP-style monthly water-quality
# NetCDF files through the ISIMIP Water Quality -> MoA -> DEB grid pipeline.
#
# This version is configured for the local files supplied by Stijn:
#   waterTemperature_monthlyAvg_1980_2019.nc
#   organic_monthlyAvg_1980_2019.nc
#   pathogen_monthlyAvg_1980_2019.nc
#   TDSload_monthlyAvg_1980_2019.nc
#   BODload_monthlyAvg_1980_2019.nc
#
# Europe plotting:
#   - Europe bbox: lon -25..45, lat 34..72.
#   - GeoMakie GeoAxis in lon/lat mode to keep raster, coastlines, and borders aligned.
#   - NaturalEarth country borders if NaturalEarth.jl is installed.
#   - Grid/graticule disabled; map border drawn around bbox.
#
# Important data note:
#   - waterTemperature is used as WT.
#   - organic is used as BOD/organic concentration proxy.
#   - pathogen is used as FC/pathogen concentration proxy.
#   - TDSload is used as a TDS/salinity pressure proxy because no salinity/TDS
#     concentration file was supplied. This is a load, not a concentration.
#   - BODload is loaded only as an optional diagnostic/export layer; the pipeline
#     uses organic as BOD because organic is the concentration-like stress proxy.
#   - Nutrient, Chemical, and Plastic are set to zero for this first real-data run.
#
# The script uses shared robust normalisation across all selected times for the
# selected bbox. This preserves temporal differences better than normalising each
# month independently.
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
using GeoMakie

# NaturalEarth is used only for country borders. Keep this optional so the script
# can still run if the package is not installed.
const HAS_NATURALEARTH = let ok = false
    try
        @eval using NaturalEarth
        ok = true
    catch err
        @warn "NaturalEarth.jl not available; country borders will be skipped. Install with: import Pkg; Pkg.add(\"NaturalEarth\")" exception=(err, catch_backtrace())
        ok = false
    end
    ok
end

# ==============================================================================
# USER SETTINGS
# ==============================================================================

# Local DynQual files supplied by the user.
WT_file       = raw"C:\Users\peete074\Downloads\waterTemperature_monthlyAvg_1980_2019.nc"
BOD_file      = raw"C:\Users\peete074\Downloads\organic_monthlyAvg_1980_2019.nc"
FC_file       = raw"C:\Users\peete074\Downloads\pathogen_monthlyAvg_1980_2019.nc"
TDS_file      = raw"C:\Users\peete074\Downloads\TDSload_monthlyAvg_1980_2019.nc"
BODload_file  = raw"C:\Users\peete074\Downloads\BODload_monthlyAvg_1980_2019.nc"

# No real files yet for these layers in the current DynQual subset.
Nutrient_file = ""
Chemical_file = ""
Plastic_file  = ""

# Variable-name candidates. The script will try these first, then auto-detect the
# main gridded variable if none of these names are found.
WT_var_candidates       = ["waterTemperature", "watertemperature", "water_temperature", "waterTemp", "watertemp", "WT", "temperature"]
BOD_var_candidates      = ["organic", "BOD", "bod", "bod_concentration", "biological_oxygen_demand"]
FC_var_candidates       = ["pathogen", "FC", "fc", "fecal_coliform", "faecal_coliform"]
TDS_var_candidates      = ["TDSload", "tdsload", "TDSLoad", "tds_load", "TDS", "tds"]
BODload_var_candidates  = ["BODload", "bodload", "BODLoad", "bod_load"]

# If TDS concentration/salinity is not available, use TDSload as a stress proxy.
# This is useful for exploratory vulnerability mapping, but it is not the same as
# salinity concentration because it depends on pollutant transport/load.
use_TDSload_as_TDS_proxy = true

# Use organic concentration as BOD-like stressor. BODload is diagnostic only by default.
use_BODload_instead_of_organic = false

# If these missing canonical variables are not supplied, fill with zero rasters.
fill_missing_variables_with_zero = true

# Time settings.
start_year = 1980
end_year   = 2019
base_year  = 1980
selected_years = [1980, 1990, 2000, 2010, 2019]
selected_month = 7
process_all_months = false

# Europe-wide spatial subset.
bbox = (lon = (-25.0, 45.0), lat = (34.0, 72.0))
# bbox = nothing

selected_profile = fish_profile()

# Use the package default fish profile by default for real data.
# If you want a more visually demonstrative sensitivity run, set this to true and
# tune the DEBAxisParams below.
use_demo_sensitive_profile = false

# Plot and raster export settings.
do_plots = get(ENV, "TTR_EXAMPLE_PLOTS", "false") == "true"
export_beamer_rasters = true
export_ascii_rasters = true
export_png_rasters = true

# GeoMakie map settings. Use lon/lat to keep heatmap, coastlines, and country borders aligned.
use_geomakie_maps = true
show_coastlines = true
show_country_borders = true
show_map_border = true
show_geo_grid = false

coastline_color = :black
coastline_linewidth = 0.7
country_border_color = (:black, 0.65)
country_border_linewidth = 0.45
map_border_color = :black
map_border_linewidth = 1.2

geo_source = "+proj=longlat +datum=WGS84"
geo_dest   = "+proj=longlat +datum=WGS84"

output_dir = joinpath(@__DIR__, "..", "output", "nc_monthly_longterm_isimip_moa_deb_inspection")
beamer_output_dir = joinpath(output_dir, "beamer_rasters_real_dynqual")
mkpath(output_dir)
mkpath(beamer_output_dir)

# Normalisation settings.
normalisation_lower_q = 0.02
normalisation_upper_q = 0.98

# ==============================================================================
# OPTIONAL DEMO-SENSITIVE PROFILE
# ==============================================================================

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
        description = profile.description * " Demo-sensitive profile for real-data inspection."
    )
end

selected_profile = use_demo_sensitive_profile ? make_demo_sensitive_profile(selected_profile) : selected_profile

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

finite_values(x) = filter(isfinite, vec(collect(skipmissing(x))))

time_index_from_year_month(year, month; base_year=1980) = (year - base_year) * 12 + month

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
        error("bbox must be nothing, a 4-tuple, or a NamedTuple with lon and lat ranges")
    end
end

function replace_missing_with_nan(A)
    return Float64.(coalesce.(A, NaN))
end

function first_existing_key(ds, candidates)
    ks = collect(keys(ds))
    lower_map = Dict(lowercase(String(k)) => String(k) for k in ks)

    for cand in candidates
        if cand in ks
            return cand
        end
        lc = lowercase(cand)
        if haskey(lower_map, lc)
            return lower_map[lc]
        end
    end

    return nothing
end

function variable_dimnames(v)
    try
        return String.(dimnames(v))
    catch
        try
            return String.(NCDatasets.dimnames(v))
        catch
            return String[]
        end
    end
end

function infer_main_variable(path, candidates)
    ds = NCDataset(path, "r")
    try
        named = first_existing_key(ds, candidates)
        if named !== nothing
            return named
        end

        coord_like = Set(["lon", "longitude", "lat", "latitude", "time", "time_bnds", "bounds", "crs"])
        possible = String[]
        for k in keys(ds)
            name = String(k)
            lname = lowercase(name)
            if lname in coord_like
                continue
            end
            v = ds[name]
            nd = try
                ndims(v)
            catch
                0
            end
            if nd >= 2
                push!(possible, name)
            end
        end

        if isempty(possible)
            error("Could not infer a gridded data variable in $path. Available variables: $(collect(keys(ds)))")
        elseif length(possible) > 1
            @warn "Multiple possible data variables found; using first candidate" path=path candidates=possible chosen=possible[1]
        end
        return possible[1]
    finally
        close(ds)
    end
end

function lon_lat_names(ds)
    lon_name = first_existing_key(ds, ["lon", "longitude", "x"])
    lat_name = first_existing_key(ds, ["lat", "latitude", "y"])
    lon_name === nothing && error("No lon/longitude/x variable found. Available variables: $(collect(keys(ds)))")
    lat_name === nothing && error("No lat/latitude/y variable found. Available variables: $(collect(keys(ds)))")
    return lon_name, lat_name
end

function selected_lon_lat_indices(lon, lat, bbox)
    br = bbox_ranges(bbox)
    if br === nothing
        return collect(1:length(lon)), collect(1:length(lat))
    end

    lon_min, lon_max, lat_min, lat_max = br.lon_min, br.lon_max, br.lat_min, br.lat_max

    if minimum(skipmissing(lon)) >= 0 && lon_min < 0
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

    isempty(lon_idx) && error("bbox selected zero longitude cells. File lon range: $(minimum(lon)) to $(maximum(lon)); bbox lon: $((lon_min, lon_max))")
    isempty(lat_idx) && error("bbox selected zero latitude cells. File lat range: $(minimum(lat)) to $(maximum(lat)); bbox lat: $((lat_min, lat_max))")

    return lon_idx, lat_idx
end

function load_nc_layer_auto(path, candidates, time_index, bbox)
    if !isfile(path)
        error("File not found: $path")
    end

    varname = infer_main_variable(path, candidates)
    ds = NCDataset(path, "r")
    try
        lon_name, lat_name = lon_lat_names(ds)
        lon = ds[lon_name][:]
        lat = ds[lat_name][:]
        lon_idx, lat_idx = selected_lon_lat_indices(lon, lat, bbox)

        v = ds[varname]
        dims = variable_dimnames(v)
        if isempty(dims)
            @warn "Could not obtain dimension names; assuming dimensions are lon, lat, time" variable=varname path=path
            data = v[lon_idx, lat_idx, time_index]
            return replace_missing_with_nan(data), lon[lon_idx], lat[lat_idx], varname
        end

        idxs = Any[]
        kept_dims = String[]
        for d in dims
            dl = lowercase(d)
            if occursin("lon", dl) || dl == "x"
                push!(idxs, lon_idx)
                push!(kept_dims, "lon")
            elseif occursin("lat", dl) || dl == "y"
                push!(idxs, lat_idx)
                push!(kept_dims, "lat")
            elseif occursin("time", dl)
                push!(idxs, time_index)
            else
                push!(idxs, Colon())
                push!(kept_dims, d)
            end
        end

        raw = Array(v[idxs...])
        data = replace_missing_with_nan(raw)

        singleton_dims = Tuple(findall(size(data) .== 1))
        if !isempty(singleton_dims) && ndims(data) > 2
            data = dropdims(data; dims=singleton_dims)
        end

        if ndims(data) != 2
            error("Loaded variable $varname from $path but result is not 2D after slicing. Size: $(size(data)), dims: $dims")
        end

        kept = [d for d in kept_dims if lowercase(d) in ["lon", "lat", "x", "y"] || occursin("lon", lowercase(d)) || occursin("lat", lowercase(d))]
        if length(kept) >= 2
            first_dim = lowercase(kept[1])
            second_dim = lowercase(kept[2])
            if (occursin("lat", first_dim) || first_dim == "y") && (occursin("lon", second_dim) || second_dim == "x")
                data = permutedims(data)
            end
        elseif size(data, 1) == length(lat_idx) && size(data, 2) == length(lon_idx)
            data = permutedims(data)
        end

        return data, lon[lon_idx], lat[lat_idx], varname
    finally
        close(ds)
    end
end

function robust_reference_values(arrays; lower_q=0.02, upper_q=0.98)
    vals = Float64[]
    for A in arrays
        append!(vals, finite_values(A))
    end
    if isempty(vals)
        return (lo=NaN, hi=NaN)
    end
    lo = quantile(vals, lower_q)
    hi = quantile(vals, upper_q)
    return (lo=lo, hi=hi)
end

function apply_reference_normalisation(A, ref)
    out = similar(A, Float64)
    lo, hi = ref.lo, ref.hi
    if !isfinite(lo) || !isfinite(hi) || isapprox(lo, hi; atol=1e-12)
        for i in eachindex(A)
            out[i] = isfinite(A[i]) ? 0.0 : NaN
        end
        return out
    end
    for i in eachindex(A)
        out[i] = isfinite(A[i]) ? clamp((A[i] - lo) / (hi - lo), 0.0, 1.0) : NaN
    end
    return out
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

function sort_for_plot(layer, lon, lat)
    # Standardise longitude to -180..180 for plotting. This matters because many
    # global NetCDF files store longitude as 0..360, while Europe is easier to
    # plot as -25..45 around Greenwich.
    lon_plot = Float64.(lon)
    lon_plot = [x > 180 ? x - 360 : x for x in lon_plot]
    lat_plot = Float64.(lat)

    lon_order = sortperm(lon_plot)
    lat_order = sortperm(lat_plot)

    lon_sorted = lon_plot[lon_order]
    lat_sorted = lat_plot[lat_order]
    layer_sorted = copy(layer)[lon_order, lat_order]

    return layer_sorted, lon_sorted, lat_sorted
end

function plot_country_borders!(ax)
    if show_country_borders && HAS_NATURALEARTH
        try
            countries = naturalearth("admin_0_countries", 110)
            p = poly!(
                ax,
                GeoMakie.to_multipoly(countries.geometry);
                color = (:white, 0.0),
                strokecolor = country_border_color,
                strokewidth = country_border_linewidth,
                overdraw = true
            )
            try
                translate!(p, 0, 0, 100)
            catch
                # translate! is useful for projected/surface plots but not always needed.
            end
        catch err
            @warn "Could not plot NaturalEarth country borders" exception=(err, catch_backtrace())
        end
    end
    return nothing
end

function plot_map_border!(ax, br)
    if show_map_border && br !== nothing
        border_lon = [br.lon_min, br.lon_max, br.lon_max, br.lon_min, br.lon_min]
        border_lat = [br.lat_min, br.lat_min, br.lat_max, br.lat_max, br.lat_min]
        lines!(ax, border_lon, border_lat; color=map_border_color, linewidth=map_border_linewidth, overdraw=true)
    end
    return nothing
end

function plot_geo_grid(layer, lon, lat; title, filename, colormap=:viridis, colorrange=nothing)
    layer_s, lon_s, lat_s = sort_for_plot(layer, lon, lat)
    fig = Figure(size=(900, 650))

    if use_geomakie_maps
        br = bbox_ranges(bbox)
        map_limits = br === nothing ? nothing : (br.lon_min, br.lon_max, br.lat_min, br.lat_max)

        ax = if map_limits === nothing
            GeoAxis(
                fig[1, 1];
                source = geo_source,
                dest = geo_dest,
                title = title,
                xlabel = "Longitude",
                ylabel = "Latitude",
                #backgroundcolor = :white,
                xgridvisible = show_geo_grid,
                ygridvisible = show_geo_grid
            )
        else
            GeoAxis(
                fig[1, 1];
                source = geo_source,
                dest = geo_dest,
                limits = map_limits,
                title = title,
                xlabel = "Longitude",
                ylabel = "Latitude",
                #backgroundcolor = :white,
                xgridvisible = show_geo_grid,
                ygridvisible = show_geo_grid
            )
        end

        hm = if colorrange === nothing
            heatmap!(ax, lon_s, lat_s, layer_s; colormap=colormap, nan_color=:transparent)
        else
            heatmap!(ax, lon_s, lat_s, layer_s; colormap=colormap, colorrange=colorrange, nan_color=:transparent)
        end

        plot_country_borders!(ax)

        if show_coastlines
            lines!(ax, GeoMakie.coastlines(); color=coastline_color, linewidth=coastline_linewidth, overdraw=true)
        end

        plot_map_border!(ax, br)

        # Re-apply limits after plotting global vector layers so they do not expand the visible domain.
        if map_limits !== nothing
            xlims!(ax, br.lon_min, br.lon_max)
            ylims!(ax, br.lat_min, br.lat_max)
        end

        Colorbar(fig[1, 2], hm)
    else
        ax = Axis(
            fig[1, 1];
            title=title,
            xlabel="Longitude",
            ylabel="Latitude",
            xgridvisible=show_geo_grid,
            ygridvisible=show_geo_grid
        )
        hm = if colorrange === nothing
            heatmap!(ax, lon_s, lat_s, layer_s, colormap=colormap, nan_color=:transparent)
        else
            heatmap!(ax, lon_s, lat_s, layer_s, colormap=colormap, colorrange=colorrange, nan_color=:transparent)
        end
        Colorbar(fig[1, 2], hm)
    end

    save(filename, fig, px_per_unit=2)
    return fig
end

function write_ascii_grid(filename, grid)
    open(filename, "w") do io
        nrows, ncols = size(grid)
        println(io, "nrows $nrows")
        println(io, "ncols $ncols")
        println(io, "nodata_value NaN")
        for j in 1:ncols
            row = [grid[i, j] for i in 1:nrows]
            println(io, join(row, " "))
        end
    end
end

function safe_name(s::AbstractString)
    return replace(s, " "=>"_", "/"=>"_", "-"=>"_", ":"=>"_", "."=>"p")
end

function export_grid(grid, lon, lat; key, group, name, title, colormap=:viridis, colorrange=nothing)
    if !export_beamer_rasters
        return nothing
    end
    group_dir = joinpath(beamer_output_dir, safe_name(group))
    mkpath(group_dir)
    base = "$(safe_name(key))__$(safe_name(group))__$(safe_name(name))"
    if export_png_rasters
        plot_geo_grid(grid, lon, lat; title=title, filename=joinpath(group_dir, base * ".png"), colormap=colormap, colorrange=colorrange)
    end
    if export_ascii_rasters
        write_ascii_grid(joinpath(group_dir, base * ".asc"), grid)
    end
end

function load_raw_layers_for_time(year, month)
    t_idx = time_index_from_year_month(year, month; base_year=base_year)
    println("Loading real DynQual layers for $year-$(lpad(month, 2, "0")) (time index $t_idx)")

    WT, lon, lat, WT_var = load_nc_layer_auto(WT_file, WT_var_candidates, t_idx, bbox)
    BOD, _, _, BOD_var = if use_BODload_instead_of_organic
        load_nc_layer_auto(BODload_file, BODload_var_candidates, t_idx, bbox)
    else
        load_nc_layer_auto(BOD_file, BOD_var_candidates, t_idx, bbox)
    end
    FC, _, _, FC_var = load_nc_layer_auto(FC_file, FC_var_candidates, t_idx, bbox)
    TDS, _, _, TDS_var = load_nc_layer_auto(TDS_file, TDS_var_candidates, t_idx, bbox)

    BODload = nothing
    BODload_var = nothing
    if isfile(BODload_file)
        try
            BODload, _, _, BODload_var = load_nc_layer_auto(BODload_file, BODload_var_candidates, t_idx, bbox)
        catch err
            @warn "Could not load BODload diagnostic layer" exception=(err, catch_backtrace())
        end
    end

    ref_shape = size(WT)
    for (nm, layer) in [("BOD", BOD), ("FC", FC), ("TDS", TDS)]
        if size(layer) != ref_shape
            error("Layer $nm has size $(size(layer)), expected $ref_shape")
        end
    end

    if !use_TDSload_as_TDS_proxy
        error("TDS_file currently points to TDSload. Set use_TDSload_as_TDS_proxy=true or provide a salinity/TDS concentration file.")
    end

    zero_layer = fill(0.0, ref_shape)

    if fill_missing_variables_with_zero
        Nutrient = zero_layer
        Chemical = zero_layer
        Plastic = zero_layer
    else
        error("Nutrient/Chemical/Plastic files are missing and fill_missing_variables_with_zero=false")
    end

    metadata = (
        WT_var = WT_var,
        BOD_var = BOD_var,
        FC_var = FC_var,
        TDS_var = TDS_var,
        BODload_var = BODload_var,
        BOD_source = use_BODload_instead_of_organic ? "BODload" : "organic",
        TDS_source = "TDSload_proxy"
    )

    return (
        WT = WT,
        BOD = BOD,
        TDS = TDS,
        FC = FC,
        Nutrient = Nutrient,
        Chemical = Chemical,
        Plastic = Plastic,
        BODload = BODload,
        lon = lon,
        lat = lat,
        metadata = metadata
    )
end

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

println("\n--- Real DynQual ISIMIP MoA DEB NetCDF Inspection Script ---")
println("Profile: $(selected_profile.name)")
println("bbox: $bbox")
println("TDS source: TDSload as proxy = $use_TDSload_as_TDS_proxy")
println("BOD source: ", use_BODload_instead_of_organic ? "BODload" : "organic")
println("Nutrient/Chemical/Plastic: filled with zeros = $fill_missing_variables_with_zero")
println("GeoMakie maps: $use_geomakie_maps, country borders: $show_country_borders, NaturalEarth available: $HAS_NATURALEARTH")

selected_times = process_all_months ? [(y, m) for y in start_year:end_year for m in 1:12] : [(y, selected_month) for y in selected_years]

# Step 1: load all selected raw layers first.
raw_by_time = Dict{String, Any}()
lon = nothing
lat = nothing
for (yr, mo) in selected_times
    key = string(yr, "-", lpad(mo, 2, "0"))
    loaded = load_raw_layers_for_time(yr, mo)
    raw_by_time[key] = loaded
    if lon === nothing
        lon = loaded.lon
        lat = loaded.lat
    end
end

first_key = first(sort(collect(keys(raw_by_time))))
println("\nDetected variable names from first load:")
println(raw_by_time[first_key].metadata)
println("Raster size: ", size(raw_by_time[first_key].WT))
println("Longitude range: ", minimum(lon), " to ", maximum(lon))
println("Latitude range: ", minimum(lat), " to ", maximum(lat))

# Step 2: shared robust normalisation across selected times.
println("\nComputing shared robust normalisation references across selected times...")
refs = (
    WT = robust_reference_values([raw_by_time[k].WT for k in keys(raw_by_time)]; lower_q=normalisation_lower_q, upper_q=normalisation_upper_q),
    BOD = robust_reference_values([raw_by_time[k].BOD for k in keys(raw_by_time)]; lower_q=normalisation_lower_q, upper_q=normalisation_upper_q),
    TDS = robust_reference_values([raw_by_time[k].TDS for k in keys(raw_by_time)]; lower_q=normalisation_lower_q, upper_q=normalisation_upper_q),
    FC = robust_reference_values([raw_by_time[k].FC for k in keys(raw_by_time)]; lower_q=normalisation_lower_q, upper_q=normalisation_upper_q),
    Nutrient = (lo=0.0, hi=1.0),
    Chemical = (lo=0.0, hi=1.0),
    Plastic = (lo=0.0, hi=1.0)
)
println("Normalisation references:")
println(refs)

norm_by_time = Dict{String, Any}()
for key in keys(raw_by_time)
    raw = raw_by_time[key]
    norm_by_time[key] = (
        WT = apply_reference_normalisation(raw.WT, refs.WT),
        BOD = apply_reference_normalisation(raw.BOD, refs.BOD),
        TDS = apply_reference_normalisation(raw.TDS, refs.TDS),
        FC = apply_reference_normalisation(raw.FC, refs.FC),
        Nutrient = raw.Nutrient,
        Chemical = raw.Chemical,
        Plastic = raw.Plastic
    )
end

# Step 3: run pipeline for each selected time.
effective_layers_by_time = Dict{String, Any}()
modes_by_time = Dict{String, Any}()
axes_by_time = Dict{String, Any}()
Z_by_time = Dict{String, Any}()
A_by_time = Dict{String, Any}()
lambda_by_time = Dict{String, Any}()
amplification_by_time = Dict{String, Any}()
summaries_by_time = Dict{String, Any}()
summary_rows = NamedTuple[]

for key in sort(collect(keys(norm_by_time)))
    println("\nRunning pipeline for $key")
    nrm = norm_by_time[key]
    layers_norm = [nrm.WT, nrm.BOD, nrm.TDS, nrm.FC, nrm.Nutrient, nrm.Chemical, nrm.Plastic]

    res = isimip_deb_pipeline_grid(
        layers_norm,
        selected_profile.exposure_filter,
        selected_profile.moa_mapping,
        selected_profile.moa_deb_mapping,
        selected_profile.deb_params
    )

    effective_layers_by_time[key] = res.effective_layers
    modes_by_time[key] = res.modes
    axes_by_time[key] = res.axes
    Z_by_time[key] = res.Z
    A_by_time[key] = res.A
    lambda_by_time[key] = res.lambda
    amplification_by_time[key] = res.amplification

    sums = (
        WT = monthly_summary(nrm.WT),
        BOD = monthly_summary(nrm.BOD),
        TDS = monthly_summary(nrm.TDS),
        FC = monthly_summary(nrm.FC),
        thermal = monthly_summary(res.modes.thermal),
        oxygen = monthly_summary(res.modes.oxygen),
        osmotic = monthly_summary(res.modes.osmotic),
        immune = monthly_summary(res.modes.immune),
        eutrophication = monthly_summary(res.modes.eutrophication),
        toxic = monthly_summary(res.modes.toxic),
        feeding = monthly_summary(res.modes.feeding),
        physical = monthly_summary(res.modes.physical),
        assimilation = monthly_summary(res.axes.assimilation),
        maintenance = monthly_summary(res.axes.maintenance),
        growth = monthly_summary(res.axes.growth),
        reproduction = monthly_summary(res.axes.reproduction),
        A = monthly_summary(res.A),
        lambda = monthly_summary(res.lambda),
        amplification = monthly_summary(res.amplification)
    )
    summaries_by_time[key] = sums

    yr = parse(Int, split(key, "-")[1])
    mo = parse(Int, split(key, "-")[2])
    push!(summary_rows, (
        year = yr,
        month = mo,
        key = key,
        WT_mean = sums.WT.mean,
        BOD_mean = sums.BOD.mean,
        TDS_mean = sums.TDS.mean,
        FC_mean = sums.FC.mean,
        thermal_mean = sums.thermal.mean,
        oxygen_mean = sums.oxygen.mean,
        osmotic_mean = sums.osmotic.mean,
        immune_mean = sums.immune.mean,
        maintenance_mean = sums.maintenance.mean,
        A_mean = sums.A.mean,
        lambda_mean = sums.lambda.mean,
        amplification_mean = sums.amplification.mean,
        amplification_median = sums.amplification.median,
        amplification_q95 = sums.amplification.q95,
        amplification_max = sums.amplification.max
    ))
end

println("\n--- Summary Table ---")
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

# Step 4: export rasters.
if export_beamer_rasters
    println("\nExporting rasters to: $beamer_output_dir")
    for key in sort(collect(keys(norm_by_time)))
        raw = raw_by_time[key]
        nrm = norm_by_time[key]

        # Raw / normalised input layers.
        for (nm, grid) in pairs((WT=raw.WT, BOD=raw.BOD, TDS=raw.TDS, FC=raw.FC))
            export_grid(grid, lon, lat; key=key, group="01_raw_inputs", name=String(nm), title="Raw $(String(nm)) ($key)", colormap=:viridis)
        end
        if raw.BODload !== nothing
            export_grid(raw.BODload, lon, lat; key=key, group="01_raw_inputs", name="BODload_diagnostic", title="Raw BODload diagnostic ($key)", colormap=:viridis)
        end
        for (nm, grid) in pairs((WT=nrm.WT, BOD=nrm.BOD, TDS=nrm.TDS, FC=nrm.FC, Nutrient=nrm.Nutrient, Chemical=nrm.Chemical, Plastic=nrm.Plastic))
            export_grid(grid, lon, lat; key=key, group="02_normalised_inputs", name=String(nm), title="Normalised $(String(nm)) ($key)", colormap=:viridis, colorrange=(0.0, 1.0))
        end

        # Modes.
        modes = modes_by_time[key]
        for nm in (:thermal, :oxygen, :osmotic, :immune, :eutrophication, :toxic, :feeding, :physical)
            export_grid(getproperty(modes, nm), lon, lat; key=key, group="04_modes", name=String(nm), title="MoA $(String(nm)) ($key)", colormap=:magma, colorrange=(0.0, 1.0))
        end

        # Axes.
        axes = axes_by_time[key]
        for nm in (:assimilation, :maintenance, :growth, :reproduction)
            export_grid(getproperty(axes, nm), lon, lat; key=key, group="05_deb_axes", name=String(nm), title="DEB axis $(String(nm)) ($key)", colormap=:plasma, colorrange=(0.0, 1.0))
        end

        # Response outputs.
        export_grid(A_by_time[key], lon, lat; key=key, group="06_response_operator", name="adaptive_margin_A", title="Adaptive margin A ($key)", colormap=:viridis)
        export_grid(lambda_by_time[key], lon, lat; key=key, group="06_response_operator", name="restoring_force_lambda", title="Restoring force lambda ($key)", colormap=:viridis)
        export_grid(amplification_by_time[key], lon, lat; key=key, group="06_response_operator", name="amplification_F", title="Amplification factor F ($key)", colormap=:inferno)
    end

    sorted_keys = sort(collect(keys(amplification_by_time)))
    if length(sorted_keys) >= 2
        k0 = first(sorted_keys)
        k1 = last(sorted_keys)
        Fdiff = amplification_by_time[k1] .- amplification_by_time[k0]
        fd = finite_values(Fdiff)
        max_abs = isempty(fd) ? 1.0 : maximum(abs.(fd))
        export_grid(Fdiff, lon, lat; key="$(k0)_to_$(k1)", group="07_differences", name="amplification_F_difference", title="Amplification F difference: $k1 minus $k0", colormap=:RdBu_11, colorrange=(-max_abs, max_abs))

        Adiff = A_by_time[k1] .- A_by_time[k0]
        ad = finite_values(Adiff)
        max_abs_A = isempty(ad) ? 1.0 : maximum(abs.(ad))
        export_grid(Adiff, lon, lat; key="$(k0)_to_$(k1)", group="07_differences", name="adaptive_margin_A_difference", title="Adaptive margin A difference: $k1 minus $k0", colormap=:RdBu_11, colorrange=(-max_abs_A, max_abs_A))
    end
end

println("\n------------------------------------------------------------")
println("Real-data run complete.")
println("Note: TDS uses TDSload as a pressure proxy, not concentration, because no salinity/TDS concentration file was supplied.")
println("Note: Nutrient, Chemical, and Plastic are zero-filled in this run.")
println("The amplification factor F remains a conditional vulnerability metric, not a disease prediction.")
println("------------------------------------------------------------\n")
