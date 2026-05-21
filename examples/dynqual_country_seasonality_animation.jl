# ==============================================================================
# dynqual_country_seasonality_animation.jl
#
# Secondary real-data script for DynQual / ISIMIP-style monthly NetCDF files.
#
# Goal:
#   1. Compute a monthly Europe-wide K grid for every month from 1980-2019.
#   2. Aggregate K by selected countries and plot seasonality / monthly time series.
#   3. Generate a Europe animation of K through time using Makie / GeoMakie.
#
# What is K here?
#   By default, K is the amplification factor F from the MoA -> DEB pipeline:
#       K = F = lambda_control / lambda_background
#   You can change K_metric below to :lambda, :A, :maintenance, :oxygen, etc.
#
# No raster grids are saved. Only plots and animation are written.
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
using GeoInterface
using CSV

# NaturalEarth is used for country polygons / country borders.
# Install once if needed:
#   import Pkg; Pkg.add("NaturalEarth")
using NaturalEarth
using Tables

# ==============================================================================
# USER SETTINGS
# ==============================================================================

# Local DynQual files.
WT_file       = raw"C:\Users\peete074\Downloads\waterTemperature_monthlyAvg_1980_2019.nc"
BOD_file      = raw"C:\Users\peete074\Downloads\organic_monthlyAvg_1980_2019.nc"
FC_file       = raw"C:\Users\peete074\Downloads\pathogen_monthlyAvg_1980_2019.nc"
TDS_file      = raw"C:\Users\peete074\Downloads\TDSload_monthlyAvg_1980_2019.nc"
BODload_file  = raw"C:\Users\peete074\Downloads\BODload_monthlyAvg_1980_2019.nc"

WT_var_candidates  = ["waterTemperature", "watertemperature", "water_temperature", "waterTemp", "watertemp", "WT", "temperature"]
BOD_var_candidates = ["organic", "BOD", "bod", "bod_concentration", "biological_oxygen_demand"]
FC_var_candidates  = ["pathogen", "FC", "fc", "fecal_coliform", "faecal_coliform"]
TDS_var_candidates = ["TDSload", "tdsload", "TDSLoad", "tds_load", "TDS", "tds"]

# Europe bbox.
bbox = (lon = (-25.0, 45.0), lat = (34.0, 72.0))

# Time range.
start_year = 1980
end_year = 2019
base_year = 1980
months = 1:12
all_times = [(y, m) for y in start_year:end_year for m in months]

# Countries to summarize. Names are matched fuzzily against NaturalEarth fields.
selected_countries = [
    "Netherlands",
    "Germany",
    "France",
    "Spain",
    "Italy",
    "Poland",
    "United Kingdom",
    "Sweden"
]

# K metric to aggregate and animate.
# Options implemented in get_K_grid(...):
#   :amplification, :lambda, :A,
#   :thermal, :oxygen, :osmotic, :immune, :eutrophication, :toxic, :feeding, :physical,
#   :assimilation, :maintenance, :growth, :reproduction
K_metric = :amplification
K_label = "Amplification factor F"

# Country statistic.
country_statistic = :mean  # :mean or :median

# Normalisation uses a first-pass sample over all months in the bbox.
normalisation_lower_q = 0.02
normalisation_upper_q = 0.98
normalisation_sample_stride = 12  # increase for faster approximate references

# Profile.
selected_profile = fish_profile()

# Output.
output_dir = joinpath(@__DIR__, "..", "output", "dynqual_country_seasonality_animation")
mkpath(output_dir)

make_timeseries_plot = true
make_seasonality_plot = true
make_animation = true
animation_filename = joinpath(output_dir, "europe_K_$(String(K_metric))_1980_2019.mp4")
animation_framerate = 12

# Plot style.
show_coastlines = true
show_country_borders = true
show_map_border = true
show_geo_grid = false
coastline_color = :black
coastline_linewidth = 0.65
country_border_color = (:black, 0.55)
country_border_linewidth = 0.35
map_border_color = :black
map_border_linewidth = 1.2
geo_source = "+proj=longlat +datum=WGS84"
geo_dest = "+proj=longlat +datum=WGS84"

# Fixed color range for animation; if nothing, it is inferred after the first pass.
K_colormap = K_metric == :amplification ? :inferno : :viridis
K_colorrange = nothing

# ==============================================================================
# NETCDF HELPERS
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
        named !== nothing && return named

        coord_like = Set(["lon", "longitude", "lat", "latitude", "time", "time_bnds", "bounds", "crs"])
        possible = String[]
        for k in keys(ds)
            name = String(k)
            if lowercase(name) in coord_like
                continue
            end
            v = ds[name]
            nd = try ndims(v) catch; 0 end
            nd >= 2 && push!(possible, name)
        end
        isempty(possible) && error("Could not infer gridded variable in $path. Variables: $(collect(keys(ds)))")
        length(possible) > 1 && @warn "Multiple possible data variables found; using first" path=path candidates=possible chosen=possible[1]
        return possible[1]
    finally
        close(ds)
    end
end

function lon_lat_names(ds)
    lon_name = first_existing_key(ds, ["lon", "longitude", "x"])
    lat_name = first_existing_key(ds, ["lat", "latitude", "y"])
    lon_name === nothing && error("No lon/longitude/x variable found. Variables: $(collect(keys(ds)))")
    lat_name === nothing && error("No lat/latitude/y variable found. Variables: $(collect(keys(ds)))")
    return lon_name, lat_name
end

function selected_lon_lat_indices(lon, lat, bbox)
    br = bbox_ranges(bbox)
    br === nothing && return collect(1:length(lon)), collect(1:length(lat))

    lon_min, lon_max, lat_min, lat_max = br.lon_min, br.lon_max, br.lat_min, br.lat_max

    if minimum(skipmissing(lon)) >= 0 && lon_min < 0
        lon_min_mod = mod(lon_min, 360)
        lon_max_mod = mod(lon_max, 360)
        lon_idx = lon_min_mod > lon_max_mod ?
            findall(x -> (x >= lon_min_mod) || (x <= lon_max_mod), lon) :
            findall(x -> lon_min_mod <= x <= lon_max_mod, lon)
    else
        lon_idx = findall(x -> lon_min <= x <= lon_max, lon)
    end
    lat_idx = findall(x -> lat_min <= x <= lat_max, lat)

    isempty(lon_idx) && error("bbox selected zero longitude cells. File lon range $(minimum(lon))..$(maximum(lon)); bbox $((lon_min, lon_max))")
    isempty(lat_idx) && error("bbox selected zero latitude cells. File lat range $(minimum(lat))..$(maximum(lat)); bbox $((lat_min, lat_max))")
    return lon_idx, lat_idx
end

function load_nc_layer_auto(path, candidates, time_index, bbox)
    isfile(path) || error("File not found: $path")
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
            @warn "No dimension names; assuming lon, lat, time" variable=varname path=path
            data = v[lon_idx, lat_idx, time_index]
            return replace_missing_with_nan(data), lon[lon_idx], lat[lat_idx], varname
        end

        idxs = Any[]
        kept_dims = String[]
        for d in dims
            dl = lowercase(d)
            if occursin("lon", dl) || dl == "x"
                push!(idxs, lon_idx); push!(kept_dims, "lon")
            elseif occursin("lat", dl) || dl == "y"
                push!(idxs, lat_idx); push!(kept_dims, "lat")
            elseif occursin("time", dl)
                push!(idxs, time_index)
            else
                push!(idxs, Colon()); push!(kept_dims, d)
            end
        end

        data = replace_missing_with_nan(Array(v[idxs...]))
        singleton_dims = Tuple(findall(size(data) .== 1))
        if !isempty(singleton_dims) && ndims(data) > 2
            data = dropdims(data; dims=singleton_dims)
        end
        ndims(data) != 2 && error("Loaded $varname from $path but result is not 2D. Size: $(size(data)), dims: $dims")

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

function sort_for_plot(layer, lon, lat)
    lon_plot = Float64.(lon)
    lon_plot = [x > 180 ? x - 360 : x for x in lon_plot]
    lat_plot = Float64.(lat)
    lon_order = sortperm(lon_plot)
    lat_order = sortperm(lat_plot)
    return copy(layer)[lon_order, lat_order], lon_plot[lon_order], lat_plot[lat_order]
end

# ==============================================================================
# NORMALISATION AND PIPELINE
# ==============================================================================

function sample_values!(store::Vector{Float64}, A; stride=12)
    vals = finite_values(A)
    isempty(vals) && return store
    append!(store, vals[1:stride:end])
    return store
end

function robust_ref_from_store(vals; lower_q=0.02, upper_q=0.98)
    isempty(vals) && return (lo=NaN, hi=NaN)
    return (lo=quantile(vals, lower_q), hi=quantile(vals, upper_q))
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

function load_raw_layers_for_time(year, month)
    t_idx = time_index_from_year_month(year, month; base_year=base_year)
    WT, lon, lat, WT_var = load_nc_layer_auto(WT_file, WT_var_candidates, t_idx, bbox)
    BOD, _, _, BOD_var = load_nc_layer_auto(BOD_file, BOD_var_candidates, t_idx, bbox)
    FC, _, _, FC_var = load_nc_layer_auto(FC_file, FC_var_candidates, t_idx, bbox)
    TDS, _, _, TDS_var = load_nc_layer_auto(TDS_file, TDS_var_candidates, t_idx, bbox)

    ref_shape = size(WT)
    for (nm, layer) in [("BOD", BOD), ("FC", FC), ("TDS", TDS)]
        size(layer) != ref_shape && error("Layer $nm has size $(size(layer)), expected $ref_shape")
    end

    zero_layer = fill(0.0, ref_shape)
    return (
        WT=WT, BOD=BOD, TDS=TDS, FC=FC,
        Nutrient=zero_layer, Chemical=zero_layer, Plastic=zero_layer,
        lon=lon, lat=lat,
        metadata=(WT_var=WT_var, BOD_var=BOD_var, FC_var=FC_var, TDS_var=TDS_var)
    )
end

function compute_normalisation_refs()
    println("\n--- First pass: sampling data for shared robust normalisation refs ---")
    stores = Dict(:WT=>Float64[], :BOD=>Float64[], :TDS=>Float64[], :FC=>Float64[])
    lon = nothing
    lat = nothing
    metadata = nothing

    for (idx, (yr, mo)) in enumerate(all_times)
        idx % 24 == 1 && println("Sampling refs: $yr-$(lpad(mo, 2, "0")) [$idx / $(length(all_times))]")
        raw = load_raw_layers_for_time(yr, mo)
        sample_values!(stores[:WT], raw.WT; stride=normalisation_sample_stride)
        sample_values!(stores[:BOD], raw.BOD; stride=normalisation_sample_stride)
        sample_values!(stores[:TDS], raw.TDS; stride=normalisation_sample_stride)
        sample_values!(stores[:FC], raw.FC; stride=normalisation_sample_stride)
        lon === nothing && (lon = raw.lon)
        lat === nothing && (lat = raw.lat)
        metadata === nothing && (metadata = raw.metadata)
    end

    refs = (
        WT = robust_ref_from_store(stores[:WT]; lower_q=normalisation_lower_q, upper_q=normalisation_upper_q),
        BOD = robust_ref_from_store(stores[:BOD]; lower_q=normalisation_lower_q, upper_q=normalisation_upper_q),
        TDS = robust_ref_from_store(stores[:TDS]; lower_q=normalisation_lower_q, upper_q=normalisation_upper_q),
        FC = robust_ref_from_store(stores[:FC]; lower_q=normalisation_lower_q, upper_q=normalisation_upper_q),
        Nutrient = (lo=0.0, hi=1.0),
        Chemical = (lo=0.0, hi=1.0),
        Plastic = (lo=0.0, hi=1.0)
    )

    println("Detected variables: ", metadata)
    println("Normalisation refs: ", refs)
    return refs, lon, lat
end

function compute_K_grid(year, month, refs)
    raw = load_raw_layers_for_time(year, month)
    nrm = (
        WT = apply_reference_normalisation(raw.WT, refs.WT),
        BOD = apply_reference_normalisation(raw.BOD, refs.BOD),
        TDS = apply_reference_normalisation(raw.TDS, refs.TDS),
        FC = apply_reference_normalisation(raw.FC, refs.FC),
        Nutrient = raw.Nutrient,
        Chemical = raw.Chemical,
        Plastic = raw.Plastic
    )
    layers_norm = [nrm.WT, nrm.BOD, nrm.TDS, nrm.FC, nrm.Nutrient, nrm.Chemical, nrm.Plastic]
    res = isimip_deb_pipeline_grid(
        layers_norm,
        selected_profile.exposure_filter,
        selected_profile.moa_mapping,
        selected_profile.moa_deb_mapping,
        selected_profile.deb_params
    )
    return get_K_grid(res), raw.lon, raw.lat
end

function get_K_grid(res)
    if K_metric == :amplification
        return res.amplification
    elseif K_metric == :lambda
        return res.lambda
    elseif K_metric == :A
        return res.A
    elseif K_metric in (:thermal, :oxygen, :osmotic, :immune, :eutrophication, :toxic, :feeding, :physical)
        return getproperty(res.modes, K_metric)
    elseif K_metric in (:assimilation, :maintenance, :growth, :reproduction)
        return getproperty(res.axes, K_metric)
    else
        error("Unknown K_metric: $K_metric")
    end
end

# ==============================================================================
# COUNTRY MASKS
# ==============================================================================

function row_name_candidates(row)
    names = String[]
    for p in propertynames(row)
        ps = String(p)
        pls = lowercase(ps)
        if pls in ["name", "name_en", "admin", "sovereignt", "brk_name", "formal_en", "geounit"]
            try
                val = getproperty(row, p)
                if val !== missing && val !== nothing
                    push!(names, String(val))
                end
            catch
            end
        end
    end
    return names
end

function find_country_rows(countries, wanted_names)
    rows = collect(Tables.rows(countries))
    found = Dict{String, Any}()
    for wanted in wanted_names
        wl = lowercase(wanted)
        matches = Any[]
        for row in rows
            vals = row_name_candidates(row)
            if any(v -> occursin(wl, lowercase(v)) || occursin(lowercase(v), wl), vals)
                push!(matches, row)
            end
        end
        if isempty(matches)
            @warn "No NaturalEarth country match found" wanted=wanted
        elseif length(matches) > 1
            @warn "Multiple NaturalEarth country matches; using first" wanted=wanted matches=[row_name_candidates(m) for m in matches]
            found[wanted] = matches[1]
        else
            found[wanted] = matches[1]
        end
    end
    return found
end

is_pointlike(x) = try
    length(x) >= 2 && x[1] isa Real && x[2] isa Real
catch
    false
end

function getxy(p)
    return Float64(p[1]), Float64(p[2])
end

function point_in_ring(x, y, ring)
    inside = false
    n = length(ring)
    n < 3 && return false
    xj, yj = getxy(ring[end])
    for i in 1:n
        xi, yi = getxy(ring[i])
        intersects = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / ((yj - yi) + eps()) + xi)
        if intersects
            inside = !inside
        end
        xj, yj = xi, yi
    end
    return inside
end

function polygon_contains(x, y, polygon)
    isempty(polygon) && return false
    outer = polygon[1]
    point_in_ring(x, y, outer) || return false
    # Holes: if point is inside a hole, it is outside the polygon.
    for h in polygon[2:end]
        point_in_ring(x, y, h) && return false
    end
    return true
end

function coords_as_polygons(coords)
    # Returns vector of polygons, where each polygon is vector of rings.
    if isempty(coords)
        return []
    end
    if is_pointlike(coords[1])
        # A single ring was supplied; wrap as one polygon.
        return [[coords]]
    elseif !isempty(coords[1]) && is_pointlike(coords[1][1])
        # Polygon: vector of rings.
        return [coords]
    else
        # MultiPolygon: vector of polygons.
        return coords
    end
end

function geometry_contains_lonlat(geom, lon, lat)
    coords = GeoInterface.coordinates(geom)
    polygons = coords_as_polygons(coords)
    for poly in polygons
        polygon_contains(lon, lat, poly) && return true
    end
    return false
end

function build_country_masks(lon, lat)
    println("\n--- Building country masks from NaturalEarth ---")
    countries = naturalearth("admin_0_countries", 110)
    country_rows = find_country_rows(countries, selected_countries)

    layer_dummy = zeros(length(lon), length(lat))
    _, lon_s, lat_s = sort_for_plot(layer_dummy, lon, lat)

    masks = Dict{String, BitMatrix}()
    for country in selected_countries
        if !haskey(country_rows, country)
            continue
        end
        row = country_rows[country]
        geom = getproperty(row, :geometry)
        mask = falses(length(lon_s), length(lat_s))
        for i in eachindex(lon_s), j in eachindex(lat_s)
            mask[i, j] = geometry_contains_lonlat(geom, lon_s[i], lat_s[j])
        end
        ncells = count(mask)
        println("Country mask: $country -> $ncells grid cells")
        if ncells == 0
            @warn "Country mask has zero cells; check country name or bbox" country=country
        end
        masks[country] = mask
    end
    return masks, lon_s, lat_s
end

function country_stat(K_sorted, mask)
    vals = K_sorted[mask]
    vals = filter(isfinite, vec(vals))
    isempty(vals) && return NaN
    if country_statistic == :median
        return median(vals)
    else
        return mean(vals)
    end
end

# ==============================================================================
# MAP PLOTTING / ANIMATION
# ==============================================================================

function plot_country_borders!(ax)
    if show_country_borders
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
            end
        catch err
            @warn "Could not plot country borders" exception=(err, catch_backtrace())
        end
    end
end

function plot_map_border!(ax, br)
    if show_map_border && br !== nothing
        border_lon = [br.lon_min, br.lon_max, br.lon_max, br.lon_min, br.lon_min]
        border_lat = [br.lat_min, br.lat_min, br.lat_max, br.lat_max, br.lat_min]
        lines!(ax, border_lon, border_lat; color=map_border_color, linewidth=map_border_linewidth, overdraw=true)
    end
end

function setup_animation_figure(K0, lon_s, lat_s; colorrange)
    br = bbox_ranges(bbox)
    fig = Figure(size=(1000, 700))
    title_obs = Observable("$(K_label): initialising")
    Kobs = Observable(K0)

    ax = GeoAxis(
        fig[1, 1];
        source = geo_source,
        dest = geo_dest,
        limits = (br.lon_min, br.lon_max, br.lat_min, br.lat_max),
        title = title_obs,
        xlabel = "Longitude",
        ylabel = "Latitude",
        #backgroundcolor = :white,
        xgridvisible = show_geo_grid,
        ygridvisible = show_geo_grid
    )

    hm = heatmap!(ax, lon_s, lat_s, Kobs; colormap=K_colormap, colorrange=colorrange, nan_color=:transparent)
    plot_country_borders!(ax)
    if show_coastlines
        lines!(ax, GeoMakie.coastlines(); color=coastline_color, linewidth=coastline_linewidth, overdraw=true)
    end
    plot_map_border!(ax, br)
    xlims!(ax, br.lon_min, br.lon_max)
    ylims!(ax, br.lat_min, br.lat_max)
    Colorbar(fig[1, 2], hm; label=K_label)
    return fig, Kobs, title_obs
end

# ==============================================================================
# TIME SERIES PLOTS
# ==============================================================================

function month_name(m)
    return ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][m]
end

function plot_country_timeseries(rows)
    fig = Figure(size=(1200, 650))
    ax = Axis(fig[1, 1]; title="Country $(String(country_statistic)) $(K_label) over time", xlabel="Year", ylabel=K_label)
    for country in selected_countries
        vals = [r.K for r in rows if r.country == country]
        years = [r.year + (r.month - 1) / 12 for r in rows if r.country == country]
        isempty(vals) && continue
        lines!(ax, years, vals; label=country)
    end
    axislegend(ax; position=:rt)
    save(joinpath(output_dir, "country_K_timeseries_$(String(K_metric)).png"), fig, px_per_unit=2)
    return fig
end

function plot_country_seasonality(rows)
    fig = Figure(size=(1000, 650))
    ax = Axis(fig[1, 1]; title="Monthly seasonality of country $(String(country_statistic)) $(K_label), $(start_year)-$(end_year)", xlabel="Month", ylabel=K_label, xticks=(1:12, month_name.(1:12)))
    for country in selected_countries
        ys = Float64[]
        for m in 1:12
            vals = [r.K for r in rows if r.country == country && r.month == m && isfinite(r.K)]
            push!(ys, isempty(vals) ? NaN : mean(vals))
        end
        all(!isfinite, ys) && continue
        lines!(ax, 1:12, ys; label=country)
        scatter!(ax, 1:12, ys; markersize=6)
    end
    axislegend(ax; position=:rt)
    save(joinpath(output_dir, "country_K_seasonality_$(String(K_metric)).png"), fig, px_per_unit=2)
    return fig
end

# ==============================================================================
# MAIN
# ==============================================================================

println("\n--- DynQual country seasonality + Europe animation ---")
println("K metric: $K_metric")
println("Countries: ", selected_countries)
println("bbox: $bbox")
println("Output: $output_dir")

refs, lon, lat = compute_normalisation_refs()

# Build country masks using sorted plotting lon/lat order.
country_masks, lon_s, lat_s = build_country_masks(lon, lat)

# Initial K for color range and animation setup.
K0_raw, _, _ = compute_K_grid(start_year, first(months), refs)
K0, lon_s2, lat_s2 = sort_for_plot(K0_raw, lon, lat)

if K_colorrange === nothing
    # Sample K over all months to get robust animation color range.
    println("\n--- Sampling K for animation color range ---")
    kvals = Float64[]
    for (idx, (yr, mo)) in enumerate(all_times)
        idx % 24 == 1 && println("Sampling K color range: $yr-$(lpad(mo, 2, "0")) [$idx / $(length(all_times))]")
        Kraw, _, _ = compute_K_grid(yr, mo, refs)
        Ksorted, _, _ = sort_for_plot(Kraw, lon, lat)
        sample_values!(kvals, Ksorted; stride=normalisation_sample_stride)
    end
    if isempty(kvals)
        K_colorrange = (0.0, 1.0)
    else
        if K_metric == :amplification
            K_colorrange = (quantile(kvals, 0.02), quantile(kvals, 0.98))
        else
            K_colorrange = (quantile(kvals, 0.02), quantile(kvals, 0.98))
        end
    end
end
println("K color range: ", K_colorrange)

country_rows = NamedTuple[]

function process_one_time!(yr, mo, Kobs=nothing, title_obs=nothing)
    Kraw, _, _ = compute_K_grid(yr, mo, refs)
    Ksorted, _, _ = sort_for_plot(Kraw, lon, lat)

    for country in selected_countries
        haskey(country_masks, country) || continue
        Kc = country_stat(Ksorted, country_masks[country])
        push!(country_rows, (year=yr, month=mo, country=country, K=Kc, metric=String(K_metric), statistic=String(country_statistic)))
    end

    if Kobs !== nothing
        Kobs[] = Ksorted
        title_obs[] = "$(K_label): $(yr)-$(lpad(mo, 2, "0"))"
    end
    return Ksorted
end

if make_animation
    fig, Kobs, title_obs = setup_animation_figure(K0, lon_s, lat_s; colorrange=K_colorrange)
    println("\n--- Recording animation: $animation_filename ---")
    record(fig, animation_filename, 1:length(all_times); framerate=animation_framerate) do frame
        yr, mo = all_times[frame]
        println("Animation frame $frame / $(length(all_times)): $yr-$(lpad(mo, 2, "0"))")
        process_one_time!(yr, mo, Kobs, title_obs)
    end
else
    println("\n--- Computing country time series without animation ---")
    for (idx, (yr, mo)) in enumerate(all_times)
        idx % 12 == 1 && println("Processing $yr [$idx / $(length(all_times))]")
        process_one_time!(yr, mo)
    end
end

if make_timeseries_plot
    plot_country_timeseries(country_rows)
end
if make_seasonality_plot
    plot_country_seasonality(country_rows)
end

# Save a lightweight CSV summary if CSV/DataFrames are available. This is not a grid;
# it is useful for checking the plotted country K values.
try
    @eval using DataFrames
    @eval using CSV
    df = DataFrame(country_rows)
    CSV.write(joinpath(output_dir, "country_K_values_$(String(K_metric)).csv"), df)
catch err
    @info "CSV/DataFrames not available; skipping CSV output" exception=(err, catch_backtrace())
end

println("\nDone.")
println("Plots written to: $output_dir")
println("Animation: ", make_animation ? animation_filename : "not generated")
