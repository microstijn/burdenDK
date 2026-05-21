# ==============================================================================
# dynqual_country_seasonality_multipanel_animation.jl
#
# Secondary real-data script for DynQual / ISIMIP-style monthly NetCDF files.
#
# Goals:
#   1. Compute monthly Europe-wide MoA/DEB outputs for 1980-2019.
#   2. Aggregate K by selected countries and plot country seasonality/time series.
#   3. Generate a multi-panel Europe animation:
#        - small panels: non-zero input layers WT, BOD/organic, TDSload, pathogen
#        - large panels: K/amplification and adaptive margin A
#        - extra informative panels: lambda/recovery and maintenance burden
#
# K is defined by default as the amplification factor F:
#   K = F = lambda_control / lambda_background
#
# No raster grids are saved. Only plots, optional CSV, and animation are written.
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
using NaturalEarth
using Tables

# ==============================================================================
# USER SETTINGS
# ==============================================================================

WT_file       = raw"C:\Users\peete074\Downloads\waterTemperature_monthlyAvg_1980_2019.nc"
BOD_file      = raw"C:\Users\peete074\Downloads\organic_monthlyAvg_1980_2019.nc"
FC_file       = raw"C:\Users\peete074\Downloads\pathogen_monthlyAvg_1980_2019.nc"
TDS_file      = raw"C:\Users\peete074\Downloads\TDSload_monthlyAvg_1980_2019.nc"
BODload_file  = raw"C:\Users\peete074\Downloads\BODload_monthlyAvg_1980_2019.nc"

WT_var_candidates  = ["waterTemperature", "watertemperature", "water_temperature", "waterTemp", "watertemp", "WT", "temperature", "triver"]
BOD_var_candidates = ["organic", "BOD", "bod", "bod_concentration", "biological_oxygen_demand"]
FC_var_candidates  = ["pathogen", "FC", "fc", "fecal_coliform", "faecal_coliform"]
TDS_var_candidates = ["TDSload", "tdsload", "TDSLoad", "tds_load", "TDS", "tds"]

bbox = (lon = (-25.0, 45.0), lat = (34.0, 72.0))

start_year = 1980
end_year = 2019
base_year = 1980
months = 1:12
all_times = [(y, m) for y in start_year:end_year for m in months]

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

# Country summary metric. The map animation will always show K, A, lambda,
# maintenance, and input layers.
K_metric = :amplification
K_label = "K = amplification factor F"
country_statistic = :mean  # :mean or :median

normalisation_lower_q = 0.02
normalisation_upper_q = 0.98
normalisation_sample_stride = 12

selected_profile = fish_profile()

output_dir = joinpath(@__DIR__, "..", "output", "dynqual_country_seasonality_multipanel_animation")
mkpath(output_dir)

make_timeseries_plot = true
make_seasonality_plot = true
make_animation = true
animation_filename = joinpath(output_dir, "europe_multipanel_K_A_inputs_1980_2019.mp4")
animation_framerate = 12

# Plot style.
show_coastlines = true
show_country_borders = true
show_map_border = true
show_geo_grid = false
coastline_color = :black
coastline_linewidth = 0.45
country_border_color = (:black, 0.45)
country_border_linewidth = 0.25
map_border_color = :black
map_border_linewidth = 0.8
geo_source = "+proj=longlat +datum=WGS84"
geo_dest = "+proj=longlat +datum=WGS84"

# Fixed animation ranges. Inputs are normalised to 0..1. A/lambda/maintenance are
# usually in or near 0..1. K is inferred robustly unless set manually.
input_colormap = :viridis
K_colormap = :inferno
A_colormap = :viridis
lambda_colormap = :viridis
maintenance_colormap = :plasma
K_colorrange = nothing
A_colorrange = (0.0, 1.0)
lambda_colorrange = (0.0, 1.0)
maintenance_colorrange = (0.0, 1.0)
input_colorrange = (0.0, 1.0)

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

replace_missing_with_nan(A) = Float64.(coalesce.(A, NaN))

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
            lowercase(name) in coord_like && continue
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
    isempty(lon_idx) && error("bbox selected zero longitude cells. File lon range $(minimum(lon))..$(maximum(lon)); bbox ($(lon_min), $(lon_max))")
    isempty(lat_idx) && error("bbox selected zero latitude cells. File lat range $(minimum(lat))..$(maximum(lat)); bbox ($(lat_min), $(lat_max))")
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
# NORMALISATION / PIPELINE
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

function run_pipeline_for_time(year, month, refs)
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
    return nrm, res, raw.lon, raw.lat
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

const PRIORITY_NAME_FIELDS = [:admin, :ADMIN, :name, :NAME, :name_en, :NAME_EN, :brk_name, :BRK_NAME, :geounit, :GEOUNIT]
const SECONDARY_NAME_FIELDS = [:formal_en, :FORMAL_EN, :name_long, :NAME_LONG]
const SOVEREIGNTY_FIELDS = [:sovereignt, :SOVEREIGNT]

function safe_get(row, p)
    try
        if p in propertynames(row)
            v = getproperty(row, p)
            if v !== missing && v !== nothing
                return String(v)
            end
        end
    catch
    end
    return nothing
end

function vals_for_fields(row, fields)
    out = String[]
    for f in fields
        v = safe_get(row, f)
        v === nothing || push!(out, v)
    end
    return unique(out)
end

function row_name_candidates(row)
    return unique(vcat(vals_for_fields(row, PRIORITY_NAME_FIELDS), vals_for_fields(row, SECONDARY_NAME_FIELDS), vals_for_fields(row, SOVEREIGNTY_FIELDS)))
end

function country_match_score(row, wanted)
    wl = lowercase(wanted)
    priority = vals_for_fields(row, PRIORITY_NAME_FIELDS)
    secondary = vals_for_fields(row, SECONDARY_NAME_FIELDS)
    sovereign = vals_for_fields(row, SOVEREIGNTY_FIELDS)
    any(v -> lowercase(v) == wl, priority) && return 100
    any(v -> occursin(wl, lowercase(v)) || occursin(lowercase(v), wl), priority) && return 80
    any(v -> lowercase(v) == wl, secondary) && return 70
    any(v -> occursin(wl, lowercase(v)) || occursin(lowercase(v), wl), secondary) && return 60
    any(v -> lowercase(v) == wl, sovereign) && return 10
    any(v -> occursin(wl, lowercase(v)) || occursin(lowercase(v), wl), sovereign) && return 5
    return 0
end

function find_country_rows(countries, wanted_names)
    rows = collect(Tables.rows(countries))
    found = Dict{String, Any}()
    for wanted in wanted_names
        scored = [(score=country_match_score(row, wanted), row=row, names=row_name_candidates(row)) for row in rows]
        scored = filter(x -> x.score > 0, scored)
        if isempty(scored)
            @warn "No NaturalEarth country match found" wanted=wanted
            continue
        end
        sort!(scored, by=x -> x.score, rev=true)
        best = scored[1]
        println("Matched country '$wanted' to NaturalEarth row: ", best.names, " with score ", best.score)
        found[wanted] = best.row
    end
    return found
end

is_pointlike(x) = try
    length(x) >= 2 && x[1] isa Real && x[2] isa Real
catch
    false
end
getxy(p) = (Float64(p[1]), Float64(p[2]))

function point_in_ring(x, y, ring)
    inside = false
    n = length(ring)
    n < 3 && return false
    xj, yj = getxy(ring[end])
    for i in 1:n
        xi, yi = getxy(ring[i])
        intersects = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / ((yj - yi) + eps()) + xi)
        intersects && (inside = !inside)
        xj, yj = xi, yi
    end
    return inside
end

function polygon_contains(x, y, polygon)
    isempty(polygon) && return false
    outer = polygon[1]
    point_in_ring(x, y, outer) || return false
    for h in polygon[2:end]
        point_in_ring(x, y, h) && return false
    end
    return true
end

function coords_as_polygons(coords)
    isempty(coords) && return []
    if is_pointlike(coords[1])
        return [[coords]]
    elseif !isempty(coords[1]) && is_pointlike(coords[1][1])
        return [coords]
    else
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
        haskey(country_rows, country) || continue
        row = country_rows[country]
        geom = getproperty(row, :geometry)
        mask = falses(length(lon_s), length(lat_s))
        for i in eachindex(lon_s), j in eachindex(lat_s)
            mask[i, j] = geometry_contains_lonlat(geom, lon_s[i], lat_s[j])
        end
        ncells = count(mask)
        println("Country mask: $country -> $ncells grid cells")
        ncells == 0 && @warn "Country mask has zero cells" country=country chosen=row_name_candidates(row)
        masks[country] = mask
    end
    return masks, lon_s, lat_s
end

function country_stat(K_sorted, mask)
    vals = K_sorted[mask]
    vals = filter(isfinite, vec(vals))
    isempty(vals) && return NaN
    country_statistic == :median ? median(vals) : mean(vals)
end

# ==============================================================================
# MAP HELPERS
# ==============================================================================

function plot_country_borders!(ax)
    show_country_borders || return nothing
    try
        countries = naturalearth("admin_0_countries", 110)
        p = poly!(ax, GeoMakie.to_multipoly(countries.geometry); color=(:white, 0.0), strokecolor=country_border_color, strokewidth=country_border_linewidth, overdraw=true)
        try
            translate!(p, 0, 0, 100)
        catch
        end
    catch err
        @warn "Could not plot country borders" exception=(err, catch_backtrace())
    end
end

function plot_map_border!(ax, br)
    if show_map_border && br !== nothing
        border_lon = [br.lon_min, br.lon_max, br.lon_max, br.lon_min, br.lon_min]
        border_lat = [br.lat_min, br.lat_min, br.lat_max, br.lat_max, br.lat_min]
        lines!(ax, border_lon, border_lat; color=map_border_color, linewidth=map_border_linewidth, overdraw=true)
    end
end

function add_map_panel!(fig, cell, data_obs, lon_s, lat_s; title, colormap, colorrange, label, title_obs=nothing)
    br = bbox_ranges(bbox)
    ax = GeoAxis(
        fig[cell...];
        source=geo_source,
        dest=geo_dest,
        limits=(br.lon_min, br.lon_max, br.lat_min, br.lat_max),
        title=title_obs === nothing ? title : title_obs,
        xlabel="",
        ylabel="",
        #backgroundcolor=:white,
        xgridvisible=show_geo_grid,
        ygridvisible=show_geo_grid,
        xticklabelsvisible=false,
        yticklabelsvisible=false
    )
    hm = heatmap!(ax, lon_s, lat_s, data_obs; colormap=colormap, colorrange=colorrange, nan_color=:transparent)
    plot_country_borders!(ax)
    show_coastlines && lines!(ax, GeoMakie.coastlines(); color=coastline_color, linewidth=coastline_linewidth, overdraw=true)
    plot_map_border!(ax, br)
    xlims!(ax, br.lon_min, br.lon_max)
    ylims!(ax, br.lat_min, br.lat_max)
    Colorbar(fig[cell[1], cell[2] + 1], hm; label=label, width=10, ticklabelsize=8, labelsize=9)
    return ax, hm
end

function prepare_panels(year, month, refs, lon, lat)
    nrm, res, _, _ = run_pipeline_for_time(year, month, refs)
    K = get_K_grid(res)
    panels = (
        WT = nrm.WT,
        BOD = nrm.BOD,
        TDS = nrm.TDS,
        FC = nrm.FC,
        K = K,
        A = res.A,
        lambda = res.lambda,
        maintenance = res.axes.maintenance
    )
    sorted = map(x -> sort_for_plot(x, lon, lat)[1], panels)
    return sorted
end

function sample_K_colorrange(refs, lon, lat)
    if K_colorrange !== nothing
        return K_colorrange
    end
    println("\n--- Sampling K for animation color range ---")
    kvals = Float64[]
    for (idx, (yr, mo)) in enumerate(all_times)
        idx % 24 == 1 && println("Sampling K color range: $yr-$(lpad(mo, 2, "0")) [$idx / $(length(all_times))]")
        _, res, _, _ = run_pipeline_for_time(yr, mo, refs)
        K = get_K_grid(res)
        Ksorted, _, _ = sort_for_plot(K, lon, lat)
        sample_values!(kvals, Ksorted; stride=normalisation_sample_stride)
    end
    return isempty(kvals) ? (1.0, 2.0) : (quantile(kvals, 0.02), quantile(kvals, 0.98))
end

# ==============================================================================
# TIME SERIES PLOTS
# ==============================================================================

month_name(m) = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][m]

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
end

# ==============================================================================
# MAIN
# ==============================================================================

println("\n--- DynQual country seasonality + multi-panel Europe animation ---")
println("K metric: $K_metric")
println("Countries: ", selected_countries)
println("bbox: $bbox")
println("Output: $output_dir")

refs, lon, lat = compute_normalisation_refs()
country_masks, lon_s, lat_s = build_country_masks(lon, lat)
K_range = sample_K_colorrange(refs, lon, lat)
println("K color range: ", K_range)

initial = prepare_panels(start_year, first(months), refs, lon, lat)
obs = Dict(
    :WT => Observable(initial.WT),
    :BOD => Observable(initial.BOD),
    :TDS => Observable(initial.TDS),
    :FC => Observable(initial.FC),
    :K => Observable(initial.K),
    :A => Observable(initial.A),
    :lambda => Observable(initial.lambda),
    :maintenance => Observable(initial.maintenance)
)
title_obs = Observable("Europe vulnerability dashboard: $(start_year)-$(lpad(first(months), 2, "0"))")

country_rows = NamedTuple[]

function process_one_time!(yr, mo; update_observables=false)
    panels = prepare_panels(yr, mo, refs, lon, lat)

    for country in selected_countries
        haskey(country_masks, country) || continue
        Kc = country_stat(panels.K, country_masks[country])
        push!(country_rows, (year=yr, month=mo, country=country, K=Kc, metric=String(K_metric), statistic=String(country_statistic)))
    end

    if update_observables
        obs[:WT][] = panels.WT
        obs[:BOD][] = panels.BOD
        obs[:TDS][] = panels.TDS
        obs[:FC][] = panels.FC
        obs[:K][] = panels.K
        obs[:A][] = panels.A
        obs[:lambda][] = panels.lambda
        obs[:maintenance][] = panels.maintenance
        title_obs[] = "Europe vulnerability dashboard: $(yr)-$(lpad(mo, 2, "0"))"
    end
    return panels
end

if make_animation
    fig = Figure(size=(1700, 1100))
    Label(fig[0, 1:10], title_obs; fontsize=24, tellwidth=false)

    # Four small input panels. Each panel occupies two columns: map + colorbar.
    add_map_panel!(fig, (1, 1), obs[:WT], lon_s, lat_s; title="Water temperature", colormap=input_colormap, colorrange=input_colorrange, label="norm.")
    add_map_panel!(fig, (1, 3), obs[:BOD], lon_s, lat_s; title="Organic / BOD", colormap=input_colormap, colorrange=input_colorrange, label="norm.")
    add_map_panel!(fig, (1, 5), obs[:TDS], lon_s, lat_s; title="TDSload proxy", colormap=input_colormap, colorrange=input_colorrange, label="norm.")
    add_map_panel!(fig, (1, 7), obs[:FC], lon_s, lat_s; title="Pathogen / FC", colormap=input_colormap, colorrange=input_colorrange, label="norm.")

    # Two large panels: K and adaptive margin.
    add_map_panel!(fig, (2, 1), obs[:K], lon_s, lat_s; title="K: amplification", colormap=K_colormap, colorrange=K_range, label="K")
    add_map_panel!(fig, (2, 5), obs[:A], lon_s, lat_s; title="Adaptive margin A", colormap=A_colormap, colorrange=A_colorrange, label="A")

    # Extra informative panels: lambda and maintenance. These explain K mechanistically.
    add_map_panel!(fig, (3, 1), obs[:lambda], lon_s, lat_s; title="Recovery λ", colormap=lambda_colormap, colorrange=lambda_colorrange, label="λ")
    add_map_panel!(fig, (3, 5), obs[:maintenance], lon_s, lat_s; title="Maintenance burden", colormap=maintenance_colormap, colorrange=maintenance_colorrange, label="s_M")

    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 8)

    println("\n--- Recording multi-panel animation: $animation_filename ---")
    record(fig, animation_filename, 1:length(all_times); framerate=animation_framerate) do frame
        yr, mo = all_times[frame]
        println("Animation frame $frame / $(length(all_times)): $yr-$(lpad(mo, 2, "0"))")
        process_one_time!(yr, mo; update_observables=true)
    end
else
    println("\n--- Computing country time series without animation ---")
    for (idx, (yr, mo)) in enumerate(all_times)
        idx % 12 == 1 && println("Processing $yr [$idx / $(length(all_times))]")
        process_one_time!(yr, mo; update_observables=false)
    end
end

make_timeseries_plot && plot_country_timeseries(country_rows)
make_seasonality_plot && plot_country_seasonality(country_rows)

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
