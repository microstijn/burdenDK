# ==============================================================================
# dynqual_K_A_country_trails_animation.jl
#
# Focused monthly animation for DynQual / ISIMIP-style monthly NetCDF files.
#
# Animation layout:
#   Panel 1: Europe map of K = amplification factor F
#   Panel 2: Europe map of adaptive margin A
#   Panel 3: Country seasonal trails, with x = month-of-year and y = country K
#
# Country trails:
#   - Each country has its own color.
#   - Previous completed years are drawn as faint trajectories in that country color.
#   - The currently animated year is drawn as a stronger trajectory.
#   - The current month is marked with a larger dot.
#
# Country aggregation defaults to the median K within the country mask. Median is
# used because K can be spatially spiky when recovery λ becomes small; the median
# better reflects the typical country-level background than the mean. Change
# country_statistic below if a hotspot-sensitive measure is wanted.
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

WT_var_candidates  = ["waterTemperature", "watertemperature", "water_temperature", "waterTemp", "watertemp", "WT", "temperature", "triver"]
BOD_var_candidates = ["organic", "BOD", "bod", "bod_concentration", "biological_oxygen_demand"]
FC_var_candidates  = ["pathogen", "FC", "fc", "fecal_coliform", "faecal_coliform"]
TDS_var_candidates = ["TDSload", "tdsload", "TDSLoad", "tds_load", "TDS", "tds"]

bbox = (lon = (-25.0, 45.0), lat = (34.0, 72.0))

start_year = 1980
end_year = 2019
base_year = 1980
years = collect(start_year:end_year)
months = collect(1:12)
all_times = [(y, m) for y in years for m in months]

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

# Country K statistic.
# Recommended default: :median because K can be hotspot/spike dominated.
# Other options: :mean, :q75, :q90.
country_statistic = :median

selected_profile = fish_profile()

# Normalisation settings: shared over all months and years in the Europe bbox.
normalisation_lower_q = 0.02
normalisation_upper_q = 0.98
normalisation_sample_stride = 12

# Output.
output_dir = joinpath(@__DIR__, "..", "output", "dynqual_K_A_country_trails_animation")
mkpath(output_dir)
animation_filename = joinpath(output_dir, "europe_K_A_country_trails_1980_2019.mp4")
animation_framerate = 12
make_animation = true
make_country_csv = true
make_static_country_plot = true

# Map style.
show_coastlines = true
show_country_borders = true
show_map_border = true
show_geo_grid = false
coastline_color = :black
coastline_linewidth = 0.55
country_border_color = (:black, 0.50)
country_border_linewidth = 0.30
map_border_color = :black
map_border_linewidth = 1.0
geo_source = "+proj=longlat +datum=WGS84"
geo_dest = "+proj=longlat +datum=WGS84"

# Color ranges. K is inferred robustly unless fixed here.
K_colormap = :inferno
A_colormap = :viridis
K_colorrange = nothing
A_colorrange = (0.0, 1.0)

# Country trail colors.
country_colors = Dict(
    "Netherlands" => :dodgerblue3,
    "Germany" => :black,
    "France" => :royalblue4,
    "Spain" => :darkorange2,
    "Italy" => :forestgreen,
    "Poland" => :crimson,
    "United Kingdom" => :purple4,
    "Sweden" => :deepskyblue4
)

# ==============================================================================
# NETCDF HELPERS
# ==============================================================================

finite_values(x) = filter(isfinite, vec(collect(skipmissing(x))))
time_index_from_year_month(year, month; base_year=1980) = (year - base_year) * 12 + month
month_name(m) = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][m]

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
    if country_statistic == :mean
        return mean(vals)
    elseif country_statistic == :q75
        return quantile(vals, 0.75)
    elseif country_statistic == :q90
        return quantile(vals, 0.90)
    else
        return median(vals)
    end
end

# ==============================================================================
# MAP AND TRAIL HELPERS
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

function add_map_panel!(fig, cell, data_obs, lon_s, lat_s; title, colormap, colorrange, label)
    br = bbox_ranges(bbox)
    ax = GeoAxis(
        fig[cell...];
        source=geo_source,
        dest=geo_dest,
        limits=(br.lon_min, br.lon_max, br.lat_min, br.lat_max),
        title=title,
        xlabel="Longitude",
        ylabel="Latitude",
        #backgroundcolor=:white,
        xgridvisible=show_geo_grid,
        ygridvisible=show_geo_grid
    )
    hm = heatmap!(ax, lon_s, lat_s, data_obs; colormap=colormap, colorrange=colorrange, nan_color=:transparent)
    plot_country_borders!(ax)
    show_coastlines && lines!(ax, GeoMakie.coastlines(); color=coastline_color, linewidth=coastline_linewidth, overdraw=true)
    plot_map_border!(ax, br)
    xlims!(ax, br.lon_min, br.lon_max)
    ylims!(ax, br.lat_min, br.lat_max)
    Colorbar(fig[cell[1], cell[2] + 1], hm; label=label, width=14, ticklabelsize=10, labelsize=11)
    return ax, hm
end

function compute_K_A(year, month, refs)
    _, res, lon, lat = run_pipeline_for_time(year, month, refs)
    K = res.amplification
    A = res.A
    K_s, lon_s, lat_s = sort_for_plot(K, lon, lat)
    A_s, _, _ = sort_for_plot(A, lon, lat)
    return K_s, A_s, lon_s, lat_s
end

function sample_K_colorrange(refs, lon, lat)
    K_colorrange !== nothing && return K_colorrange
    println("\n--- Sampling K for animation color range ---")
    kvals = Float64[]
    for (idx, (yr, mo)) in enumerate(all_times)
        idx % 24 == 1 && println("Sampling K color range: $yr-$(lpad(mo, 2, "0")) [$idx / $(length(all_times))]")
        _, res, _, _ = run_pipeline_for_time(yr, mo, refs)
        K_s, _, _ = sort_for_plot(res.amplification, lon, lat)
        sample_values!(kvals, K_s; stride=normalisation_sample_stride)
    end
    return isempty(kvals) ? (1.0, 2.0) : (quantile(kvals, 0.02), quantile(kvals, 0.98))
end

function trail_arrays(country_values, country, current_year, current_month)
    M = country_values[country]
    year_idx = current_year - start_year + 1

    # Background: all completed years before current_year, separated by NaNs.
    xb = Float64[]
    yb = Float64[]
    if year_idx > 1
        for yi in 1:(year_idx - 1)
            append!(xb, Float64.(months)); push!(xb, NaN)
            append!(yb, M[yi, :]); push!(yb, NaN)
        end
    end

    # Current year: Jan to current_month.
    xc = Float64.(1:current_month)
    yc = M[year_idx, 1:current_month]
    xd = [Float64(current_month)]
    yd = [M[year_idx, current_month]]
    return xb, yb, xc, yc, xd, yd
end

function update_trails!(trail_obs, country_values, current_year, current_month)
    for country in selected_countries
        haskey(country_values, country) || continue
        xb, yb, xc, yc, xd, yd = trail_arrays(country_values, country, current_year, current_month)
        trail_obs[country][:background_x][] = xb
        trail_obs[country][:background_y][] = yb
        trail_obs[country][:current_x][] = xc
        trail_obs[country][:current_y][] = yc
        trail_obs[country][:dot_x][] = xd
        trail_obs[country][:dot_y][] = yd
    end
end

# ==============================================================================
# STATIC OUTPUTS
# ==============================================================================

function plot_country_timeseries(country_values)
    fig = Figure(size=(1200, 650))
    ax = Axis(fig[1, 1]; title="Country $(String(country_statistic)) K over time", xlabel="Year", ylabel="K")
    for country in selected_countries
        haskey(country_values, country) || continue
        M = country_values[country]
        xs = Float64[]
        ys = Float64[]
        for yi in eachindex(years), m in months
            push!(xs, years[yi] + (m - 1) / 12)
            push!(ys, M[yi, m])
        end
        lines!(ax, xs, ys; color=country_colors[country], label=country)
    end
    axislegend(ax; position=:rt)
    save(joinpath(output_dir, "country_K_timeseries_$(String(country_statistic)).png"), fig, px_per_unit=2)
end

function plot_country_seasonality(country_values)
    fig = Figure(size=(1000, 650))
    ax = Axis(fig[1, 1]; title="Monthly country K seasonality, $(start_year)-$(end_year)", xlabel="Month", ylabel="K", xticks=(1:12, month_name.(1:12)))
    for country in selected_countries
        haskey(country_values, country) || continue
        M = country_values[country]
        ys = [mean(filter(isfinite, M[:, m])) for m in months]
        lines!(ax, months, ys; color=country_colors[country], label=country)
        scatter!(ax, months, ys; color=country_colors[country], markersize=6)
    end
    axislegend(ax; position=:rt)
    save(joinpath(output_dir, "country_K_seasonality_$(String(country_statistic)).png"), fig, px_per_unit=2)
end

function save_country_csv(country_values)
    make_country_csv || return nothing
    try
        @eval using DataFrames
        @eval using CSV
        rows = NamedTuple[]
        for country in selected_countries
            haskey(country_values, country) || continue
            M = country_values[country]
            for yi in eachindex(years), m in months
                push!(rows, (country=country, year=years[yi], month=m, K=M[yi, m], statistic=String(country_statistic)))
            end
        end
        CSV.write(joinpath(output_dir, "country_K_values_$(String(country_statistic)).csv"), DataFrame(rows))
    catch err
        @info "CSV/DataFrames not available; skipping CSV output" exception=(err, catch_backtrace())
    end
end

# ==============================================================================
# MAIN
# ==============================================================================

println("\n--- DynQual K/A maps + country seasonal trails animation ---")
println("Countries: ", selected_countries)
println("Country statistic: $country_statistic")
println("bbox: $bbox")
println("Output: $output_dir")

refs, lon, lat = compute_normalisation_refs()
country_masks, lon_s, lat_s = build_country_masks(lon, lat)
K_range = sample_K_colorrange(refs, lon, lat)
println("K color range: ", K_range)

# Precompute country values and initial map panels in one pass.
country_values = Dict(country => fill(NaN, length(years), 12) for country in selected_countries)
K_initial = nothing
A_initial = nothing

println("\n--- Precomputing monthly country K values ---")
for (frame, (yr, mo)) in enumerate(all_times)
    frame % 24 == 1 && println("Country values: $yr-$(lpad(mo, 2, "0")) [$frame / $(length(all_times))]")
    K_s, A_s, _, _ = compute_K_A(yr, mo, refs)
    if K_initial === nothing
        K_initial = K_s
        A_initial = A_s
    end
    yi = yr - start_year + 1
    for country in selected_countries
        haskey(country_masks, country) || continue
        country_values[country][yi, mo] = country_stat(K_s, country_masks[country])
    end
end

make_static_country_plot && plot_country_timeseries(country_values)
make_static_country_plot && plot_country_seasonality(country_values)
using CSV
save_country_csv(country_values)

if make_animation
    Kobs = Observable(K_initial)
    Aobs = Observable(A_initial)
    title_obs = Observable("Europe K and adaptive margin: $(start_year)-01")

    fig = Figure(size=(1500, 950))
    Label(fig[0, 1:2], title_obs; fontsize=24, tellwidth=false)

    # Two large map panels.
    add_map_panel!(fig, (1, 1), Kobs, lon_s, lat_s; title="K: amplification", colormap=K_colormap, colorrange=K_range, label="K")
    add_map_panel!(fig, (1, 3), Aobs, lon_s, lat_s; title="Adaptive margin A", colormap=A_colormap, colorrange=A_colorrange, label="A")

    # Country seasonal trail panel.
    ax_trail = Axis(
        fig[2, 1:4];
        title="Country seasonal K trajectories: previous years faint, current year highlighted",
        xlabel="Month within year",
        ylabel="Country K",
        xticks=(1:12, month_name.(1:12))
    )

    # Robust y-range from country values.
    all_country_vals = Float64[]
    for country in selected_countries
        haskey(country_values, country) || continue
        append!(all_country_vals, finite_values(country_values[country]))
    end
    if !isempty(all_country_vals)
        ymin = quantile(all_country_vals, 0.02)
        ymax = quantile(all_country_vals, 0.98)
        pad = 0.05 * (ymax - ymin + eps())
        ylims!(ax_trail, ymin - pad, ymax + pad)
    end
    xlims!(ax_trail, 1, 12)

    trail_obs = Dict{String, Dict{Symbol, Observable}}()
    for country in selected_countries
        haskey(country_values, country) || continue
        col = country_colors[country]
        xb, yb, xc, yc, xd, yd = trail_arrays(country_values, country, start_year, 1)
        trail_obs[country] = Dict(
            :background_x => Observable(xb),
            :background_y => Observable(yb),
            :current_x => Observable(xc),
            :current_y => Observable(yc),
            :dot_x => Observable(xd),
            :dot_y => Observable(yd)
        )
        lines!(ax_trail, trail_obs[country][:background_x], trail_obs[country][:background_y]; color=(col, 0.18), linewidth=1.1)
        lines!(ax_trail, trail_obs[country][:current_x], trail_obs[country][:current_y]; color=col, linewidth=3.0, label=country)
        scatter!(ax_trail, trail_obs[country][:dot_x], trail_obs[country][:dot_y]; color=col, strokecolor=:white, strokewidth=1, markersize=13)
    end
    axislegend(ax_trail; position=:rt, nbanks=2)

    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 12)

    println("\n--- Recording animation: $animation_filename ---")
    record(fig, animation_filename, 1:length(all_times); framerate=animation_framerate) do frame
        yr, mo = all_times[frame]
        println("Animation frame $frame / $(length(all_times)): $yr-$(lpad(mo, 2, "0"))")
        K_s, A_s, _, _ = compute_K_A(yr, mo, refs)
        Kobs[] = K_s
        Aobs[] = A_s
        update_trails!(trail_obs, country_values, yr, mo)
        title_obs[] = "Europe K and adaptive margin: $(yr)-$(lpad(mo, 2, "0"))"
    end
end

println("\nDone.")
println("Outputs written to: $output_dir")
println("Animation: ", make_animation ? animation_filename : "not generated")
