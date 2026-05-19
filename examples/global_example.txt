using NCDatasets
using Statistics
using CairoMakie

# Change this if your module is capitalized.
using TwoTimescaleResilience

############################################################
# User settings
############################################################

pathogen_file = raw"C:\Users\peete074\Downloads\pathogen_monthlyAvg_1980_2019.nc"
organic_file  = raw"C:\Users\peete074\Downloads\organic_monthlyAvg_1980_2019.nc"

pathogen_var = "pathogen"
organic_var  = "organic"

# Time index:
# January 1980 = 1
# December 2019 = 480
year = 2010
month = 7
time_index = (year - 1980) * 12 + month

# Start with Europe. Set bbox = nothing for global.
# Format: (lon_min, lon_max, lat_min, lat_max)
bbox = (-12.0, 35.0, 34.0, 72.0)

output_dir = joinpath(@__DIR__, "..", "output", "nc_real_raster_demo")
mkpath(output_dir)

############################################################
# Helpers
############################################################

function finite_values(x)
    return vec(x[isfinite.(x)])
end

function robust_normalise(x; lower_q=0.02, upper_q=0.98, log_transform=false)
    y = copy(x)

    if log_transform
        y .= log10.(1 .+ max.(y, 0.0))
    end

    vals = finite_values(y)
    lo = quantile(vals, lower_q)
    hi = quantile(vals, upper_q)

    z = similar(y, Float64)

    for i in eachindex(y)
        if !isfinite(y[i])
            z[i] = NaN
        else
            z[i] = clamp((y[i] - lo) / (hi - lo), 0.0, 1.0)
        end
    end

    return z
end

function replace_missing_with_nan(A)
    out = Array{Float64}(undef, size(A))
    for i in eachindex(A)
        out[i] = ismissing(A[i]) ? NaN : Float64(A[i])
    end
    return out
end

function load_nc_slice(path, varname; time_index::Int, bbox=nothing)
    ds = NCDataset(path, "r")

    lon = collect(ds["longitude"][:])
    lat = collect(ds["latitude"][:])

    if bbox === nothing
        lon_idx = eachindex(lon)
        lat_idx = eachindex(lat)
    else
        lon_min, lon_max, lat_min, lat_max = bbox
        lon_idx = findall(x -> lon_min <= x <= lon_max, lon)
        lat_idx = findall(y -> lat_min <= y <= lat_max, lat)
    end

    raw = ds[varname][lon_idx, lat_idx, time_index]

    layer = replace_missing_with_nan(raw)

    lon_sel = lon[lon_idx]
    lat_sel = lat[lat_idx]

    close(ds)

    return layer, lon_sel, lat_sel
end

function sort_for_plot(layer, lon, lat)
    L = layer
    x = lon
    y = lat

    if x[1] > x[end]
        x = reverse(x)
        L = reverse(L, dims=1)
    end

    if y[1] > y[end]
        y = reverse(y)
        L = reverse(L, dims=2)
    end

    return L, x, y
end

function plot_geo_grid(layer, lon, lat; title, filename, colormap=:viridis)
    L, x, y = sort_for_plot(layer, lon, lat)

    fig = Figure(size=(1100, 750), backgroundcolor=:white)
    ax = Axis(
        fig[1, 1],
        title=title,
        xlabel="longitude",
        ylabel="latitude",
        backgroundcolor=:white
    )

    hm = heatmap!(ax, x, y, L, colormap=colormap)
    Colorbar(fig[1, 2], hm)

    save(filename, fig, px_per_unit=2)
    return filename
end

# ESRI ASCII expects rows × columns = latitude × longitude.
# NetCDF layer is longitude × latitude.
function grid_for_ascii(layer_lonlat, lat)
    G = permutedims(layer_lonlat)  # latitude × longitude

    # ESRI ASCII convention: first data row is north.
    # If latitude is south-to-north, reverse rows.
    if lat[1] < lat[end]
        G = reverse(G, dims=1)
    end

    return G
end

function ascii_georef(lon, lat)
    dx = abs(lon[2] - lon[1])
    dy = abs(lat[2] - lat[1])

    @assert isapprox(dx, dy; rtol=1e-6) "longitude/latitude resolution differs; ASCII grid assumes square cells"

    xllcorner = minimum(lon) - dx / 2
    yllcorner = minimum(lat) - dy / 2
    cellsize = dx

    return xllcorner, yllcorner, cellsize
end

############################################################
# Load real NetCDF rasters
############################################################

println("Loading pathogen raster...")
pathogen_raw, lon, lat = load_nc_slice(
    pathogen_file,
    pathogen_var;
    time_index=time_index,
    bbox=bbox
)

println("Loading organic raster...")
organic_raw, lon2, lat2 = load_nc_slice(
    organic_file,
    organic_var;
    time_index=time_index,
    bbox=bbox
)

@assert lon == lon2
@assert lat == lat2
@assert size(pathogen_raw) == size(organic_raw)

println("Loaded raster size: ", size(pathogen_raw))
println("Longitude range: ", minimum(lon), " to ", maximum(lon))
println("Latitude range:  ", minimum(lat), " to ", maximum(lat))

############################################################
# Normalise to dimensionless stressor layers
############################################################

# Pathogen has large dynamic range, so log transform is usually helpful.
pathogen_norm = robust_normalise(
    pathogen_raw;
    lower_q=0.02,
    upper_q=0.98,
    log_transform=true
)

# Organic pollution/BOD can also be skewed, but start without log.
organic_norm = robust_normalise(
    organic_raw;
    lower_q=0.02,
    upper_q=0.98,
    log_transform=false
)

############################################################
# Combine into background burden
############################################################

layers = [pathogen_norm, organic_norm]

# First defensible simple weighting.
weights = [0.50, 0.50]

# More-than-additive background interaction:
# pathogen × organic pollution.
interaction = zeros(2, 2)
interaction[1, 2] = 0.25

Bgrid = compute_background_index_grid(layers, weights; interaction=interaction)

############################################################
# Background burden -> margin -> restoring force -> amplification
############################################################

params = BackgroundParams(
    A0 = 1.0,
    alpha = 1.15,
    KB = 0.55,
    hill = 2.0,
    use_saturating_phi = true,
    lambda_min = 0.04,
    lambda_max = 1.0,
    KA = 0.30
)

Agrid      = adaptive_margin_grid(Bgrid, params)
lambdagrid = restoring_force_grid(Bgrid, params)
Fgrid      = amplification_factor_grid(Bgrid, params)

############################################################
# Write ASCII grids
############################################################

xllcorner, yllcorner, cellsize = ascii_georef(lon, lat)

write_ascii_grid(
    joinpath(output_dir, "pathogen_normalised.asc"),
    grid_for_ascii(pathogen_norm, lat);
    xllcorner=xllcorner,
    yllcorner=yllcorner,
    cellsize=cellsize
)

write_ascii_grid(
    joinpath(output_dir, "organic_normalised.asc"),
    grid_for_ascii(organic_norm, lat);
    xllcorner=xllcorner,
    yllcorner=yllcorner,
    cellsize=cellsize
)

write_ascii_grid(
    joinpath(output_dir, "background_index_B.asc"),
    grid_for_ascii(Bgrid, lat);
    xllcorner=xllcorner,
    yllcorner=yllcorner,
    cellsize=cellsize
)

write_ascii_grid(
    joinpath(output_dir, "adaptive_margin_A.asc"),
    grid_for_ascii(Agrid, lat);
    xllcorner=xllcorner,
    yllcorner=yllcorner,
    cellsize=cellsize
)

write_ascii_grid(
    joinpath(output_dir, "restoring_force_lambda.asc"),
    grid_for_ascii(lambdagrid, lat);
    xllcorner=xllcorner,
    yllcorner=yllcorner,
    cellsize=cellsize
)

write_ascii_grid(
    joinpath(output_dir, "amplification_factor.asc"),
    grid_for_ascii(Fgrid, lat);
    xllcorner=xllcorner,
    yllcorner=yllcorner,
    cellsize=cellsize
)

############################################################
# Plot real rasters
############################################################

plot_geo_grid(
    pathogen_norm,
    lon,
    lat;
    title="Normalised pathogen stress, $(year)-$(month)",
    filename=joinpath(output_dir, "pathogen_normalised.png")
)

plot_geo_grid(
    organic_norm,
    lon,
    lat;
    title="Normalised organic stress, $(year)-$(month)",
    filename=joinpath(output_dir, "organic_normalised.png")
)

plot_geo_grid(
    Bgrid,
    lon,
    lat;
    title="Combined background burden B, $(year)-$(month)",
    filename=joinpath(output_dir, "background_index_B.png")
)

plot_geo_grid(
    Agrid,
    lon,
    lat;
    title="Adaptive margin A(B), $(year)-$(month)",
    filename=joinpath(output_dir, "adaptive_margin_A.png"),
    colormap=:balance
)

plot_geo_grid(
    lambdagrid,
    lon,
    lat;
    title="Restoring force λ(B), $(year)-$(month)",
    filename=joinpath(output_dir, "restoring_force_lambda.png")
)

plot_geo_grid(
    Fgrid,
    lon,
    lat;
    title="Amplification factor λ(0)/λ(B), $(year)-$(month)",
    filename=joinpath(output_dir, "amplification_factor.png"),
    colormap=:magma
)

############################################################
# Print summaries
############################################################

println("\n================ Raster summary ================")
println("year:       ", year)
println("month:      ", month)
println("time index: ", time_index)

for (name, grid) in [
    ("pathogen_norm", pathogen_norm),
    ("organic_norm", organic_norm),
    ("B", Bgrid),
    ("A", Agrid),
    ("lambda", lambdagrid),
    ("amplification", Fgrid)
]
    vals = finite_values(grid)

    println("\n", name)
    println("  min:    ", minimum(vals))
    println("  mean:   ", mean(vals))
    println("  median: ", median(vals))
    println("  max:    ", maximum(vals))
end

println("\nOutputs written to:")
println(output_dir)