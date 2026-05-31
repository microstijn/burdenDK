# dynqual_synthetic_isimip_pressure_demo.jl

using TwoTimescaleResilience
using Statistics
using Dates
using JSON
using CSV
using DataFrames

# 1. Unconditional dependency for memory-conscious NetCDF reading
using NCDatasets

function get_env_or_fallback(env_var::String, fallback::String)
    val = get(ENV, env_var, fallback)
    return isempty(val) ? fallback : val
end

# Memory-conscious local helpers for reading NetCDF files
function open_dynqual_dataset(file_path::String)
    return NCDataset(file_path, "r")
end

function detect_lon_lat_time_vars(ds)
    lon_var = haskey(ds, "lon") ? "lon" : (haskey(ds, "longitude") ? "longitude" : "")
    lat_var = haskey(ds, "lat") ? "lat" : (haskey(ds, "latitude") ? "latitude" : "")
    time_var = haskey(ds, "time") ? "time" : ""
    
    if lon_var == "" || lat_var == ""
        error("Could not identify lon/lat variables in dataset")
    end
    return lon_var, lat_var, time_var
end

function detect_main_variable(ds, preferred_var_name::String)
    if haskey(ds, preferred_var_name)
        return preferred_var_name
    end
    lon_var, lat_var, time_var = detect_lon_lat_time_vars(ds)
    for v in keys(ds)
        if v ∉ [lon_var, lat_var, time_var, "crs", "time_bnds", "lat_bnds", "lon_bnds"]
            return v
        end
    end
    error("Could not identify main variable; preferred '$preferred_var_name' not found.")
end

function subset_indices(lons, lats, bbox, spatial_stride::Int)
    lon_idx = findall(x -> bbox.lon[1] <= x <= bbox.lon[2], lons)
    lat_idx = findall(x -> bbox.lat[1] <= x <= bbox.lat[2], lats)
    
    if isempty(lon_idx) || isempty(lat_idx)
        error("Bbox subset resulted in empty grid.")
    end
    
    if spatial_stride > 1
        lon_idx = lon_idx[1:spatial_stride:end]
        lat_idx = lat_idx[1:spatial_stride:end]
    end
    
    return lon_idx, lat_idx
end

function read_month_slice(ds, var_name::String, lon_idx::Vector{Int}, lat_idx::Vector{Int}, t::Int)
    dim_names = dimnames(ds[var_name])
    
    # We must construct the slicing tuple dynamically based on dimension names
    # typically ("lon", "lat", "time") or ("longitude", "latitude", "time")
    
    idx_tuple = Any[Colon() for _ in 1:length(dim_names)]
    for (i, dname) in enumerate(dim_names)
        if dname == "lon" || dname == "longitude"
            idx_tuple[i] = lon_idx
        elseif dname == "lat" || dname == "latitude"
            idx_tuple[i] = lat_idx
        elseif dname == "time"
            idx_tuple[i] = t
        end
    end
    
    data = ds[var_name][idx_tuple...]
    
    # Force data into 2D [lon, lat] explicitly to handle different dimension orders safely
    lon_dim_idx = findfirst(x -> x == "lon" || x == "longitude", dim_names)
    lat_dim_idx = findfirst(x -> x == "lat" || x == "latitude", dim_names)
    
    # If the file had [time, lat, lon], data is a Matrix but indices are swapped.
    # We always return data_f as [lon_idx, lat_idx]
    
    data_f = Array{Float32, 2}(undef, length(lon_idx), length(lat_idx))
    
    for i in 1:length(lon_idx)
        for j in 1:length(lat_idx)
            v = if lon_dim_idx == 1 && lat_dim_idx == 2
                data[i, j]
            elseif lon_dim_idx == 2 && lat_dim_idx == 1
                data[j, i]
            elseif length(dim_names) == 3 && lon_dim_idx == 3 && lat_dim_idx == 2
                # e.g., time, lat, lon
                data[j, i]
            elseif length(dim_names) == 3 && lon_dim_idx == 2 && lat_dim_idx == 3
                # e.g., time, lon, lat
                data[i, j]
            else
                # Fallback assuming data was fetched in order [lon, lat]
                data[i, j]
            end
            
            if ismissing(v) || isnothing(v) || typeof(v) <: Dates.TimeType || (v isa Number && !isfinite(v))
                data_f[i, j] = NaN32
            else
                data_f[i, j] = Float32(v)
            end
        end
    end
    
    return data_f
end

function estimate_log_quantiles_sampled(
    file_path::String,
    preferred_var_name::String,
    lon_idx::Vector{Int},
    lat_idx::Vector{Int};
    time_stride::Int = 6,
    max_samples::Int = 2_000_000
)
    ds = open_dynqual_dataset(file_path)
    var_name = detect_main_variable(ds, preferred_var_name)
    _, _, time_var = detect_lon_lat_time_vars(ds)
    n_time = length(ds[time_var])
    
    sampled_vals = Float64[]
    n_total_seen = 0
    n_finite_seen = 0
    
    time_indices = 1:time_stride:n_time
    # Adaptive spatial striding based on expected samples
    cells_per_slice = length(lon_idx) * length(lat_idx)
    expected_finite = cells_per_slice * length(time_indices) * 0.9 # assume 90% finite
    
    spatial_stride_inner = 1
    if expected_finite > max_samples
        spatial_stride_inner = ceil(Int, expected_finite / max_samples)
    end
    
    for t in time_indices
        slice = read_month_slice(ds, var_name, lon_idx, lat_idx, t)
        
        # Subsample spatially
        if spatial_stride_inner > 1
            slice_1d = slice[1:spatial_stride_inner:end]
        else
            slice_1d = slice
        end
        
        n_total_seen += length(slice)
        
        for v in slice_1d
            if !isnan(v)
                n_finite_seen += spatial_stride_inner # Approximate count
                v_pos = max(0.0, Float64(v))
                log_val = log1p(v_pos)
                if length(sampled_vals) < max_samples
                    push!(sampled_vals, log_val)
                end
            end
        end
    end
    close(ds)
    
    fraction_missing_estimate = 1.0 - (n_finite_seen / max(1, n_total_seen))
    
    if isempty(sampled_vals)
        error("No valid values found for $preferred_var_name during sampling.")
    end
    
    p02 = quantile(sampled_vals, 0.02)
    p98 = quantile(sampled_vals, 0.98)
    
    if p98 <= p02
        @warn "Robust scaling degenerate for $preferred_var_name: p02=$p02, p98=$p98"
        p98 = p02 + 1.0
    end
    
    return p02, p98, n_total_seen, n_finite_seen, length(sampled_vals), fraction_missing_estimate
end

# Helper to scale a slice robustly
function robust_log_scale_slice(slice::Matrix{Float32}, p02::Float64, p98::Float64)
    scaled_slice = similar(slice)
    for i in eachindex(slice)
        v = slice[i]
        if isnan(v)
            scaled_slice[i] = 0.0f0
        else
            v_pos = max(0.0f0, v)
            log_val = log1p(v_pos)
            s = clamp((log_val - p02) / (p98 - p02), 0.0, 1.0)
            scaled_slice[i] = Float32(s)
        end
    end
    return scaled_slice
end

# Plotting helpers
function make_dynqual_raw_climatology_figure(fig_dir, lons, lats, raw_clim_maps)
    println("  -> Figure 1: Raw climatology")
    fig1 = Figure(size = (1200, 800))
    Label(fig1[0, :], "Real DynQual water-quality fields over Europe, 1980–2019 climatology", font=:bold, fontsize=20)
    
    for (i, p_info) in enumerate(raw_clim_maps)
        row = cld(i, 2)
        col = mod1(i, 2) * 2 - 1
        ax = Axis(fig1[row, col], title=p_info.title, aspect=DataAspect())
        hm = heatmap!(ax, lons, lats, p_info.data, colormap=:viridis)
        Colorbar(fig1[row, col+1], hm, label="log1p")
    end
    save(joinpath(fig_dir, "dynqual_raw_climatology_maps.png"), fig1)
end

function make_dynqual_derived_pressure_layers_figure(fig_dir, lons, lats, derived_clim_maps)
    println("  -> Figure 2: Derived pressure archetypes")
    fig2 = Figure(size = (1500, 1000))
    Label(fig2[0, :], "DynQual-derived synthetic pressure archetypes (Climatological Mean)", font=:bold, fontsize=20)
    
    for (i, p_info) in enumerate(derived_clim_maps)
        row = cld(i, 3)
        col = mod1(i, 3) * 2 - 1
        ax = Axis(fig2[row, col], title=p_info.title, aspect=DataAspect())
        hm = heatmap!(ax, lons, lats, p_info.data, colormap=:plasma, colorrange=(0.0, 1.0))
        Colorbar(fig2[row, col+1], hm, label="Relative Impairment")
    end
    save(joinpath(fig_dir, "dynqual_derived_pressure_layers.png"), fig2)
end

function make_dynqual_vulnerability_regime_maps_figure(fig_dir, lons, lats, nx, ny, k_clusters, map_base, map_rec)
    println("  -> Figure 3: Vulnerability regimes")
    fig3 = Figure(size = (1200, 600))
    Label(fig3[0, :], "Threshold-free vulnerability regimes from DynQual-derived pressures", font=:bold, fontsize=20)
    
    ax_base = Axis(fig3[1, 1], title="Baseline Tranche (1980s) Vulnerability Regimes", aspect=DataAspect())
    map_base_2d = reshape(map_base, nx, ny)
    hm_base = heatmap!(ax_base, lons, lats, map_base_2d, colormap=:Set1_5, colorrange=(0.5, k_clusters + 0.5))
    Colorbar(fig3[1, 2], hm_base, ticks=1:k_clusters)
    
    ax_rec = Axis(fig3[1, 3], title="Recent Tranche (2010s) Vulnerability Regimes", aspect=DataAspect())
    map_rec_2d = reshape(map_rec, nx, ny)
    hm_rec = heatmap!(ax_rec, lons, lats, map_rec_2d, colormap=:Set1_5, colorrange=(0.5, k_clusters + 0.5))
    Colorbar(fig3[1, 4], hm_rec, ticks=1:k_clusters)
    
    save(joinpath(fig_dir, "dynqual_vulnerability_regime_maps.png"), fig3)
end

function make_dynqual_regime_explanation_heatmap_figure(fig_dir, k_clusters, baseline_kept_names, baseline_centroids)
    println("  -> Figure 4: Regime explanation heatmap")
    fig4 = Figure(size = (1000, 800))
    Label(fig4[0, :], "What distinguishes the vulnerability regimes?", font=:bold, fontsize=20)
    
    ax_heat = Axis(fig4[1, 1],
        xticks = (1:k_clusters, ["Cluster $i" for i in 1:k_clusters]),
        yticks = (1:length(baseline_kept_names), baseline_kept_names),
        xticklabelrotation = pi/4
    )
    
    hm_heat = heatmap!(ax_heat, 1:k_clusters, 1:length(baseline_kept_names), baseline_centroids', colormap=:RdBu, colorrange=(-3.0, 3.0))
    Colorbar(fig4[1, 2], hm_heat, label="Standardised Feature Value")
    
    for c in 1:k_clusters
        for f in 1:length(baseline_kept_names)
            val = baseline_centroids[c, f]
            text!(ax_heat, string(round(val, digits=2)), position=(c, f), 
                  align=(:center, :center), color=abs(val) > 1.5 ? :white : :black, fontsize=12)
        end
    end
    
    save(joinpath(fig_dir, "dynqual_regime_explanation_heatmap.png"), fig4)
end

function run_dynqual_synthetic_isimip_pressure_demo(;
    output_dir = joinpath(dirname(@__DIR__), "output", "dynqual_synthetic_isimip_pressure_demo"),
    bbox = (lon = (-25.0, 45.0), lat = (34.0, 72.0))
)
    println("--- Starting Memory-Conscious DynQual Synthetic ISIMIP Pressure Demo ---")
    mkpath(output_dir)

    # 1. Paths and environment variables
    default_bod_file = raw"C:\Users\peete074\Downloads\organic_monthlyAvg_1980_2019.nc"
    default_fc_file = raw"C:\Users\peete074\Downloads\pathogen_monthlyAvg_1980_2019.nc"
    default_tds_file = raw"C:\Users\peete074\Downloads\TDSload_monthlyAvg_1980_2019.nc"
    default_bodload_file = raw"C:\Users\peete074\Downloads\BODload_monthlyAvg_1980_2019.nc"

    bod_file = get_env_or_fallback("TTR_DYNQUAL_BOD_FILE", default_bod_file)
    fc_file = get_env_or_fallback("TTR_DYNQUAL_FC_FILE", default_fc_file)
    tds_file = get_env_or_fallback("TTR_DYNQUAL_TDS_FILE", default_tds_file)
    bodload_file = get_env_or_fallback("TTR_DYNQUAL_BODLOAD_FILE", default_bodload_file)

    for (name, path) in [("BOD", bod_file), ("FC", fc_file), ("TDS", tds_file), ("BODload", bodload_file)]
        if !isfile(path)
            error("Required file for $name not found at path: $path. Please set the appropriate TTR_DYNQUAL_*_FILE environment variable.")
        end
    end

    lon_min = parse(Float64, get_env_or_fallback("TTR_DYNQUAL_LON_MIN", string(bbox.lon[1])))
    lon_max = parse(Float64, get_env_or_fallback("TTR_DYNQUAL_LON_MAX", string(bbox.lon[2])))
    lat_min = parse(Float64, get_env_or_fallback("TTR_DYNQUAL_LAT_MIN", string(bbox.lat[1])))
    lat_max = parse(Float64, get_env_or_fallback("TTR_DYNQUAL_LAT_MAX", string(bbox.lat[2])))
    bbox_override = (lon = (lon_min, lon_max), lat = (lat_min, lat_max))

    n_species = parse(Int, get_env_or_fallback("TTR_DYNQUAL_N_SPECIES", "8"))
    k_clusters = parse(Int, get_env_or_fallback("TTR_DYNQUAL_CLUSTER_K", "5"))
    make_plots = get_env_or_fallback("TTR_DYNQUAL_MAKE_PLOTS", "true") == "true"
    write_netcdf = get_env_or_fallback("TTR_DYNQUAL_WRITE_NETCDF", "false") == "true"

    spatial_stride = parse(Int, get_env_or_fallback("TTR_DYNQUAL_SPATIAL_STRIDE", "1"))
    time_stride = parse(Int, get_env_or_fallback("TTR_DYNQUAL_QUANTILE_TIME_STRIDE", "6"))
    max_samples = parse(Int, get_env_or_fallback("TTR_DYNQUAL_MAX_QUANTILE_SAMPLES", "2000000"))

    if make_plots
        try
            @eval using CairoMakie
        catch
            error("CairoMakie is required for plotting but is not available. Install it or set TTR_DYNQUAL_MAKE_PLOTS=false.")
        end
    end

    println("Output Directory: ", output_dir)
    println("BBox: ", bbox_override)
    println("Spatial Stride: ", spatial_stride)
    println("Plotting Enabled: ", make_plots)

    # 2. Get subset indices
    ds_bod = open_dynqual_dataset(bod_file)
    lon_var, lat_var, time_var = detect_lon_lat_time_vars(ds_bod)
    lons_full = ds_bod[lon_var][:]
    lats_full = ds_bod[lat_var][:]
    n_time = length(ds_bod[time_var])
    close(ds_bod)

    if n_time != 480
        @warn "Expected 480 months (1980-2019) but got $n_time."
    end

    lon_idx, lat_idx = subset_indices(lons_full, lats_full, bbox_override, spatial_stride)
    sub_lons = lons_full[lon_idx]
    sub_lats = lats_full[lat_idx]
    
    nx = length(lon_idx)
    ny = length(lat_idx)
    println("Grid size after subset and stride: nx=$nx, ny=$ny, time=$n_time")

    # 3. Robust scaling sampled
    println("Estimating robust log-scaling quantiles...")
    scaling_summary = NamedTuple[]
    scaling_params = Dict{String, Tuple{Float64, Float64}}()

    for (fname, path, preferred) in [
        ("BOD", bod_file, "organic_monthlyAvg_1980_2019"),
        ("FC", fc_file, "pathogen_monthlyAvg_1980_2019"),
        ("TDS", tds_file, "TDSload_monthlyAvg_1980_2019"),
        ("BODload", bodload_file, "BODload_monthlyAvg_1980_2019")
    ]
        println("  -> Sampling $fname ...")
        p02, p98, nt, nf, ns, frac_miss = estimate_log_quantiles_sampled(
            path, preferred, lon_idx, lat_idx;
            time_stride=time_stride, max_samples=max_samples
        )
        scaling_params[fname] = (p02, p98)
        push!(scaling_summary, (
            variable_name = fname,
            preferred_var = preferred,
            n_total_seen = nt,
            n_finite_seen = nf,
            n_sampled = ns,
            fraction_missing_estimate = frac_miss,
            p02_log1p = p02,
            p98_log1p = p98
        ))
    end
    CSV.write(joinpath(output_dir, "dynqual_scaling_summary.csv"), scaling_summary)

    # Free up memory before streaming
    GC.gc()

    # 4. Model Chain & Pressure mappings
    pressure_mapping_weights = [
        (pressure_proxy = "organic_oxygen_demand_proxy", source_variables = "BOD", transform = "robust log 0-1", interpretation = "organic pollution / oxygen-demand stress", primary_deb_axis = "maintenance", secondary_deb_axis = "growth", memory_rho = 0.30, internal_magnification_K = 1.0, notes = "", w_A = 0.0, w_M = 0.8, w_G = 0.2, w_R = 0.0),
        (pressure_proxy = "pathogen_exposure_proxy", source_variables = "FC/pathogen", transform = "robust log 0-1", interpretation = "microbial/pathogen pressure", primary_deb_axis = "maintenance", secondary_deb_axis = "reproduction", memory_rho = 0.10, internal_magnification_K = 1.0, notes = "", w_A = 0.0, w_M = 0.6, w_G = 0.0, w_R = 0.4),
        (pressure_proxy = "ionic_salinity_proxy", source_variables = "TDSload", transform = "robust log 0-1", interpretation = "dissolved-solids / ionic / osmotic pressure proxy", primary_deb_axis = "maintenance", secondary_deb_axis = "assimilation", memory_rho = 0.70, internal_magnification_K = 1.5, notes = "Load, not local concentration.", w_A = 0.3, w_M = 0.7, w_G = 0.0, w_R = 0.0),
        (pressure_proxy = "wastewater_source_proxy", source_variables = "BODload", transform = "robust log 0-1", interpretation = "organic source/load pressure", primary_deb_axis = "growth", secondary_deb_axis = "reproduction", memory_rho = 0.40, internal_magnification_K = 1.0, notes = "", w_A = 0.0, w_M = 0.0, w_G = 0.5, w_R = 0.5),
        (pressure_proxy = "low_flow_concentration_proxy", source_variables = "BOD + BODload", transform = "derived heuristic", interpretation = "relative low-flow / dilution proxy", primary_deb_axis = "maintenance", secondary_deb_axis = "", memory_rho = 0.50, internal_magnification_K = 1.2, notes = "Not real discharge.", w_A = 0.0, w_M = 1.0, w_G = 0.0, w_R = 0.0),
        (pressure_proxy = "combined_wastewater_proxy", source_variables = "composite", transform = "weighted sum", interpretation = "composite wastewater pressure", primary_deb_axis = "maintenance", secondary_deb_axis = "", memory_rho = 0.40, internal_magnification_K = 1.0, notes = "", w_A = 0.2, w_M = 0.4, w_G = 0.2, w_R = 0.2)
    ]
    CSV.write(joinpath(output_dir, "dynqual_pressure_mapping.csv"), pressure_mapping_weights)
    n_proxies = length(pressure_mapping_weights)

    amp_lib = load_amp_species_library()
    species_keys = collect(keys(amp_lib))[1:min(n_species, length(amp_lib))]
    
    selected_species_summary = NamedTuple[]
    species_params = []
    
    for (i, sp_key) in enumerate(species_keys)
        rec = amp_lib[sp_key]
        p = amp_record_to_deb_params(rec)
        push!(species_params, p)
        push!(selected_species_summary, (
            species_key = sp_key,
            species_name = replace(sp_key, "_" => " "),
            archetype_labels = get(rec, "archetype_labels", ""),
            A0 = p.A0,
            lambda_min = p.lambda_min,
            lambda_max = p.lambda_max,
            KA = p.KA,
            alpha_assimilation = p.alpha_axes[1],
            alpha_maintenance = p.alpha_axes[2],
            alpha_growth = p.alpha_axes[3],
            alpha_reproduction = p.alpha_axes[4]
        ))
    end
    CSV.write(joinpath(output_dir, "selected_dynqual_demo_species.csv"), selected_species_summary)

    # 5. Tranches Definitions
    n_months_per_tranche = 120
    actual_n_tranches = min(4, floor(Int, n_time / n_months_per_tranche))
    
    tranches = NamedTuple[]
    for h in 1:actual_n_tranches
        start_m = (h - 1) * n_months_per_tranche + 1
        end_m = h * n_months_per_tranche
        year_start = 1980 + (h - 1) * 10
        year_end = 1980 + h * 10 - 1
        push!(tranches, (
            tranche_id = h,
            label = "tranche_$h",
            year_start = year_start,
            year_end = year_end,
            month_start = start_m,
            month_end = end_m,
            n_months = end_m - start_m + 1
        ))
    end
    CSV.write(joinpath(output_dir, "dynqual_tranche_definitions.csv"), tranches)

    # 6. Streaming Tranche Processing
    println("Streaming data and computing tranche summaries...")

    # Climatology Accumulators
    sum_raw_bod = zeros(Float32, nx, ny)
    sum_raw_fc = zeros(Float32, nx, ny)
    sum_raw_tds = zeros(Float32, nx, ny)
    sum_raw_bodload = zeros(Float32, nx, ny)

    sum_org_p = zeros(Float32, nx, ny)
    sum_pat_p = zeros(Float32, nx, ny)
    sum_ion_p = zeros(Float32, nx, ny)
    sum_was_p = zeros(Float32, nx, ny)
    sum_low_p = zeros(Float32, nx, ny)
    sum_com_p = zeros(Float32, nx, ny)
    
    count_clim = zeros(Int32, nx, ny)
    
    # Pressure Memory State
    B_state = zeros(Float32, nx, ny, n_proxies)
    
    # Baseline Clustering Objects
    baseline_kept_indices = Int[]
    baseline_kept_names = String[]
    baseline_means = Float64[]
    baseline_stds = Float64[]
    baseline_centroids = zeros(Float64, 0, 0)
    tranche_clusters = Dict{Int, Any}()

    # We need to compute an overarching p98 for low_flow pressure to scale it properly.
    # Since low flow pressure requires streaming to know, we'll use a heuristic for scaling it per tranche, 
    # or just clamp it. The instructions say "if lf_p98 > 0, scale it." We will stream and scale dynamically per tranche 
    # for simplicity, or just use it raw clamped to [0,1]
    
    for h in 1:actual_n_tranches
        tranche = tranches[h]
        start_m = tranche.month_start
        end_m = tranche.month_end
        n_m = tranche.n_months
        
        println("  -> Processing Tranche $h ($start_m to $end_m)...")
        
        # Open datasets
        ds_bod = open_dynqual_dataset(bod_file)
        var_bod = detect_main_variable(ds_bod, "organic_monthlyAvg_1980_2019")
        
        ds_fc = open_dynqual_dataset(fc_file)
        var_fc = detect_main_variable(ds_fc, "pathogen_monthlyAvg_1980_2019")
        
        ds_tds = open_dynqual_dataset(tds_file)
        var_tds = detect_main_variable(ds_tds, "TDSload_monthlyAvg_1980_2019")
        
        ds_bodload = open_dynqual_dataset(bodload_file)
        var_bodload = detect_main_variable(ds_bodload, "BODload_monthlyAvg_1980_2019")
        
        # Tranche Accumulators for Features
        sum_E_a = zeros(Float32, nx, ny)
        sum_E_m = zeros(Float32, nx, ny)
        sum_E_g = zeros(Float32, nx, ny)
        sum_E_r = zeros(Float32, nx, ny)
        
        # For Q and F exact p95 within tranche, we need n_m * n_species values per cell.
        # This is n_cells * 120 * 8 * 4 bytes = ~3.8MB, totally safe to allocate per tranche.
        Q_tranche = zeros(Float32, nx, ny, n_m, n_species)
        F_tranche = zeros(Float32, nx, ny, n_m, n_species)
        
        for (m_idx, t) in enumerate(start_m:end_m)
            raw_bod = read_month_slice(ds_bod, var_bod, lon_idx, lat_idx, t)
            raw_fc = read_month_slice(ds_fc, var_fc, lon_idx, lat_idx, t)
            raw_tds = read_month_slice(ds_tds, var_tds, lon_idx, lat_idx, t)
            raw_bodload = read_month_slice(ds_bodload, var_bodload, lon_idx, lat_idx, t)
            
            scaled_bod = robust_log_scale_slice(raw_bod, scaling_params["BOD"]...)
            scaled_fc = robust_log_scale_slice(raw_fc, scaling_params["FC"]...)
            scaled_tds = robust_log_scale_slice(raw_tds, scaling_params["TDS"]...)
            scaled_bodload = robust_log_scale_slice(raw_bodload, scaling_params["BODload"]...)
            
            # Update Raw Climatology sums
            for x in 1:nx, y in 1:ny
                # Only add if raw bod is valid (proxy for water cell)
                if !isnan(raw_bod[x,y])
                    sum_raw_bod[x,y] += log1p(max(0.0f0, raw_bod[x,y]))
                    sum_raw_fc[x,y] += log1p(max(0.0f0, raw_fc[x,y]))
                    sum_raw_tds[x,y] += log1p(max(0.0f0, raw_tds[x,y]))
                    sum_raw_bodload[x,y] += log1p(max(0.0f0, raw_bodload[x,y]))
                    count_clim[x,y] += 1
                end
            end
            
            # Low flow proxy computation
            low_flow_pressure = zeros(Float32, nx, ny)
            for x in 1:nx, y in 1:ny
                bod_val = max(0.0f0, raw_bod[x,y])
                bodload_val = max(0.0f0, raw_bodload[x,y])
                if isnan(bod_val) || isnan(bodload_val)
                    low_flow_pressure[x,y] = 0.0f0
                else
                    q_proxy = bodload_val / (bod_val + 1f-6)
                    low_flow_pressure[x,y] = clamp(scaled_bod[x,y] * exp(-0.1f0 * q_proxy), 0.0f0, 1.0f0)
                end
            end
            
            combined_pressure = zeros(Float32, nx, ny)
            for x in 1:nx, y in 1:ny
                combined_pressure[x,y] = clamp(0.4f0 * scaled_bod[x,y] + 0.3f0 * scaled_fc[x,y] + 0.2f0 * scaled_bodload[x,y] + 0.1f0 * low_flow_pressure[x,y], 0.0f0, 1.0f0)
            end
            
            # Update Derived Climatology sums
            for x in 1:nx, y in 1:ny
                if !isnan(raw_bod[x,y])
                    sum_org_p[x,y] += scaled_bod[x,y]
                    sum_pat_p[x,y] += scaled_fc[x,y]
                    sum_ion_p[x,y] += scaled_tds[x,y]
                    sum_was_p[x,y] += scaled_bodload[x,y]
                    sum_low_p[x,y] += low_flow_pressure[x,y]
                    sum_com_p[x,y] += combined_pressure[x,y]
                end
            end
            
            # Route through memory & DEB
            proxies_slices = [scaled_bod, scaled_fc, scaled_tds, scaled_bodload, low_flow_pressure, combined_pressure]
            
            for x in 1:nx, y in 1:ny
                if isnan(raw_bod[x,y]) continue end
                
                # Memory update
                for p in 1:n_proxies
                    rho = pressure_mapping_weights[p].memory_rho
                    K_mag = pressure_mapping_weights[p].internal_magnification_K
                    P_pt = proxies_slices[p][x, y]
                    B_state[x, y, p] = rho * B_state[x, y, p] + (1.0f0 - rho) * K_mag * P_pt
                end
                
                # E_axis
                E_a, E_m, E_g, E_r = 0.0f0, 0.0f0, 0.0f0, 0.0f0
                for p in 1:n_proxies
                    E_proxy = clamp(B_state[x, y, p], 0.0f0, 1.0f0)
                    E_a += pressure_mapping_weights[p].w_A * E_proxy
                    E_m += pressure_mapping_weights[p].w_M * E_proxy
                    E_g += pressure_mapping_weights[p].w_G * E_proxy
                    E_r += pressure_mapping_weights[p].w_R * E_proxy
                end
                
                E_a = clamp(E_a, 0.0f0, 1.0f0)
                E_m = clamp(E_m, 0.0f0, 1.0f0)
                E_g = clamp(E_g, 0.0f0, 1.0f0)
                E_r = clamp(E_r, 0.0f0, 1.0f0)
                
                sum_E_a[x, y] += E_a
                sum_E_m[x, y] += E_m
                sum_E_g[x, y] += E_g
                sum_E_r[x, y] += E_r
                
                axes_imp = (assimilation=Float64(E_a), maintenance=Float64(E_m), growth=Float64(E_g), reproduction=Float64(E_r))
                
                for sp in 1:n_species
                    resp = compute_adaptive_margin_response_from_impairment(axes_imp, species_params[sp]; mixture_effect_model="grouped_ca_then_ia_axis_effects")
                    Q_tranche[x, y, m_idx, sp] = Float32(resp.Q_t)
                    F_tranche[x, y, m_idx, sp] = Float32(resp.amplification)
                end
            end
        end
        
        close(ds_bod)
        close(ds_fc)
        close(ds_tds)
        close(ds_bodload)
        GC.gc()
        
        # Build compact feature matrix for this tranche
        feature_names = [
            "mean_E_assimilation_grouped", "mean_E_maintenance_grouped", "mean_E_growth_grouped", "mean_E_reproduction_grouped",
            "mean_Q_grouped", "p95_Q_grouped", "mean_F_grouped", "p95_F_grouped", "max_F_grouped", "axis_entropy"
        ]
        n_features = length(feature_names)
        
        # We only keep cells that are valid (water)
        valid_cells = CartesianIndices((1:nx, 1:ny))
        flat_valid = filter(idx -> count_clim[idx[1], idx[2]] > 0, vec(valid_cells))
        n_valid = length(flat_valid)
        
        feat_matrix = zeros(Float64, n_valid, n_features)
        
        for (i, idx) in enumerate(flat_valid)
            x, y = idx[1], idx[2]
            
            # Means
            feat_matrix[i, 1] = sum_E_a[x, y] / n_m
            feat_matrix[i, 2] = sum_E_m[x, y] / n_m
            feat_matrix[i, 3] = sum_E_g[x, y] / n_m
            feat_matrix[i, 4] = sum_E_r[x, y] / n_m
            
            # Q & F aggregations
            Q_vals = vec(Q_tranche[x, y, :, :])
            F_vals = vec(F_tranche[x, y, :, :])
            
            feat_matrix[i, 5] = mean(Q_vals)
            feat_matrix[i, 6] = quantile(Q_vals, 0.95)
            feat_matrix[i, 7] = mean(F_vals)
            feat_matrix[i, 8] = quantile(F_vals, 0.95)
            feat_matrix[i, 9] = maximum(F_vals)
            
            # Entropy
            E_total = feat_matrix[i, 1] + feat_matrix[i, 2] + feat_matrix[i, 3] + feat_matrix[i, 4]
            entropy = 0.0
            if E_total > 1e-6
                for f in 1:4
                    p = feat_matrix[i, f] / E_total
                    if p > 0
                        entropy -= p * log(p)
                    end
                end
            end
            feat_matrix[i, 10] = entropy
        end
        
        # Free memory
        Q_tranche = zeros(Float32,0)
        F_tranche = zeros(Float32,0)
        GC.gc()
        
        # Clustering
        if h == 1
            standardized = standardize_threshold_free_vulnerability_features(
                feat_matrix,
                feature_names
            )
            
            baseline_kept_indices = standardized.kept_feature_indices
            baseline_kept_names = standardized.standardized_feature_names
            baseline_means = standardized.means
            baseline_stds = standardized.stds
            
            cluster_res = cluster_threshold_free_vulnerability_regimes(
                standardized.standardized_features;
                k = k_clusters,
                feature_names = baseline_kept_names
            )
            
            baseline_centroids = cluster_res.centroids_standardized
            
            # Map back to 2D
            cluster_map = fill(NaN, nx, ny)
            for (i, idx) in enumerate(flat_valid)
                cluster_map[idx[1], idx[2]] = cluster_res.cluster_id[i]
            end
            tranche_clusters[h] = cluster_map
            
            summary = summarize_threshold_free_vulnerability_clusters(
                cluster_res, standardized.standardized_features;
                feature_names = baseline_kept_names
            )
            CSV.write(joinpath(output_dir, "dynqual_vulnerability_cluster_summary.csv"), summary.cluster_table)
        else
            # Project onto baseline
            f_mat = feat_matrix[:, baseline_kept_indices]
            s_mat = zeros(Float64, size(f_mat))
            for out_idx in 1:length(baseline_kept_names)
                col = f_mat[:, out_idx]
                m = baseline_means[out_idx]
                s = baseline_stds[out_idx]
                for i in 1:size(f_mat, 1)
                    s_mat[i, out_idx] = (col[i] - m) / s
                end
            end
            
            assignments = zeros(Int, size(s_mat, 1))
            for i in 1:size(s_mat, 1)
                best_k = 1
                min_dist = Inf
                for k in 1:k_clusters
                    d = sum(abs2.(s_mat[i, :] .- baseline_centroids[k, :]))
                    if d < min_dist
                        min_dist = d
                        best_k = k
                    end
                end
                assignments[i] = best_k
            end
            
            cluster_map = fill(NaN, nx, ny)
            for (i, idx) in enumerate(flat_valid)
                cluster_map[idx[1], idx[2]] = assignments[i]
            end
            tranche_clusters[h] = cluster_map
        end
    end

    # 7. Metadata
    metadata = Dict(
        "generated_by" => "TwoTimescaleResilience",
        "example_script" => "dynqual_synthetic_isimip_pressure_demo.jl",
        "created_at" => Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ"),
        "bbox" => bbox_override,
        "input_files" => Dict("BOD" => bod_file, "FC" => fc_file, "TDS" => tds_file, "BODload" => bodload_file),
        "input_files_present" => true,
        "memory_optimized" => true,
        "streaming_reader" => true,
        "quantile_estimation" => "deterministic_sampled_log_quantiles",
        "quantile_time_stride" => time_stride,
        "max_quantile_samples" => max_samples,
        "spatial_stride" => spatial_stride,
        "write_netcdf" => write_netcdf,
        "full_monthly_arrays_materialized" => false,
        "memory_initialization" => "zero",
        "pressure_proxies" => [p.pressure_proxy for p in pressure_mapping_weights],
        "tranches" => tranches,
        "n_species" => n_species,
        "cluster_k" => k_clusters,
        "synthetic_pressure_proxies" => true,
        "threshold_free_features" => true
    )
    open(joinpath(output_dir, "dynqual_synthetic_isimip_metadata.json"), "w") do f
        JSON.print(f, metadata, 4)
    end

    # 8. Plotting
    if make_plots
        println("Generating plots...")
        fig_dir = joinpath(output_dir, "figures")
        mkpath(fig_dir)
        
        # Prepare valid mean maps
        safe_div(a, b) = b > 0 ? a / b : NaN32
        
        mean_raw_bod = safe_div.(sum_raw_bod, count_clim)
        mean_raw_fc = safe_div.(sum_raw_fc, count_clim)
        mean_raw_tds = safe_div.(sum_raw_tds, count_clim)
        mean_raw_bodload = safe_div.(sum_raw_bodload, count_clim)
        
        raw_clim_maps = [
            (title="BOD / organic", data=mean_raw_bod),
            (title="Pathogen / FC proxy", data=mean_raw_fc),
            (title="TDS load", data=mean_raw_tds),
            (title="BOD load", data=mean_raw_bodload)
        ]
        make_dynqual_raw_climatology_figure(fig_dir, sub_lons, sub_lats, raw_clim_maps)
        
        mean_org_p = safe_div.(sum_org_p, count_clim)
        mean_pat_p = safe_div.(sum_pat_p, count_clim)
        mean_ion_p = safe_div.(sum_ion_p, count_clim)
        mean_was_p = safe_div.(sum_was_p, count_clim)
        mean_low_p = safe_div.(sum_low_p, count_clim)
        mean_com_p = safe_div.(sum_com_p, count_clim)
        
        derived_clim_maps = [
            (title="Organic O2 Demand Proxy", data=mean_org_p),
            (title="Pathogen Exposure Proxy", data=mean_pat_p),
            (title="Ionic / Salinity Proxy", data=mean_ion_p),
            (title="Wastewater Source Proxy", data=mean_was_p),
            (title="Low-flow Concentration Proxy", data=mean_low_p),
            (title="Combined Wastewater Proxy", data=mean_com_p)
        ]
        make_dynqual_derived_pressure_layers_figure(fig_dir, sub_lons, sub_lats, derived_clim_maps)
        
        if actual_n_tranches >= 2
            make_dynqual_vulnerability_regime_maps_figure(fig_dir, sub_lons, sub_lats, nx, ny, k_clusters, tranche_clusters[1], tranche_clusters[actual_n_tranches])
        end
        
        make_dynqual_regime_explanation_heatmap_figure(fig_dir, k_clusters, baseline_kept_names, baseline_centroids)
    end
    
    println("--- Demo Complete ---")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_dynqual_synthetic_isimip_pressure_demo()
end
