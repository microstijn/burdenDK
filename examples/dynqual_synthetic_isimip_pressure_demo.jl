# dynqual_synthetic_isimip_pressure_demo.jl

using TwoTimescaleResilience
using Statistics
using Dates
using JSON
using CSV
using DataFrames

# Load optional dependencies
try
    using NCDatasets
catch e
    @warn "NCDatasets is required to read NetCDF files. Please add it to your environment."
end

function get_env_or_fallback(env_var::String, fallback::String)
    val = get(ENV, env_var, fallback)
    return isempty(val) ? fallback : val
end

function run_dynqual_synthetic_isimip_pressure_demo(;
    output_dir = joinpath(dirname(@__DIR__), "output", "dynqual_synthetic_isimip_pressure_demo"),
    bbox = (lon = (-25.0, 45.0), lat = (34.0, 72.0))
)
    println("--- Starting DynQual Synthetic ISIMIP Pressure Demo ---")
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

    if make_plots
        try
            @eval using CairoMakie
        catch
            error("CairoMakie is required for plotting but is not available. Install it or set TTR_DYNQUAL_MAKE_PLOTS=false.")
        end
    end

    println("Output Directory: ", output_dir)
    println("BBox: ", bbox_override)
    println("Plotting Enabled: ", make_plots)

    # 2. Local Helper for Reading & Subsetting
    function read_and_subset(file_path::String, var_name::String)
        println("Reading $var_name from $file_path...")
        ds = NCDataset(file_path, "r")

        # Identify variables
        lon_var = haskey(ds, "lon") ? "lon" : (haskey(ds, "longitude") ? "longitude" : "")
        lat_var = haskey(ds, "lat") ? "lat" : (haskey(ds, "latitude") ? "latitude" : "")
        time_var = haskey(ds, "time") ? "time" : ""

        if lon_var == "" || lat_var == ""
            close(ds)
            error("Could not identify lon/lat variables in $file_path")
        end

        lons = ds[lon_var][:]
        lats = ds[lat_var][:]

        # Subset indices
        lon_idx = findall(x -> bbox_override.lon[1] <= x <= bbox_override.lon[2], lons)
        lat_idx = findall(x -> bbox_override.lat[1] <= x <= bbox_override.lat[2], lats)

        if isempty(lon_idx) || isempty(lat_idx)
            close(ds)
            error("Bbox subset resulted in empty grid.")
        end

        sub_lons = lons[lon_idx]
        sub_lats = lats[lat_idx]

        # Find the main variable
        target_var = var_name
        if !haskey(ds, target_var)
            # Try to guess
            for v in keys(ds)
                if v ∉ [lon_var, lat_var, time_var, "crs"]
                    target_var = v
                    break
                end
            end
        end

        data = ds[target_var][lon_idx, lat_idx, :]

        close(ds)

        # Ensure array is Float64 and handle missing/fill values
        data_f = Array{Float64, 3}(undef, size(data))
        for i in eachindex(data)
            v = data[i]
            if ismissing(v) || isnothing(v) || (v isa Number && !isfinite(v))
                data_f[i] = NaN
            else
                data_f[i] = Float64(v)
            end
        end
        return data_f, sub_lons, sub_lats
    end

    raw_bod, lons, lats = read_and_subset(bod_file, "organic")
    raw_fc, _, _ = read_and_subset(fc_file, "pathogen")
    raw_tds, _, _ = read_and_subset(tds_file, "TDSload")
    raw_bodload, _, _ = read_and_subset(bodload_file, "BODload")

    nx, ny, n_time = size(raw_bod)
    println("Data size after subset: nx=$nx, ny=$ny, time=$n_time")

    if n_time != 480
        @warn "Expected 480 months (1980-2019) but got $n_time. Adjusting tranche logic."
    end

    # 3. Robust scaling helper
    scaling_summary = []

    function robust_log_scale_01(data::Array{Float64, 3}, var_name::String)
        valid_vals = Float64[]
        n_total = length(data)
        n_nan = 0

        log_data = similar(data)
        for i in eachindex(data)
            v = data[i]
            if isnan(v)
                n_nan += 1
                log_data[i] = NaN
            else
                v_pos = max(0.0, v)
                log_val = log1p(v_pos)
                log_data[i] = log_val
                push!(valid_vals, log_val)
            end
        end

        frac_missing = n_nan / n_total
        if frac_missing > 0.05
            @warn "Variable $var_name has $(round(frac_missing*100, digits=1))% missing/non-finite values."
        end

        if isempty(valid_vals)
            error("No valid values for variable $var_name")
        end

        p02 = quantile(valid_vals, 0.02)
        p98 = quantile(valid_vals, 0.98)

        if p98 <= p02
            @warn "Robust scaling degenerate for $var_name: p02=$p02, p98=$p98"
            p98 = p02 + 1.0
        end

        scaled_data = similar(data)
        scaled_vals = Float64[]
        for i in eachindex(log_data)
            v = log_data[i]
            if isnan(v)
                scaled_data[i] = 0.0
            else
                s = clamp((v - p02) / (p98 - p02), 0.0, 1.0)
                scaled_data[i] = s
                push!(scaled_vals, s)
            end
        end

        push!(scaling_summary, (
            variable_name = var_name,
            n_total = n_total,
            n_finite = n_total - n_nan,
            n_missing_or_nonfinite = n_nan,
            fraction_missing = frac_missing,
            p02_log1p = p02,
            p98_log1p = p98,
            min_scaled = isempty(scaled_vals) ? NaN : minimum(scaled_vals),
            max_scaled = isempty(scaled_vals) ? NaN : maximum(scaled_vals),
            mean_scaled = isempty(scaled_vals) ? NaN : mean(scaled_vals)
        ))

        return scaled_data, log_data
    end

    scaled_bod, _ = robust_log_scale_01(raw_bod, "organic_monthlyAvg_1980_2019")
    scaled_fc, _ = robust_log_scale_01(raw_fc, "pathogen_monthlyAvg_1980_2019")
    scaled_tds, _ = robust_log_scale_01(raw_tds, "TDSload_monthlyAvg_1980_2019")
    scaled_bodload, _ = robust_log_scale_01(raw_bodload, "BODload_monthlyAvg_1980_2019")

    CSV.write(joinpath(output_dir, "dynqual_scaling_summary.csv"), scaling_summary)

    # 4. Pressure Proxies
    organic_pressure = scaled_bod
    pathogen_pressure = scaled_fc
    ionic_pressure = scaled_tds
    wastewater_pressure = scaled_bodload

    low_flow_pressure = similar(organic_pressure)
    for i in eachindex(organic_pressure)
        bod_val = max(0.0, raw_bod[i])
        bodload_val = max(0.0, raw_bodload[i])
        if isnan(bod_val) || isnan(bodload_val)
            low_flow_pressure[i] = 0.0
        else
            q_proxy = bodload_val / (bod_val + 1e-6)
            low_flow_pressure[i] = clamp(scaled_bod[i] * exp(-0.1 * q_proxy), 0.0, 1.0)
        end
    end

    valid_lf = filter(x -> !isnan(x) && x > 0, low_flow_pressure)
    if !isempty(valid_lf)
        lf_p98 = quantile(valid_lf, 0.98)
        if lf_p98 > 0
            low_flow_pressure .= clamp.(low_flow_pressure ./ lf_p98, 0.0, 1.0)
        end
    end

    combined_pressure = similar(organic_pressure)
    for i in eachindex(combined_pressure)
        combined_pressure[i] = clamp(0.4 * organic_pressure[i] + 0.3 * pathogen_pressure[i] + 0.2 * wastewater_pressure[i] + 0.1 * low_flow_pressure[i], 0.0, 1.0)
    end

    pressure_mapping_weights = [
        (pressure_proxy = "organic_oxygen_demand_proxy", source_variables = "BOD", transform = "robust log 0-1", interpretation = "organic pollution / oxygen-demand stress", primary_deb_axis = "maintenance", secondary_deb_axis = "growth", memory_rho = 0.30, internal_magnification_K = 1.0, notes = "", w_A = 0.0, w_M = 0.8, w_G = 0.2, w_R = 0.0),
        (pressure_proxy = "pathogen_exposure_proxy", source_variables = "FC/pathogen", transform = "robust log 0-1", interpretation = "microbial/pathogen pressure", primary_deb_axis = "maintenance", secondary_deb_axis = "reproduction", memory_rho = 0.10, internal_magnification_K = 1.0, notes = "", w_A = 0.0, w_M = 0.6, w_G = 0.0, w_R = 0.4),
        (pressure_proxy = "ionic_salinity_proxy", source_variables = "TDSload", transform = "robust log 0-1", interpretation = "dissolved-solids / ionic / osmotic pressure proxy", primary_deb_axis = "maintenance", secondary_deb_axis = "assimilation", memory_rho = 0.70, internal_magnification_K = 1.5, notes = "Load, not local concentration.", w_A = 0.3, w_M = 0.7, w_G = 0.0, w_R = 0.0),
        (pressure_proxy = "wastewater_source_proxy", source_variables = "BODload", transform = "robust log 0-1", interpretation = "organic source/load pressure", primary_deb_axis = "growth", secondary_deb_axis = "reproduction", memory_rho = 0.40, internal_magnification_K = 1.0, notes = "", w_A = 0.0, w_M = 0.0, w_G = 0.5, w_R = 0.5),
        (pressure_proxy = "low_flow_concentration_proxy", source_variables = "BOD + BODload", transform = "derived heuristic", interpretation = "relative low-flow / dilution proxy", primary_deb_axis = "maintenance", secondary_deb_axis = "", memory_rho = 0.50, internal_magnification_K = 1.2, notes = "Not real discharge.", w_A = 0.0, w_M = 1.0, w_G = 0.0, w_R = 0.0),
        (pressure_proxy = "combined_wastewater_proxy", source_variables = "composite", transform = "weighted sum", interpretation = "composite wastewater pressure", primary_deb_axis = "maintenance", secondary_deb_axis = "", memory_rho = 0.40, internal_magnification_K = 1.0, notes = "", w_A = 0.2, w_M = 0.4, w_G = 0.2, w_R = 0.2)
    ]
    CSV.write(joinpath(output_dir, "dynqual_pressure_mapping.csv"), pressure_mapping_weights)

    # 5. Model Chain
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

    # To avoid 30 GB memory usage, pre-allocate efficiently
    Q_out = zeros(Float32, nx, ny, n_time, n_species, 1)
    F_out = zeros(Float32, nx, ny, n_time, n_species, 1)

    E_assimilation = zeros(Float32, nx, ny, n_time, 1)
    E_maintenance = zeros(Float32, nx, ny, n_time, 1)
    E_growth = zeros(Float32, nx, ny, n_time, 1)
    E_reproduction = zeros(Float32, nx, ny, n_time, 1)

    println("Running DEB routing and memory...")

    proxies = [
        (layer=organic_pressure, mapping=pressure_mapping_weights[1]),
        (layer=pathogen_pressure, mapping=pressure_mapping_weights[2]),
        (layer=ionic_pressure, mapping=pressure_mapping_weights[3]),
        (layer=wastewater_pressure, mapping=pressure_mapping_weights[4]),
        (layer=low_flow_pressure, mapping=pressure_mapping_weights[5]),
        (layer=combined_pressure, mapping=pressure_mapping_weights[6])
    ]
    n_proxies = length(proxies)

    Threads.@threads for idx in CartesianIndices((1:nx, 1:ny))
        x = idx[1]
        y = idx[2]

        B_state = zeros(Float64, n_proxies)

        for t in 1:n_time
            for p in 1:n_proxies
                rho = proxies[p].mapping.memory_rho
                K_mag = proxies[p].mapping.internal_magnification_K
                P_pt = proxies[p].layer[x, y, t]
                if isnan(P_pt) P_pt = 0.0 end
                B_state[p] = rho * B_state[p] + (1.0 - rho) * K_mag * P_pt
            end

            E_a, E_m, E_g, E_r = 0.0, 0.0, 0.0, 0.0
            for p in 1:n_proxies
                E_proxy = min(1.0, max(0.0, B_state[p]))
                E_a += proxies[p].mapping.w_A * E_proxy
                E_m += proxies[p].mapping.w_M * E_proxy
                E_g += proxies[p].mapping.w_G * E_proxy
                E_r += proxies[p].mapping.w_R * E_proxy
            end

            E_a = min(1.0, max(0.0, E_a))
            E_m = min(1.0, max(0.0, E_m))
            E_g = min(1.0, max(0.0, E_g))
            E_r = min(1.0, max(0.0, E_r))

            E_assimilation[x, y, t, 1] = E_a
            E_maintenance[x, y, t, 1] = E_m
            E_growth[x, y, t, 1] = E_g
            E_reproduction[x, y, t, 1] = E_r

            axes_imp = (assimilation=E_a, maintenance=E_m, growth=E_g, reproduction=E_r)

            for sp in 1:n_species
                resp = compute_adaptive_margin_response_from_impairment(axes_imp, species_params[sp]; mixture_effect_model="grouped_ca_then_ia_axis_effects")
                Q_out[x, y, t, sp, 1] = Float32(resp.Q_t)
                F_out[x, y, t, sp, 1] = Float32(resp.amplification)
            end
        end
    end

    # 6. Tranches & Features
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
            n_months = n_months_per_tranche
        ))
    end
    CSV.write(joinpath(output_dir, "dynqual_tranche_definitions.csv"), tranches)

    println("Extracting threshold-free features and clustering...")
    baseline_kept_indices = Int[]
    baseline_kept_names = String[]
    baseline_means = Float64[]
    baseline_stds = Float64[]
    baseline_centroids = zeros(Float64, 0, 0)
    cluster_labels = String[]

    tranche_clusters = Dict{Int, Any}()

    for h in 1:actual_n_tranches
        start_m = tranches[h].month_start
        end_m = tranches[h].month_end

        Q_slice = Float64.(Q_out[:, :, start_m:end_m, :, :])
        F_slice = Float64.(F_out[:, :, start_m:end_m, :, :])
        E_a_slice = Float64.(E_assimilation[:, :, start_m:end_m, :])
        E_m_slice = Float64.(E_maintenance[:, :, start_m:end_m, :])
        E_g_slice = Float64.(E_growth[:, :, start_m:end_m, :])
        E_r_slice = Float64.(E_reproduction[:, :, start_m:end_m, :])

        resp_tuple = (
            Q_t = Q_slice, F_t = F_slice,
            E_assimilation = E_a_slice, E_maintenance = E_m_slice, E_growth = E_g_slice, E_reproduction = E_r_slice
        )

        feat_res = build_threshold_free_vulnerability_features(
            resp_tuple;
            mixture_model_names = ["grouped_ca_then_ia_axis_effects"],
            preferred_mixture_model = "grouped_ca_then_ia_axis_effects",
            month_values = collect(start_m:end_m)
        )

        if h == 1
            mask = .!occursin.("month_of_max", feat_res.feature_names)
            kept_for_clustering = findall(mask)

            excluded = feat_res.feature_names[.!mask]
            CSV.write(joinpath(output_dir, "dynqual_excluded_from_fixed_reference_clustering.csv"),
                      [(feature_name=f, reason="absolute_time") for f in excluded])

            standardized = standardize_threshold_free_vulnerability_features(
                feat_res.feature_matrix[:, kept_for_clustering],
                feat_res.feature_names[kept_for_clustering]
            )

            baseline_kept_indices = kept_for_clustering[standardized.kept_feature_indices]
            baseline_kept_names = standardized.standardized_feature_names
            baseline_means = standardized.means
            baseline_stds = standardized.stds

            cluster_res = cluster_threshold_free_vulnerability_regimes(
                standardized.standardized_features;
                k = k_clusters,
                feature_names = baseline_kept_names
            )

            baseline_centroids = cluster_res.centroids_standardized
            cluster_labels = label_threshold_free_vulnerability_regimes(baseline_centroids, baseline_kept_names)

            tranche_clusters[h] = cluster_res.cluster_id

            summary = summarize_threshold_free_vulnerability_clusters(
                cluster_res, standardized.standardized_features;
                feature_names = baseline_kept_names
            )
            CSV.write(joinpath(output_dir, "dynqual_vulnerability_cluster_summary.csv"), summary.cluster_table)
        else
            f_mat = feat_res.feature_matrix[:, baseline_kept_indices]
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
            tranche_clusters[h] = assignments
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
        "pressure_proxies" => [p.pressure_proxy for p in pressure_mapping_weights],
        "pressure_mapping_weights" => pressure_mapping_weights,
        "tranches" => tranches,
        "species_selection_mode" => "first_N_from_AmP_Library",
        "n_species" => n_species,
        "cluster_k" => k_clusters,
        "synthetic_pressure_proxies" => true,
        "real_dynqual_spatiotemporal_patterns" => true,
        "real_raster_ingestion_core_feature" => false,
        "threshold_free_features" => true,
        "safe_unsafe_classes" => false,
        "physiological_Z_t" => false,
        "DEBtox_D_t" => false,
        "synergism_antagonism" => false,
        "fitted_interactions" => false
    )
    open(joinpath(output_dir, "dynqual_synthetic_isimip_metadata.json"), "w") do f
        JSON.print(f, metadata, 4)
    end

    # 8. Plotting
    if make_plots
        println("Generating plots...")
        fig_dir = joinpath(output_dir, "figures")
        mkpath(fig_dir)

        @eval begin
            println("  -> Figure 1: Raw climatology")
            fig1 = Figure(size = (1200, 800))
            Label(fig1[0, :], "Real DynQual water-quality fields over Europe, 1980–2019 climatology", font=:bold, fontsize=20)

            axes_data = [
                (title="BOD / organic", data=dropdims(mean($raw_bod, dims=3), dims=3)),
                (title="Pathogen / FC proxy", data=dropdims(mean($raw_fc, dims=3), dims=3)),
                (title="TDS load", data=dropdims(mean($raw_tds, dims=3), dims=3)),
                (title="BOD load", data=dropdims(mean($raw_bodload, dims=3), dims=3))
            ]

            for (i, p_info) in enumerate(axes_data)
                row = cld(i, 2)
                col = mod1(i, 2) * 2 - 1
                ax = Axis(fig1[row, col], title=p_info.title, aspect=DataAspect())
                plot_data = log1p.(max.(0.0, p_info.data))
                hm = heatmap!(ax, $lons, $lats, plot_data, colormap=:viridis)
                Colorbar(fig1[row, col+1], hm, label="log1p")
            end
            save(joinpath($fig_dir, "dynqual_raw_climatology_maps.png"), fig1)

            println("  -> Figure 2: Derived pressure archetypes")
            fig2 = Figure(size = (1500, 1000))
            Label(fig2[0, :], "DynQual-derived synthetic pressure archetypes (Climatological Mean)", font=:bold, fontsize=20)

            der_data = [
                (title="Organic O2 Demand Proxy", data=dropdims(mean($organic_pressure, dims=3), dims=3)),
                (title="Pathogen Exposure Proxy", data=dropdims(mean($pathogen_pressure, dims=3), dims=3)),
                (title="Ionic / Salinity Proxy", data=dropdims(mean($ionic_pressure, dims=3), dims=3)),
                (title="Wastewater Source Proxy", data=dropdims(mean($wastewater_pressure, dims=3), dims=3)),
                (title="Low-flow Concentration Proxy", data=dropdims(mean($low_flow_pressure, dims=3), dims=3)),
                (title="Combined Wastewater Proxy", data=dropdims(mean($combined_pressure, dims=3), dims=3))
            ]

            for (i, p_info) in enumerate(der_data)
                row = cld(i, 3)
                col = mod1(i, 3) * 2 - 1
                ax = Axis(fig2[row, col], title=p_info.title, aspect=DataAspect())
                hm = heatmap!(ax, $lons, $lats, p_info.data, colormap=:plasma, colorrange=(0.0, 1.0))
                Colorbar(fig2[row, col+1], hm, label="Relative Impairment")
            end
            save(joinpath($fig_dir, "dynqual_derived_pressure_layers.png"), fig2)

            if $actual_n_tranches >= 2
                println("  -> Figure 3: Vulnerability regimes")
                fig3 = Figure(size = (1200, 600))
                Label(fig3[0, :], "Threshold-free vulnerability regimes from DynQual-derived pressures", font=:bold, fontsize=20)

                ax_base = Axis(fig3[1, 1], title="Baseline Tranche (1980s) Vulnerability Regimes", aspect=DataAspect())
                map_base = reshape($tranche_clusters[1], $nx, $ny)
                hm_base = heatmap!(ax_base, $lons, $lats, map_base, colormap=:Set1_5, colorrange=(0.5, $k_clusters + 0.5))
                Colorbar(fig3[1, 2], hm_base, ticks=1:$k_clusters)

                ax_rec = Axis(fig3[1, 3], title="Recent Tranche (2010s) Vulnerability Regimes", aspect=DataAspect())
                map_rec = reshape($tranche_clusters[$actual_n_tranches], $nx, $ny)
                hm_rec = heatmap!(ax_rec, $lons, $lats, map_rec, colormap=:Set1_5, colorrange=(0.5, $k_clusters + 0.5))
                Colorbar(fig3[1, 4], hm_rec, ticks=1:$k_clusters)

                save(joinpath($fig_dir, "dynqual_vulnerability_regime_maps.png"), fig3)
            end

            println("  -> Figure 4: Regime explanation heatmap")
            fig4 = Figure(size = (1000, 800))
            Label(fig4[0, :], "What distinguishes the vulnerability regimes?", font=:bold, fontsize=20)

            ax_heat = Axis(fig4[1, 1],
                xticks = (1:$k_clusters, ["Cluster $i" for i in 1:$k_clusters]),
                yticks = (1:length($baseline_kept_names), $baseline_kept_names),
                xticklabelrotation = pi/4
            )

            hm_heat = heatmap!(ax_heat, 1:$k_clusters, 1:length($baseline_kept_names), $baseline_centroids', colormap=:RdBu, colorrange=(-3.0, 3.0))
            Colorbar(fig4[1, 2], hm_heat, label="Standardised Feature Value")

            for c in 1:$k_clusters
                for f in 1:length($baseline_kept_names)
                    val = $baseline_centroids[c, f]
                    text!(ax_heat, string(round(val, digits=2)), position=(c, f),
                          align=(:center, :center), color=abs(val) > 1.5 ? :white : :black, fontsize=12)
                end
            end

            save(joinpath($fig_dir, "dynqual_regime_explanation_heatmap.png"), fig4)
        end
    end

    println("--- Demo Complete ---")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_dynqual_synthetic_isimip_pressure_demo()
end
