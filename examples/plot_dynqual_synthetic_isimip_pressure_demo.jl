# plot_dynqual_synthetic_isimip_pressure_demo.jl

using Statistics
using Dates
using JSON
using CSV
using DataFrames
using NCDatasets
using Downloads
using GeoJSON

# --- Helper Functions ---

function finite_max_abs(A)
    filtered = filter(x -> !ismissing(x) && isfinite(x), A)
    return isempty(filtered) ? NaN : maximum(abs, filtered)
end

function feature_index_by_name(df_meta, feature_name)
    idx = findfirst(df_meta.feature_name .== feature_name)
    if isnothing(idx)
        return nothing
    end
    return df_meta.feature_index[idx]
end

function select_available_features(df_meta, preferred_names)
    indices = Int[]
    names = String[]
    descriptors = String[]

    for name in preferred_names
        idx = feature_index_by_name(df_meta, name)
        if !isnothing(idx)
            push!(indices, idx)
            push!(names, name)

            row_idx = findfirst(df_meta.feature_name .== name)
            desc = df_meta.feature_descriptor[row_idx]
            push!(descriptors, desc)
        end
    end
    return indices, names, descriptors
end

function transition_matrix(cluster_base, cluster_recent, k; valid_mask=nothing)
    T = zeros(Float64, k, k)
    counts = zeros(Int, k, k)

    for i in eachindex(cluster_base)
        if isnothing(valid_mask) || valid_mask[i]
            cb = cluster_base[i]
            cr = cluster_recent[i]
            if isfinite(cb) && isfinite(cr) && cb > 0 && cr > 0 && cb <= k && cr <= k
                counts[Int(cb), Int(cr)] += 1
            end
        end
    end

    for i in 1:k
        row_sum = sum(counts[i, :])
        if row_sum > 0
            T[i, :] .= counts[i, :] ./ row_sum
        end
    end

    return T, counts
end

function changed_regime_map(cluster_base, cluster_recent; valid_mask=nothing)
    changed = fill(NaN, size(cluster_base))
    total_valid = 0
    total_changed = 0

    for i in eachindex(cluster_base)
        if isnothing(valid_mask) || valid_mask[i]
            cb = cluster_base[i]
            cr = cluster_recent[i]
            if isfinite(cb) && isfinite(cr) && cb > 0 && cr > 0
                total_valid += 1
                if cb != cr
                    changed[i] = 1.0
                    total_changed += 1
                else
                    changed[i] = 0.0
                end
            end
        end
    end

    fraction_changed = total_valid > 0 ? total_changed / total_valid : 0.0
    return changed, fraction_changed
end

function feature_delta_map(feature_map, feature_index, baseline_idx, recent_idx)
    base_slice = feature_map[:, :, feature_index, baseline_idx]
    recent_slice = feature_map[:, :, feature_index, recent_idx]
    return recent_slice .- base_slice
end

function clusterwise_feature_delta(feature_map, cluster_recent, feature_indices, baseline_idx, recent_idx, k; valid_mask=nothing)
    delta_matrix = zeros(Float64, length(feature_indices), k)

    for (f_i, feat_idx) in enumerate(feature_indices)
        delta_map = feature_delta_map(feature_map, feat_idx, baseline_idx, recent_idx)

        for c in 1:k
            cell_deltas = Float64[]
            for i in eachindex(cluster_recent)
                if isnothing(valid_mask) || valid_mask[i]
                    if isfinite(cluster_recent[i]) && cluster_recent[i] == c
                        d = delta_map[i]
                        if !ismissing(d) && isfinite(d)
                            push!(cell_deltas, d)
                        end
                    end
                end
            end

            if !isempty(cell_deltas)
                delta_matrix[f_i, c] = mean(cell_deltas)
            else
                delta_matrix[f_i, c] = 0.0
            end
        end
    end

    return delta_matrix
end

# --- Plotting Functions ---

function draw_geojson_ring!(ax, ring; color=:black, linewidth=0.45)
    n = length(ring)
    n < 2 && return

    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)

    for i in 1:n
        xs[i] = Float64(ring[i][1])
        ys[i] = Float64(ring[i][2])
    end

    lines!(ax, xs, ys; color=color, linewidth=linewidth)
end

function add_borders!(ax, countries)
    if countries isa AbstractDict && haskey(countries, "features")
        for feature in countries["features"]
            geom = get(feature, "geometry", nothing)
            isnothing(geom) && continue

            geom_type = get(geom, "type", "")
            coords = get(geom, "coordinates", nothing)
            isnothing(coords) && continue

            if geom_type == "Polygon"
                for ring in coords
                    draw_geojson_ring!(ax, ring)
                end

            elseif geom_type == "MultiPolygon"
                for polygon in coords
                    for ring in polygon
                        draw_geojson_ring!(ax, ring)
                    end
                end
            end
        end
    else
        lines!(ax, GeoMakie.coastlines(); color=:black, linewidth=0.8)
    end
end

function make_dynqual_raw_climatology_figure(fig_dir, lons, lats, lon_lims, lat_lims, raw_clim_maps, countries)
    try
        @eval using CairoMakie, GeoMakie
    catch
        error("CairoMakie and GeoMakie are required for plotting but are not available.")
    end

    println("  -> Figure 1: Raw climatology")
    fig1 = Figure(size = (1200, 800))
    Label(fig1[0, :], "Real DynQual water-quality fields over Europe, 1980–2019 climatology", font=:bold, fontsize=20)

    for (i, p_info) in enumerate(raw_clim_maps)
        row = cld(i, 2)
        col = mod1(i, 2) * 2 - 1
        ax = GeoAxis(fig1[row, col], title=p_info.title, dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
        hm = heatmap!(ax, lons, lats, p_info.data, colormap=:viridis)
        add_borders!(ax, countries)
        Colorbar(fig1[row, col+1], hm, label="log1p")
        colsize!(fig1.layout, col+1, Auto())
    end
    save(joinpath(fig_dir, "dynqual_raw_climatology_maps.png"), fig1)
end

function make_dynqual_derived_pressure_layers_figure(fig_dir, lons, lats, lon_lims, lat_lims, derived_clim_maps, countries)
    try
        @eval using CairoMakie, GeoMakie
    catch
        error("CairoMakie and GeoMakie are required for plotting but are not available.")
    end

    println("  -> Figure 2: Derived pressure archetypes")
    fig2 = Figure(size = (1500, 1000))
    Label(fig2[0, :], "DynQual-derived synthetic pressure archetypes (Climatological Mean)", font=:bold, fontsize=20)

    for (i, p_info) in enumerate(derived_clim_maps)
        row = cld(i, 3)
        col = mod1(i, 3) * 2 - 1
        ax = GeoAxis(fig2[row, col], title=p_info.title, dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
        hm = heatmap!(ax, lons, lats, p_info.data, colormap=:plasma, colorrange=(0.0, 1.0))
        add_borders!(ax, countries)
        Colorbar(fig2[row, col+1], hm, label="Relative Impairment")
        colsize!(fig2.layout, col+1, Auto())
    end
    save(joinpath(fig_dir, "dynqual_derived_pressure_layers.png"), fig2)
end

function make_dynqual_vulnerability_regime_maps_figure(fig_dir, lons, lats, nx, ny, lon_lims, lat_lims, k_clusters, map_base, map_rec, countries)
    try
        @eval using CairoMakie, GeoMakie
    catch
        error("CairoMakie and GeoMakie are required for plotting but are not available.")
    end

    println("  -> Figure 3: Vulnerability regimes")
    fig3 = Figure(size = (1200, 600))
    Label(fig3[0, :], "Threshold-free vulnerability regimes from DynQual-derived pressures", font=:bold, fontsize=20)

    ax_base = GeoAxis(fig3[1, 1], title="Baseline Tranche (1980s) Vulnerability Regimes", dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
    map_base_2d = ndims(map_base) == 1 ? reshape(map_base, nx, ny) : map_base
    hm_base = heatmap!(ax_base, lons, lats, map_base_2d, colormap=:Set1_5, colorrange=(0.5, k_clusters + 0.5))
    add_borders!(ax_base, countries)
    Colorbar(fig3[1, 2], hm_base, ticks=1:k_clusters)
    colsize!(fig3.layout, 2, Auto())

    ax_rec = GeoAxis(fig3[1, 3], title="Recent Tranche (2010s) Vulnerability Regimes", dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
    map_rec_2d = ndims(map_rec) == 1 ? reshape(map_rec, nx, ny) : map_rec
    hm_rec = heatmap!(ax_rec, lons, lats, map_rec_2d, colormap=:Set1_5, colorrange=(0.5, k_clusters + 0.5))
    add_borders!(ax_rec, countries)
    Colorbar(fig3[1, 4], hm_rec, ticks=1:k_clusters)
    colsize!(fig3.layout, 4, Auto())

    save(joinpath(fig_dir, "dynqual_vulnerability_regime_maps.png"), fig3)
end

function make_dynqual_baseline_recent_comparison_figure(fig_dir, lons, lats, nx, ny, lon_lims, lat_lims, k_clusters, cluster_base, cluster_recent, df_meta, feature_map, valid_cell_mask, countries)
    try
        @eval using CairoMakie, GeoMakie
    catch
        error("CairoMakie and GeoMakie are required for plotting but are not available.")
    end

    println("  -> Figure A: Baseline vs Recent Comparison")
    fig = Figure(size = (2800, 1600))
    Label(fig[0, 1], "Baseline vs Recent Vulnerability Regimes Comparison", font=:bold, fontsize=24)

    top_grid = fig[1, 1] = GridLayout()
    bottom_grid = fig[2, 1] = GridLayout()

    # --- Top Row: Baseline Map, Recent Map, Changed Regime Map ---
    ax_base = GeoAxis(top_grid[1, 1], title="Baseline Tranche", dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
    cb_2d = ndims(cluster_base) == 1 ? reshape(cluster_base, nx, ny) : cluster_base
    hm_base = heatmap!(ax_base, lons, lats, cb_2d, colormap=:Set1_5, colorrange=(0.5, k_clusters + 0.5))
    add_borders!(ax_base, countries)
    Colorbar(top_grid[1, 2], hm_base, ticks=1:k_clusters)

    ax_rec = GeoAxis(top_grid[1, 3], title="Recent Tranche", dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
    cr_2d = ndims(cluster_recent) == 1 ? reshape(cluster_recent, nx, ny) : cluster_recent
    hm_rec = heatmap!(ax_rec, lons, lats, cr_2d, colormap=:Set1_5, colorrange=(0.5, k_clusters + 0.5))
    add_borders!(ax_rec, countries)
    Colorbar(top_grid[1, 4], hm_rec, ticks=1:k_clusters)

    changed_map, frac_changed = changed_regime_map(cluster_base, cluster_recent; valid_mask=valid_cell_mask)
    changed_2d = ndims(changed_map) == 1 ? reshape(changed_map, nx, ny) : changed_map
    ax_change = GeoAxis(top_grid[1, 5], title="Changed regime cells\n$(round(frac_changed * 100, digits=1))% of valid cells changed", dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
    cmap_change = cgrad([:lightgrey, :purple], 2, categorical=true)
    hm_change = heatmap!(ax_change, lons, lats, changed_2d, colormap=cmap_change, colorrange=(-0.5, 1.5))
    add_borders!(ax_change, countries)
    Colorbar(top_grid[1, 6], hm_change, ticks=( [0.0, 1.0], ["Unchanged", "Changed"] ))

    colsize!(top_grid, 2, Auto())
    colsize!(top_grid, 4, Auto())
    colsize!(top_grid, 6, Auto())

    # --- Bottom Row: Transition Heatmap, Feature Delta Map 1, Feature Delta Map 2 ---
    T, counts = transition_matrix(cluster_base, cluster_recent, k_clusters; valid_mask=valid_cell_mask)
    
    # Enforce strict String array to prevent Makie Vector{Any} StackOverflow
    x_ticks_labels = String["Recent $i" for i in 1:k_clusters]
    y_ticks_labels = String["Baseline $i" for i in 1:k_clusters]
    
    ax_trans = Axis(bottom_grid[1, 1], title="Cluster transitions: baseline to recent",
        xticks=(1:k_clusters, x_ticks_labels),
        yticks=(1:k_clusters, y_ticks_labels)
    )
    hm_trans = heatmap!(ax_trans, 1:k_clusters, 1:k_clusters, T', colormap=:viridis, colorrange=(0.0, 1.0))
    Colorbar(bottom_grid[1, 2], hm_trans, label="fraction of baseline cluster")

    if k_clusters <= 8
        for i in 1:k_clusters
            for j in 1:k_clusters
                val = T[i, j]
                if val > 0
                    text!(ax_trans, string(round(val * 100, digits=1), "%"), position=(j, i),
                          align=(:center, :center), color=val > 0.5 ? :white : :black, fontsize=12)
                end
            end
        end
    end

    ax_feat1 = Axis(bottom_grid[1, 3], title="No first feature-delta map available")
    hidespines!(ax_feat1)
    hidedecorations!(ax_feat1)

    ax_feat2 = Axis(bottom_grid[1, 5], title="No second feature-delta map available")
    hidespines!(ax_feat2)
    hidedecorations!(ax_feat2)

    if !isnothing(feature_map)
        group1_prefs = ["p95_F_grouped", "mean_F_grouped", "p95_Q_grouped", "mean_Q_grouped",
                        "mean_E_maintenance_grouped", "mean_E_growth_grouped", "mean_E_reproduction_grouped", "mean_E_assimilation_grouped", "axis_entropy"]
        idx1_arr, name1_arr, desc1_arr = select_available_features(df_meta, group1_prefs)

        if !isempty(idx1_arr)
            feat1_idx = idx1_arr[1]
            feat1_desc = desc1_arr[1]

            baseline_idx = 1
            recent_idx = size(feature_map, 4)
            delta1 = feature_delta_map(feature_map, feat1_idx, baseline_idx, recent_idx)

            max_v = finite_max_abs(delta1)
            if isfinite(max_v) && max_v > 0
                delta1_2d = ndims(delta1) == 1 ? reshape(delta1, nx, ny) : delta1
                ax_feat1_geo = GeoAxis(bottom_grid[1, 3], title="Δ $(feat1_desc)", dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
                hm_feat1 = heatmap!(ax_feat1_geo, lons, lats, delta1_2d, colormap=:balance, colorrange=(-max_v, max_v))
                add_borders!(ax_feat1_geo, countries)
                Colorbar(bottom_grid[1, 4], hm_feat1)
            end

            group2_prefs = filter(x -> x != name1_arr[1], [
                "mean_E_maintenance_grouped", "mean_E_growth_grouped", "mean_E_reproduction_grouped", "mean_E_assimilation_grouped",
                "p95_F_grouped", "mean_F_grouped", "p95_Q_grouped", "mean_Q_grouped", "axis_entropy"
            ])
            idx2_arr, name2_arr, desc2_arr = select_available_features(df_meta, group2_prefs)

            if !isempty(idx2_arr)
                feat2_idx = idx2_arr[1]
                feat2_desc = desc2_arr[1]

                delta2 = feature_delta_map(feature_map, feat2_idx, baseline_idx, recent_idx)
                max_v2 = finite_max_abs(delta2)
                if isfinite(max_v2) && max_v2 > 0
                    delta2_2d = ndims(delta2) == 1 ? reshape(delta2, nx, ny) : delta2
                    ax_feat2_geo = GeoAxis(bottom_grid[1, 5], title="Δ $(feat2_desc)", dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
                    hm_feat2 = heatmap!(ax_feat2_geo, lons, lats, delta2_2d, colormap=:balance, colorrange=(-max_v2, max_v2))
                    add_borders!(ax_feat2_geo, countries)
                    Colorbar(bottom_grid[1, 6], hm_feat2)
                end
            end
        end
    end

    colsize!(bottom_grid, 2, Auto())
    colsize!(bottom_grid, 4, Auto())
    colsize!(bottom_grid, 6, Auto())

    save(joinpath(fig_dir, "dynqual_baseline_recent_comparison.png"), fig)
end

function make_dynqual_cluster_transition_heatmap_figure(fig_dir, k_clusters, cluster_base, cluster_recent, valid_cell_mask)
    try
        @eval using CairoMakie
    catch
        error("CairoMakie is required for plotting but is not available.")
    end

    println("  -> Figure B: Cluster transition heatmap")
    fig = Figure(size = (800, 600))
    Label(fig[0, :], "Cluster transitions: baseline to recent", font=:bold, fontsize=20)

    T, counts = transition_matrix(cluster_base, cluster_recent, k_clusters; valid_mask=valid_cell_mask)

    x_ticks_labels = String["Recent $i" for i in 1:k_clusters]
    y_ticks_labels = String["Baseline $i" for i in 1:k_clusters]

    ax_trans = Axis(fig[1, 1],
        xticks=(1:k_clusters, x_ticks_labels),
        yticks=(1:k_clusters, y_ticks_labels)
    )

    hm_trans = heatmap!(ax_trans, 1:k_clusters, 1:k_clusters, T', colormap=:viridis, colorrange=(0.0, 1.0))
    Colorbar(fig[1, 2], hm_trans, label="fraction of baseline cluster")
    colsize!(fig.layout, 2, Auto())

    if k_clusters <= 8
        for i in 1:k_clusters
            for j in 1:k_clusters
                val = T[i, j]
                if val > 0
                    text!(ax_trans, string(round(val * 100, digits=1), "%"), position=(j, i),
                          align=(:center, :center), color=val > 0.5 ? :white : :black, fontsize=12)
                end
            end
        end
    end

    save(joinpath(fig_dir, "dynqual_cluster_transition_baseline_to_recent.png"), fig)
end

function make_dynqual_feature_delta_maps(fig_dir, lons, lats, nx, ny, lon_lims, lat_lims, df_meta, feature_map, skipped_figures, countries)
    try
        @eval using CairoMakie, GeoMakie
    catch
        error("CairoMakie and GeoMakie are required for plotting but are not available.")
    end

    println("  -> Figure C: Feature delta maps")
    preferred_names = [
        "p95_F_grouped",
        "mean_F_grouped",
        "p95_Q_grouped",
        "mean_Q_grouped",
        "mean_E_maintenance_grouped",
        "mean_E_growth_grouped",
        "mean_E_reproduction_grouped",
        "mean_E_assimilation_grouped",
        "axis_entropy"
    ]

    generated = String[]
    baseline_idx = 1
    recent_idx = size(feature_map, 4)

    for name in preferred_names
        idx = feature_index_by_name(df_meta, name)
        if isnothing(idx)
            @warn "Skipping $(name) delta map because it is not present in feature metadata."
            push!(skipped_figures, "dynqual_$(name)_delta_map.png: missing feature")
            continue
        end

        delta = feature_delta_map(feature_map, idx, baseline_idx, recent_idx)
        max_v = finite_max_abs(delta)

        if !isfinite(max_v) || max_v == 0
            @warn "Skipping $(name) delta map because max absolute delta is 0 or non-finite."
            push!(skipped_figures, "dynqual_$(name)_delta_map.png: no variation")
            continue
        end

        fig = Figure(size = (800, 600))
        desc = df_meta.feature_descriptor[findfirst(df_meta.feature_name .== name)]

        ax = GeoAxis(fig[1, 1], title="Δ $(desc) (Recent - Baseline)", dest="+proj=longlat", xgridvisible=false, ygridvisible=false, limits=(lon_lims, lat_lims))
        delta_2d = ndims(delta) == 1 ? reshape(delta, nx, ny) : delta

        hm = heatmap!(ax, lons, lats, delta_2d, colormap=:balance, colorrange=(-max_v, max_v))
        add_borders!(ax, countries)
        Colorbar(fig[1, 2], hm)
        colsize!(fig.layout, 2, Auto())

        filename = "dynqual_$(name)_delta_map.png"
        save(joinpath(fig_dir, filename), fig)
        push!(generated, filename)
    end

    return generated
end

function make_dynqual_feature_delta_by_recent_cluster_figure(fig_dir, k_clusters, cluster_base, cluster_recent, df_meta, feature_map, valid_cell_mask)
    try
        @eval using CairoMakie
    catch
        error("CairoMakie is required for plotting but is not available.")
    end

    println("  -> Figure D: Feature-delta heatmap by recent cluster")

    preferred_names = [
        "p95_F_grouped",
        "mean_F_grouped",
        "p95_Q_grouped",
        "mean_Q_grouped",
        "mean_E_assimilation_grouped",
        "mean_E_maintenance_grouped",
        "mean_E_growth_grouped",
        "mean_E_reproduction_grouped",
        "axis_entropy"
    ]

    indices, names, descriptors = select_available_features(df_meta, preferred_names)

    if isempty(indices)
        @warn "No valid features found for feature-delta by cluster heatmap."
        return false
    end

    baseline_idx = 1
    recent_idx = size(feature_map, 4)

    if isnothing(valid_cell_mask)
        valid_mask = map(i -> isfinite(cluster_base[i]) && isfinite(cluster_recent[i]) && cluster_base[i] > 0 && cluster_recent[i] > 0, eachindex(cluster_base))
    else
        valid_mask = valid_cell_mask
    end

    delta_matrix = clusterwise_feature_delta(feature_map, cluster_recent, indices, baseline_idx, recent_idx, k_clusters, valid_mask=valid_mask)

    fig = Figure(size = (1000, 800))
    Label(fig[0, :], "Mean feature change by Recent Cluster", font=:bold, fontsize=20)

    n_features = length(indices)

    trunc_desc = String[length(d) > 30 ? d[1:27] * "..." : d for d in descriptors]
    x_ticks_labels = String["Cluster $i" for i in 1:k_clusters]

    ax_heat = Axis(fig[1, 1],
        xticks = (1:k_clusters, x_ticks_labels),
        yticks = (1:n_features, trunc_desc),
        xticklabelrotation = pi/4
    )

    max_v = finite_max_abs(delta_matrix)
    colorrange = max_v > 0 ? (-max_v, max_v) : (-1.0, 1.0)

    hm_heat = heatmap!(ax_heat, 1:k_clusters, 1:n_features, delta_matrix', colormap=:balance, colorrange=colorrange)
    Colorbar(fig[1, 2], hm_heat, label="Mean Δ (Recent - Baseline)")
    colsize!(fig.layout, 2, Auto())

    if k_clusters * n_features <= 80
        for i in 1:n_features
            for j in 1:k_clusters
                val = delta_matrix[i, j]
                if !ismissing(val) && isfinite(val)
                    text!(ax_heat, string(round(val, digits=3)), position=(j, i),
                          align=(:center, :center), color=abs(val) > (max_v * 0.6) ? :white : :black, fontsize=12)
                end
            end
        end
    end

    save(joinpath(fig_dir, "dynqual_feature_delta_by_recent_cluster.png"), fig)
    return true
end

function make_dynqual_regime_explanation_heatmap_figure(fig_dir, k_clusters, baseline_kept_names, baseline_centroids, df_meta)
    try
        @eval using CairoMakie
    catch
        error("CairoMakie is required for plotting but is not available.")
    end

    println("  -> Figure 4: Regime explanation heatmap")
    fig4 = Figure(size = (1000, 800))
    Label(fig4[0, :], "What distinguishes the vulnerability regimes?", font=:bold, fontsize=20)

    preferred_names = [
        "mean_E_assimilation_grouped",
        "mean_E_maintenance_grouped",
        "mean_E_growth_grouped",
        "mean_E_reproduction_grouped",
        "mean_Q_grouped",
        "p95_Q_grouped",
        "mean_F_grouped",
        "p95_F_grouped",
        "axis_entropy"
    ]

    selected_indices = Int[]
    selected_descriptors = String[]

    for name in preferred_names
        idx = findfirst(df_meta.feature_name .== name)
        if !isnothing(idx) && df_meta.used_for_clustering[idx]
            desc = df_meta.feature_descriptor[idx]
            centroid_row_idx = findfirst(baseline_kept_names .== desc)
            if !isnothing(centroid_row_idx)
                push!(selected_indices, centroid_row_idx)
                trunc_desc = length(desc) > 30 ? desc[1:27] * "..." : desc
                push!(selected_descriptors, trunc_desc)
            end
        end
    end

    if isempty(selected_indices)
        selected_indices = 1:length(baseline_kept_names)
        selected_descriptors = String[baseline_kept_names...]
    end

    n_features = length(selected_indices)
    filtered_centroids = baseline_centroids[:, selected_indices]

    x_ticks_labels = String["Cluster $i" for i in 1:k_clusters]

    ax_heat = Axis(fig4[1, 1],
        xticks = (1:k_clusters, x_ticks_labels),
        yticks = (1:n_features, selected_descriptors),
        xticklabelrotation = pi/4
    )

    max_v = finite_max_abs(filtered_centroids)
    colorrange = max_v > 0 ? (-max_v, max_v) : (-1.0, 1.0)

    hm_heat = heatmap!(ax_heat, 1:k_clusters, 1:n_features, filtered_centroids, colormap=:balance, colorrange=colorrange)

    if k_clusters * n_features <= 80
        for i in 1:k_clusters
            for j in 1:n_features
                val = filtered_centroids[i, j]
                if !ismissing(val) && isfinite(val)
                    text!(ax_heat, string(round(val, digits=2)), position=(i, j),
                          align=(:center, :center), color=abs(val) > (max_v * 0.6) ? :white : :black, fontsize=12)
                end
            end
        end
    end

    Colorbar(fig4[1, 2], hm_heat, label="Standardized Feature Value")
    colsize!(fig4.layout, 2, Auto())

    save(joinpath(fig_dir, "dynqual_regime_explanation_heatmap.png"), fig4)
end

function run_plot_dynqual_synthetic_isimip_pressure_demo(;
    output_dir = get(ENV, "TTR_DYNQUAL_DEMO_OUTPUT_DIR", joinpath(dirname(@__DIR__), "output", "dynqual_synthetic_isimip_pressure_demo")),
    figures_dir = get(ENV, "TTR_DYNQUAL_FIGURES_DIR", joinpath(get(ENV, "TTR_DYNQUAL_DEMO_OUTPUT_DIR", joinpath(dirname(@__DIR__), "output", "dynqual_synthetic_isimip_pressure_demo")), "figures"))
)
    println("--- Lightweight Plotting: DynQual Synthetic ISIMIP Pressure Demo ---")

    cache_path = joinpath(output_dir, "dynqual_demo_cache.nc")
    if !isfile(cache_path)
        cache_path = joinpath(output_dir, "dynqual_synthetic_isimip_analysis_cache.nc")
    end

    if !isfile(cache_path)
        error("Cache NetCDF missing: Run examples/dynqual_synthetic_isimip_pressure_demo.jl first, or enable TTR_DYNQUAL_WRITE_CACHE=true.")
    end

    metadata_json_path = joinpath(output_dir, "dynqual_synthetic_isimip_metadata.json")
    if !isfile(metadata_json_path)
        error("Metadata JSON missing: Run examples/dynqual_synthetic_isimip_pressure_demo.jl first.")
    end

    feature_metadata_csv = joinpath(output_dir, "dynqual_feature_metadata.csv")
    if !isfile(feature_metadata_csv)
        error("Feature metadata CSV missing: Run examples/dynqual_synthetic_isimip_pressure_demo.jl first.")
    end

    mkpath(figures_dir)
    println("Reading from cache: $cache_path")

    # Download GeoJSON country borders (cached locally to prevent re-downloading)
    country_geojson_url = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson"
    country_geojson_path = joinpath(output_dir, "ne_50m_admin_0_countries.geojson")
    
    local countries = nothing
    try
        if !isfile(country_geojson_path)
            println("Downloading 50m country borders for overlay...")
            Downloads.download(country_geojson_url, country_geojson_path)
        end
        countries = GeoJSON.read(read(country_geojson_path))
    catch e
        @warn "Failed to download/read GeoJSON country borders. Falling back to simple coastlines."
    end

    # Extract tightly bound, strongly-typed variables using a NamedTuple.
    # This entirely prevents Julia from "boxing" local variables into `Any` types,
    # solving the root cause of the Makie StackOverflow.
    cache = NCDataset(cache_path, "r") do ds
        clean_map(x) = Float64.(coalesce.(x, NaN))
        
        return (
            lons = Float64.(ds["lon"][:]),
            lats = Float64.(ds["lat"][:]),
            nx = Int(length(ds["lon"])),
            ny = Int(length(ds["lat"])),
            actual_n_tranches = Int(length(ds["tranche"])),
            k_clusters = Int(length(ds["cluster"])),

            mean_raw_bod = clean_map(ds["raw_BOD_climatology"][:,:]),
            mean_raw_fc = clean_map(ds["raw_pathogen_climatology"][:,:]),
            mean_raw_tds = clean_map(ds["raw_TDSload_climatology"][:,:]),
            mean_raw_bodload = clean_map(ds["raw_BODload_climatology"][:,:]),

            mean_org_p = clean_map(ds["organic_pressure_climatology"][:,:]),
            mean_pat_p = clean_map(ds["pathogen_pressure_climatology"][:,:]),
            mean_ion_p = clean_map(ds["ionic_pressure_climatology"][:,:]),
            mean_was_p = clean_map(ds["wastewater_source_pressure_climatology"][:,:]),
            mean_low_p = clean_map(ds["low_flow_concentration_pressure_climatology"][:,:]),
            mean_com_p = clean_map(ds["combined_wastewater_pressure_climatology"][:,:]),

            cluster_id_all = clean_map(ds["cluster_id"][:, :, :]),
            baseline_centroids_all = clean_map(ds["baseline_centroids"][:, :]),

            feature_map_all = haskey(ds, "feature_map") ? clean_map(ds["feature_map"][:, :, :, :]) : nothing,
            valid_cell_mask = haskey(ds, "valid_cell_mask") ? coalesce.(ds["valid_cell_mask"][:, :], false) : nothing
        )
    end
    
    lon_lims = extrema(filter(isfinite, cache.lons))
    lat_lims = extrema(filter(isfinite, cache.lats))

    df_meta = CSV.read(feature_metadata_csv, DataFrame)
    df_clust = df_meta[df_meta.used_for_clustering .== true, :]
    baseline_kept_names = String[df_clust.feature_descriptor...]

    generated_figures = String[]
    skipped_figures = String[]

    raw_clim_maps = [
        (title="BOD / organic", data=cache.mean_raw_bod),
        (title="Pathogen / FC proxy", data=cache.mean_raw_fc),
        (title="TDS load", data=cache.mean_raw_tds),
        (title="BOD load", data=cache.mean_raw_bodload)
    ]
    

    make_dynqual_raw_climatology_figure(figures_dir, cache.lons, cache.lats, lon_lims, lat_lims, raw_clim_maps, countries)
    push!(generated_figures, "dynqual_raw_climatology_maps.png")

    derived_clim_maps = [
        (title="Organic O2 Demand Proxy", data=cache.mean_org_p),
        (title="Pathogen Exposure Proxy", data=cache.mean_pat_p),
        (title="Ionic / Salinity Proxy", data=cache.mean_ion_p),
        (title="Wastewater Source Proxy", data=cache.mean_was_p),
        (title="Low-flow Concentration Proxy", data=cache.mean_low_p),
        (title="Combined Wastewater Proxy", data=cache.mean_com_p)
    ]
    make_dynqual_derived_pressure_layers_figure(figures_dir, cache.lons, cache.lats, lon_lims, lat_lims, derived_clim_maps, countries)
    push!(generated_figures, "dynqual_derived_pressure_layers.png")

    if cache.actual_n_tranches >= 2
        cluster_base = cache.cluster_id_all[:, :, 1]
        cluster_recent = cache.cluster_id_all[:, :, cache.actual_n_tranches]
        
        make_dynqual_vulnerability_regime_maps_figure(figures_dir, cache.lons, cache.lats, cache.nx, cache.ny, lon_lims, lat_lims, cache.k_clusters, cluster_base, cluster_recent, countries)
        push!(generated_figures, "dynqual_vulnerability_regime_maps.png")

        make_dynqual_baseline_recent_comparison_figure(figures_dir, cache.lons, cache.lats, cache.nx, cache.ny, lon_lims, lat_lims, cache.k_clusters, cluster_base, cluster_recent, df_meta, cache.feature_map_all, cache.valid_cell_mask, countries)
        push!(generated_figures, "dynqual_baseline_recent_comparison.png")

        make_dynqual_cluster_transition_heatmap_figure(figures_dir, cache.k_clusters, cluster_base, cluster_recent, cache.valid_cell_mask)
        push!(generated_figures, "dynqual_cluster_transition_baseline_to_recent.png")

        if !isnothing(cache.feature_map_all)
            generated_delta_maps = make_dynqual_feature_delta_maps(figures_dir, cache.lons, cache.lats, cache.nx, cache.ny, lon_lims, lat_lims, df_meta, cache.feature_map_all, skipped_figures, countries)
            append!(generated_figures, generated_delta_maps)

            if make_dynqual_feature_delta_by_recent_cluster_figure(figures_dir, cache.k_clusters, cluster_base, cluster_recent, df_meta, cache.feature_map_all, cache.valid_cell_mask)
                push!(generated_figures, "dynqual_feature_delta_by_recent_cluster.png")
            else
                push!(skipped_figures, "dynqual_feature_delta_by_recent_cluster.png: no valid features")
            end
        else
            @warn "Skipping feature delta maps and heatmaps because feature_map is missing from cache."
            push!(skipped_figures, "Feature Delta Maps: feature_map missing")
            push!(skipped_figures, "dynqual_feature_delta_by_recent_cluster.png: feature_map missing")
        end
    end

    make_dynqual_regime_explanation_heatmap_figure(figures_dir, cache.k_clusters, baseline_kept_names, cache.baseline_centroids_all, df_meta)
    push!(generated_figures, "dynqual_regime_explanation_heatmap.png")

    println("\nGenerated figures:")
    for fig in generated_figures
        println(" - $fig")
    end

    if !isempty(skipped_figures)
        println("\nSkipped optional figures:")
        for fig in skipped_figures
            println(" - $fig")
        end
    end

    println("--- Lightweight Plotting Complete ---")
end


if abspath(PROGRAM_FILE) == @__FILE__
    run_plot_dynqual_synthetic_isimip_pressure_demo()
end