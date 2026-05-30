"""
plot_archetype_compound_memory_multitranche_grid_demo.jl

Visualization script for the multi-tranche archetype compound-memory grid demo.
Reads NetCDF and CSV outputs and generates heatmaps and summary figures.

Run examples/archetype_compound_memory_multitranche_grid_demo.jl before running this script.
"""

try
    using CairoMakie
    using NCDatasets
    using CSV
    using DataFrames
catch e
    error("Plotting dependencies (CairoMakie, NCDatasets, CSV, DataFrames) are required. Please ensure they are installed.")
end

function run_archetype_compound_memory_multitranche_grid_plots(;
    output_dir = get(ENV, "TTR_MULTITRANCHE_DEMO_OUTPUT_DIR", normpath(joinpath(@__DIR__, "..", "output", "archetype_compound_memory_multitranche_grid_demo"))),
    figures_dir = joinpath(output_dir, "figures")
)
    nc_path = joinpath(output_dir, "vulnerability_regime_multitranche_outputs.nc")
    if !isfile(nc_path)
        error("NetCDF output not found at $(nc_path). Run examples/archetype_compound_memory_multitranche_grid_demo.jl before plotting.")
    end

    if !isdir(figures_dir)
        mkpath(figures_dir)
    end

    @info "Reading NetCDF: $nc_path"
    generated_figures = String[]
    skipped_figures = String[]

    # Load available variables
    ds = NCDataset(nc_path, "r")
    var_names = keys(ds)

    if !("cluster_id" in var_names)
        close(ds)
        error("Required variable 'cluster_id' is missing from NetCDF.")
    end

    # Dimension order: [x, y, tranche]
    # cluster_id[x, y, tranche]
    cluster_map = collect(ds["cluster_id"][:, :, :])
    nx, ny, n_tranches = size(cluster_map)

    # Load optional variables
    optional_vars = [
        "p95_F_grouped",
        "p95_Q_grouped",
        "max_F_grouped",
        "max_Q_grouped",
        "axis_entropy",
        "max_delta_F_grouped_minus_TU",
        "max_delta_F_IA_minus_TU",
        "month_of_max_p95_F_grouped",
        "cluster_changed_from_baseline",
        "regime_intensity_delta_from_baseline"
    ]

    loaded_vars = Dict{String, Any}()
    for var in optional_vars
        if var in var_names
            loaded_vars[var] = collect(ds[var][:, :, :])
        else
            @warn "Optional variable '$var' is missing from NetCDF."
        end
    end

    close(ds)

    # Ensure plotting functions are fully self-contained as helper methods or part of main logic.
    # FIGURE 1: CLUSTER MAPS BY TRANCHE
    @info "Generating Figure 1: cluster_maps_by_tranche.png"
    ncols = min(n_tranches, 4)
    nrows = ceil(Int, n_tranches / ncols)

    fig_clusters = Figure(size=(ncols * 400, nrows * 300))
    Label(fig_clusters[0, :], "Threshold-free vulnerability regimes by tranche", fontsize=20, font=:bold)

    unique_clusters = unique(filter(!ismissing, cluster_map))
    max_k = isempty(unique_clusters) ? 0 : maximum(unique_clusters)

    global_hm = nothing
    for t in 1:n_tranches
        r = div(t - 1, ncols) + 1
        c = mod(t - 1, ncols) + 1
        ax = Axis(fig_clusters[r, c], title="Tranche $t", xlabel="x index", ylabel="y index")

        if max_k > 0
            hm = heatmap!(ax, cluster_map[:, :, t], colormap=:tab10, colorrange=(0.5, max_k + 0.5))
            if isnothing(global_hm)
                global_hm = hm
            end
        else
            heatmap!(ax, cluster_map[:, :, t])
        end
    end

    if !isnothing(global_hm) && max_k > 0
        Colorbar(fig_clusters[1:nrows, ncols+1], global_hm, label="cluster_id", ticks=1:max_k)
    end

    cluster_maps_path = joinpath(figures_dir, "cluster_maps_by_tranche.png")
    save(cluster_maps_path, fig_clusters)
    push!(generated_figures, "cluster_maps_by_tranche.png")

    # FIGURE 2: CLUSTER TRANSITION HEATMAP
    trans_csv = joinpath(output_dir, "tranche_cluster_transition_matrix.csv")
    if isfile(trans_csv)
        @info "Generating Figure 2: cluster_transition_heatmap_T1_to_final.png"
        df_trans = DataFrame(CSV.File(trans_csv))

        # Determine K
        k_trans = 0
        if "cluster_from" in names(df_trans) && "cluster_to" in names(df_trans)
            k_trans = max(maximum(df_trans.cluster_from, init=0), maximum(df_trans.cluster_to, init=0))
        end
        K = max(max_k, k_trans)

        if K > 0 && "tranche_from" in names(df_trans) && "tranche_to" in names(df_trans) && "fraction_of_from_cluster" in names(df_trans)
            df_subset = filter(row -> row.tranche_from == 1 && row.tranche_to == n_tranches, df_trans)

            trans_matrix = zeros(Float64, K, K)
            for row in eachrow(df_subset)
                cf = row.cluster_from
                ct = row.cluster_to
                if 1 <= cf <= K && 1 <= ct <= K
                    trans_matrix[cf, ct] = row.fraction_of_from_cluster
                end
            end

            fig_trans = Figure(size=(600, 500))
            ax_trans = Axis(fig_trans[1, 1], title="Cluster transition matrix: Tranche 1 to final tranche",
                            xlabel="target cluster", ylabel="source cluster",
                            xticks=(1:K, string.(1:K)), yticks=(1:K, string.(1:K)))

            # Note: transpose to make rows = source (y-axis) and columns = target (x-axis)
            # Makie's heatmap plots matrix M[x, y] with x on x-axis and y on y-axis.
            # So if we want rows=source, columns=target, target should be x and source should be y.
            # Thus, we pass trans_matrix' (transpose) to heatmap!
            hm_trans = heatmap!(ax_trans, trans_matrix', colormap=:viridis)
            Colorbar(fig_trans[1, 2], hm_trans, label="fraction of source cluster")

            trans_path = joinpath(figures_dir, "cluster_transition_heatmap_T1_to_final.png")
            save(trans_path, fig_trans)
            push!(generated_figures, "cluster_transition_heatmap_T1_to_final.png")
        else
            @warn "Required columns missing in $trans_csv or K=0. Skipping Figure 2."
            push!(skipped_figures, "cluster_transition_heatmap_T1_to_final.png")
        end
    else
        @warn "File not found: $trans_csv. Skipping Figure 2."
        push!(skipped_figures, "cluster_transition_heatmap_T1_to_final.png")
    end

    # FIGURE 3: FEATURE-CHANGE HEATMAP FROM BASELINE
    feat_csv = joinpath(output_dir, "tranche_feature_change_summary.csv")
    if isfile(feat_csv)
        @info "Generating Figure 3: feature_change_heatmap_from_baseline.png"
        df_feat = DataFrame(CSV.File(feat_csv))

        required_cols = ["tranche_from", "tranche_to", "feature_name"]
        if all(c -> c in names(df_feat), required_cols) && ("median_delta" in names(df_feat) || "mean_delta" in names(df_feat))
            val_col = "median_delta" in names(df_feat) ? :median_delta : :mean_delta

            # Filter from baseline
            df_feat_base = filter(row -> row.tranche_from == 1, df_feat)

            # Exclude month-related variables
            filter!(row -> !occursin("month_of_max", row.feature_name), df_feat_base)

            unique_features = unique(df_feat_base.feature_name)
            unique_targets = sort(unique(df_feat_base.tranche_to))

            if !isempty(unique_features) && !isempty(unique_targets)
                feat_matrix = zeros(Float64, length(unique_targets), length(unique_features))

                for (i, t_to) in enumerate(unique_targets)
                    for (j, feat) in enumerate(unique_features)
                        r = filter(row -> row.tranche_to == t_to && row.feature_name == feat, df_feat_base)
                        if nrow(r) > 0
                            feat_matrix[i, j] = r[1, val_col]
                        end
                    end
                end

                max_abs = maximum(abs.(feat_matrix))
                crange = max_abs > 0 ? (-max_abs, max_abs) : (-1.0, 1.0)

                fig_feat = Figure(size=(800, 600))
                target_labels = ["T1→T$t" for t in unique_targets]
                ax_feat = Axis(fig_feat[1, 1], title="Median feature change from baseline tranche (raw deltas)",
                               xticks=(1:length(unique_targets), target_labels),
                               yticks=(1:length(unique_features), unique_features))

                # Makie plots [x, y], so x=targets, y=features -> feat_matrix
                hm_feat = heatmap!(ax_feat, feat_matrix, colormap=:balance, colorrange=crange)
                Colorbar(fig_feat[1, 2], hm_feat, label="median Δ feature")

                feat_path = joinpath(figures_dir, "feature_change_heatmap_from_baseline.png")
                save(feat_path, fig_feat)
                push!(generated_figures, "feature_change_heatmap_from_baseline.png")
            else
                @warn "No valid features or targets found in $feat_csv. Skipping Figure 3."
                push!(skipped_figures, "feature_change_heatmap_from_baseline.png")
            end
        else
            @warn "Required columns missing in $feat_csv. Skipping Figure 3."
            push!(skipped_figures, "feature_change_heatmap_from_baseline.png")
        end
    else
        @warn "File not found: $feat_csv. Skipping Figure 3."
        push!(skipped_figures, "feature_change_heatmap_from_baseline.png")
    end

    # FIGURE 4: CLUSTER AREA FRACTION HEATMAP
    area_csv = joinpath(output_dir, "tranche_cluster_area_fractions.csv")
    if isfile(area_csv)
        @info "Generating Figure 4: cluster_area_fraction_heatmap.png"
        df_area = DataFrame(CSV.File(area_csv))

        if "tranche" in names(df_area) && "cluster_id" in names(df_area) && "fraction_cells" in names(df_area)
            c_ids = sort(unique(df_area.cluster_id))
            t_ids = sort(unique(df_area.tranche))

            if !isempty(c_ids) && !isempty(t_ids)
                area_matrix = zeros(Float64, length(t_ids), length(c_ids))

                for (i, t) in enumerate(t_ids)
                    for (j, c) in enumerate(c_ids)
                        r = filter(row -> row.tranche == t && row.cluster_id == c, df_area)
                        if nrow(r) > 0
                            area_matrix[i, j] = r[1, :fraction_cells]
                        end
                    end
                end

                fig_area = Figure(size=(600, 500))
                ax_area = Axis(fig_area[1, 1], title="Cluster area fractions by tranche",
                               xticks=(1:length(t_ids), string.(t_ids)), xlabel="tranche",
                               yticks=(1:length(c_ids), string.(c_ids)), ylabel="cluster_id")

                hm_area = heatmap!(ax_area, area_matrix, colormap=:viridis)
                Colorbar(fig_area[1, 2], hm_area, label="fraction of cells")

                area_path = joinpath(figures_dir, "cluster_area_fraction_heatmap.png")
                save(area_path, fig_area)
                push!(generated_figures, "cluster_area_fraction_heatmap.png")
            else
                @warn "No valid clusters or tranches in $area_csv. Skipping Figure 4."
                push!(skipped_figures, "cluster_area_fraction_heatmap.png")
            end
        else
            @warn "Required columns missing in $area_csv. Skipping Figure 4."
            push!(skipped_figures, "cluster_area_fraction_heatmap.png")
        end
    else
        @warn "File not found: $area_csv. Skipping Figure 4."
        push!(skipped_figures, "cluster_area_fraction_heatmap.png")
    end

    # FIGURE 5: CLUSTER AREA DELTA HEATMAP
    area_delta_csv = joinpath(output_dir, "tranche_cluster_area_change.csv")
    if isfile(area_delta_csv)
        @info "Generating Figure 5: cluster_area_delta_heatmap.png"
        df_area_delta = DataFrame(CSV.File(area_delta_csv))

        if "tranche_from" in names(df_area_delta) && "tranche_to" in names(df_area_delta) && "cluster_id" in names(df_area_delta) && "delta_fraction" in names(df_area_delta)
            df_ad_base = filter(row -> row.tranche_from == 1, df_area_delta)
            c_ids = sort(unique(df_ad_base.cluster_id))
            t_tos = sort(unique(df_ad_base.tranche_to))

            if !isempty(c_ids) && !isempty(t_tos)
                ad_matrix = zeros(Float64, length(t_tos), length(c_ids))

                for (i, t) in enumerate(t_tos)
                    for (j, c) in enumerate(c_ids)
                        r = filter(row -> row.tranche_to == t && row.cluster_id == c, df_ad_base)
                        if nrow(r) > 0
                            ad_matrix[i, j] = r[1, :delta_fraction]
                        end
                    end
                end

                max_ad = maximum(abs.(ad_matrix))
                crange = max_ad > 0 ? (-max_ad, max_ad) : (-1.0, 1.0)

                fig_ad = Figure(size=(600, 500))
                ax_ad = Axis(fig_ad[1, 1], title="Change in cluster area fraction from baseline",
                               xticks=(1:length(t_tos), ["T1→T$t" for t in t_tos]),
                               yticks=(1:length(c_ids), string.(c_ids)), ylabel="cluster_id")

                hm_ad = heatmap!(ax_ad, ad_matrix, colormap=:balance, colorrange=crange)
                Colorbar(fig_ad[1, 2], hm_ad, label="Δ fraction of cells")

                ad_path = joinpath(figures_dir, "cluster_area_delta_heatmap.png")
                save(ad_path, fig_ad)
                push!(generated_figures, "cluster_area_delta_heatmap.png")
            else
                @warn "No valid clusters or comparisons in $area_delta_csv. Skipping Figure 5."
                push!(skipped_figures, "cluster_area_delta_heatmap.png")
            end
        else
            @warn "Required columns missing in $area_delta_csv. Skipping Figure 5."
            push!(skipped_figures, "cluster_area_delta_heatmap.png")
        end
    else
        @warn "File not found: $area_delta_csv. Skipping Figure 5."
        push!(skipped_figures, "cluster_area_delta_heatmap.png")
    end

    # FIGURE 6: TRANCHE DISTANCE HEATMAP
    dist_csv = joinpath(output_dir, "tranche_cluster_distribution_distances.csv")
    if isfile(dist_csv)
        @info "Generating Figure 6: tranche_distance_heatmap.png"
        df_dist = DataFrame(CSV.File(dist_csv))

        if "tranche_from" in names(df_dist) && "tranche_to" in names(df_dist)
            val_col = "total_variation_distance" in names(df_dist) ? :total_variation_distance :
                      ("jensen_shannon_divergence" in names(df_dist) ? :jensen_shannon_divergence : nothing)

            if !isnothing(val_col)
                ts = sort(unique(vcat(df_dist.tranche_from, df_dist.tranche_to)))
                dist_matrix = zeros(Float64, length(ts), length(ts))

                for row in eachrow(df_dist)
                    i = findfirst(x -> x == row.tranche_from, ts)
                    j = findfirst(x -> x == row.tranche_to, ts)
                    if !isnothing(i) && !isnothing(j)
                        dist_matrix[i, j] = row[val_col]
                    end
                end

                # Make symmetric if it only contains 1-way
                for i in 1:length(ts)
                    for j in i+1:length(ts)
                        if dist_matrix[i, j] > 0 && dist_matrix[j, i] == 0
                            dist_matrix[j, i] = dist_matrix[i, j]
                        elseif dist_matrix[j, i] > 0 && dist_matrix[i, j] == 0
                            dist_matrix[i, j] = dist_matrix[j, i]
                        end
                    end
                end

                fig_dist = Figure(size=(600, 500))
                ax_dist = Axis(fig_dist[1, 1], title="Cluster-distribution distance between tranches",
                               xticks=(1:length(ts), string.(ts)), xlabel="tranche",
                               yticks=(1:length(ts), string.(ts)), ylabel="tranche")

                hm_dist = heatmap!(ax_dist, dist_matrix, colormap=:viridis)
                Colorbar(fig_dist[1, 2], hm_dist, label=string(val_col))

                dist_path = joinpath(figures_dir, "tranche_distance_heatmap.png")
                save(dist_path, fig_dist)
                push!(generated_figures, "tranche_distance_heatmap.png")
            else
                @warn "Distance columns missing in $dist_csv. Skipping Figure 6."
                push!(skipped_figures, "tranche_distance_heatmap.png")
            end
        else
            @warn "tranche_from or tranche_to missing in $dist_csv. Skipping Figure 6."
            push!(skipped_figures, "tranche_distance_heatmap.png")
        end
    else
        @warn "File not found: $dist_csv. Skipping Figure 6."
        push!(skipped_figures, "tranche_distance_heatmap.png")
    end

    # FIGURE 7: REGIME INTENSITY DELTA MAP FINAL
    if haskey(loaded_vars, "regime_intensity_delta_from_baseline")
        @info "Generating Figure 7: regime_intensity_delta_map_final.png"
        data = loaded_vars["regime_intensity_delta_from_baseline"][:, :, end]
        max_v = maximum(abs.(filter(x -> !ismissing(x) && isfinite(x), data)))
        crange = max_v > 0 ? (-max_v, max_v) : (-1.0, 1.0)

        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Regime-intensity change from baseline: final tranche", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, data, colormap=:balance, colorrange=crange)
        Colorbar(fig[1, 2], hm, label="Δ relative regime intensity")

        save(joinpath(figures_dir, "regime_intensity_delta_map_final.png"), fig)
        push!(generated_figures, "regime_intensity_delta_map_final.png")
    else
        @warn "Variable 'regime_intensity_delta_from_baseline' missing. Skipping Figure 7."
        push!(skipped_figures, "regime_intensity_delta_map_final.png")
    end

    # FIGURE 8: P95 F DELTA MAP FINAL
    if haskey(loaded_vars, "p95_F_grouped")
        @info "Generating Figure 8: p95_F_delta_map_final.png"
        data = loaded_vars["p95_F_grouped"][:, :, end] .- loaded_vars["p95_F_grouped"][:, :, 1]
        max_v = maximum(abs.(filter(x -> !ismissing(x) && isfinite(x), data)))
        crange = max_v > 0 ? (-max_v, max_v) : (-1.0, 1.0)

        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Change in p95(F), grouped model: final minus baseline", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, data, colormap=:balance, colorrange=crange)
        Colorbar(fig[1, 2], hm, label="Δ p95(F)")

        save(joinpath(figures_dir, "p95_F_delta_map_final.png"), fig)
        push!(generated_figures, "p95_F_delta_map_final.png")
    else
        @warn "Variable 'p95_F_grouped' missing. Skipping Figure 8."
        push!(skipped_figures, "p95_F_delta_map_final.png")
    end

    # FIGURE 9: P95 Q DELTA MAP FINAL
    if haskey(loaded_vars, "p95_Q_grouped")
        @info "Generating Figure 9: p95_Q_delta_map_final.png"
        data = loaded_vars["p95_Q_grouped"][:, :, end] .- loaded_vars["p95_Q_grouped"][:, :, 1]
        max_v = maximum(abs.(filter(x -> !ismissing(x) && isfinite(x), data)))
        crange = max_v > 0 ? (-max_v, max_v) : (-1.0, 1.0)

        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Change in p95(Q), grouped model: final minus baseline", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, data, colormap=:balance, colorrange=crange)
        Colorbar(fig[1, 2], hm, label="Δ p95(Q)")

        save(joinpath(figures_dir, "p95_Q_delta_map_final.png"), fig)
        push!(generated_figures, "p95_Q_delta_map_final.png")
    else
        @warn "Variable 'p95_Q_grouped' missing. Skipping Figure 9."
        push!(skipped_figures, "p95_Q_delta_map_final.png")
    end

    @info "=== SUMMARY ==="
    @info "Generated figures: " * join(generated_figures, ", ")
    if !isempty(skipped_figures)
        @info "Skipped figures: " * join(skipped_figures, ", ")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_archetype_compound_memory_multitranche_grid_plots()
end
