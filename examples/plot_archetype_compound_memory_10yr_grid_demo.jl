"""
plot_archetype_compound_memory_10yr_grid_demo.jl

Visualization script for the 10-year archetype compound-memory grid demo.
Reads NetCDF and CSV outputs and generates heatmaps and summary figures.

Run examples/archetype_compound_memory_10yr_grid_demo.jl before running this script.
"""

using Pkg

# Pre-import at the top level to avoid world age issues inside the function
try
    using CairoMakie
    using NCDatasets
    using CSV
catch e
    error("Plotting dependencies (CairoMakie, NCDatasets, CSV) are required. Please ensure they are installed.")
end

function run_archetype_compound_memory_10yr_grid_plots(;
    output_dir = get(ENV, "TTR_10YR_DEMO_OUTPUT_DIR", normpath(joinpath(@__DIR__, "..", "output", "archetype_compound_memory_10yr_grid_demo"))),
    figures_dir = joinpath(output_dir, "figures")
)
    nc_path = joinpath(output_dir, "vulnerability_regime_outputs.nc")
    if !isfile(nc_path)
        error("NetCDF output not found at $(nc_path). Run examples/archetype_compound_memory_10yr_grid_demo.jl before plotting.")
    end

    if !isdir(figures_dir)
        mkpath(figures_dir)
    end

    @info "Reading NetCDF: $nc_path"

    # Load available variables
    ds = NCDataset(nc_path, "r")
    var_names = keys(ds)

    if !("cluster_id" in var_names)
        close(ds)
        error("Required variable 'cluster_id' is missing from NetCDF.")
    end

    # Load cluster map
    cluster_map = collect(ds["cluster_id"][:, :])

    # Load optional variables
    optional_vars = [
        "p95_F_grouped",
        "p95_Q_grouped",
        "axis_entropy",
        "dominant_axis_code",
        "max_delta_F_grouped_minus_TU",
        "max_delta_F_IA_minus_TU",
        "month_of_max_p95_F_grouped"
    ]

    loaded_vars = Dict{String, Any}()
    for var in optional_vars
        if var in var_names
            loaded_vars[var] = collect(ds[var][:, :])
        else
            @warn "Optional variable '$var' is missing from NetCDF. Skipping related plots."
        end
    end

    close(ds)

    @info "Generating individual heatmaps"

    # 1. cluster_map.png
    fig_cluster = Figure(size=(800, 600))
    ax_cluster = Axis(fig_cluster[1, 1], title="Threshold-free vulnerability regimes", xlabel="x index", ylabel="y index")
    # cluster_id is categorical
    unique_clusters = unique(filter(!ismissing, cluster_map))
    n_clusters = length(unique_clusters)
    if n_clusters > 0
        hm_cluster = heatmap!(ax_cluster, cluster_map, colormap=:tab10, colorrange=(0.5, n_clusters+0.5))
        Colorbar(fig_cluster[1, 2], hm_cluster, label="cluster_id", ticks=1:n_clusters)
    else
        heatmap!(ax_cluster, cluster_map)
    end
    save(joinpath(figures_dir, "cluster_map.png"), fig_cluster)

    # 2. p95_F_grouped_map.png
    if haskey(loaded_vars, "p95_F_grouped")
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Upper-envelope amplification, grouped model: p95(F)", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, loaded_vars["p95_F_grouped"], colormap=:viridis)
        Colorbar(fig[1, 2], hm, label="p95_F_grouped")
        save(joinpath(figures_dir, "p95_F_grouped_map.png"), fig)
    end

    # 3. p95_Q_grouped_map.png
    if haskey(loaded_vars, "p95_Q_grouped")
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Upper-envelope impairment, grouped model: p95(Q)", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, loaded_vars["p95_Q_grouped"], colormap=:viridis)
        Colorbar(fig[1, 2], hm, label="p95_Q_grouped")
        save(joinpath(figures_dir, "p95_Q_grouped_map.png"), fig)
    end

    # 4. axis_entropy_map.png
    if haskey(loaded_vars, "axis_entropy")
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Axis entropy: spread of impairment across DEB-informed axes", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, loaded_vars["axis_entropy"], colormap=:viridis)
        Colorbar(fig[1, 2], hm, label="axis_entropy")
        save(joinpath(figures_dir, "axis_entropy_map.png"), fig)
    end

    # 5. mixture_sensitivity_map.png
    sens_var = ""
    if haskey(loaded_vars, "max_delta_F_grouped_minus_TU")
        sens_var = "max_delta_F_grouped_minus_TU"
    elseif haskey(loaded_vars, "max_delta_F_IA_minus_TU")
        sens_var = "max_delta_F_IA_minus_TU"
    end

    if !isempty(sens_var)
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Mixture-effect sensitivity", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, loaded_vars[sens_var], colormap=:viridis)
        Colorbar(fig[1, 2], hm, label=sens_var)
        save(joinpath(figures_dir, "mixture_sensitivity_map.png"), fig)
    end

    # 6. seasonality_map.png
    if haskey(loaded_vars, "month_of_max_p95_F_grouped")
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Month of maximum upper-envelope amplification", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, loaded_vars["month_of_max_p95_F_grouped"], colormap=:viridis)
        Colorbar(fig[1, 2], hm, label="month_of_max_p95_F_grouped")
        save(joinpath(figures_dir, "seasonality_map.png"), fig)
    end

    # 8. dominant_axis_map.png (Optional)
    if haskey(loaded_vars, "dominant_axis_code")
        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1], title="Dominant Axis Code", xlabel="x index", ylabel="y index")
        hm = heatmap!(ax, loaded_vars["dominant_axis_code"], colormap=:tab10)
        Colorbar(fig[1, 2], hm, label="dominant_axis_code")
        save(joinpath(figures_dir, "dominant_axis_map.png"), fig)
    end

    @info "Generating summary panel"

    # 7. vulnerability_regime_summary_panel.png
    panel_vars = []
    push!(panel_vars, ("cluster_id", cluster_map, :tab10))
    if haskey(loaded_vars, "p95_F_grouped")
        push!(panel_vars, ("p95_F_grouped", loaded_vars["p95_F_grouped"], :viridis))
    end
    if haskey(loaded_vars, "p95_Q_grouped")
        push!(panel_vars, ("p95_Q_grouped", loaded_vars["p95_Q_grouped"], :viridis))
    end
    if haskey(loaded_vars, "axis_entropy")
        push!(panel_vars, ("axis_entropy", loaded_vars["axis_entropy"], :viridis))
    end
    if !isempty(sens_var)
        push!(panel_vars, ("mixture sensitivity ($sens_var)", loaded_vars[sens_var], :viridis))
    end
    if haskey(loaded_vars, "month_of_max_p95_F_grouped")
        push!(panel_vars, ("month_of_max_p95_F_grouped", loaded_vars["month_of_max_p95_F_grouped"], :viridis))
    end

    n_panels = length(panel_vars)
    if n_panels > 0
        cols = ceil(Int, sqrt(n_panels))
        rows = ceil(Int, n_panels / cols)

        fig_summary = Figure(size=(cols * 400, rows * 300))

        for (i, (var_name, data, cmap)) in enumerate(panel_vars)
            r = div(i - 1, cols) + 1
            c = mod(i - 1, cols) + 1

            # Use nested grid to avoid overlap
            ga = fig_summary[r, c] = GridLayout()
            ax = Axis(ga[1, 1], title=var_name)

            if var_name == "cluster_id" && n_clusters > 0
                hm = heatmap!(ax, data, colormap=cmap, colorrange=(0.5, n_clusters+0.5))
                Colorbar(ga[1, 2], hm, ticks=1:n_clusters, label=var_name, width=10)
            else
                hm = heatmap!(ax, data, colormap=cmap)
                Colorbar(ga[1, 2], hm, label=var_name, width=10)
            end
        end
        save(joinpath(figures_dir, "vulnerability_regime_summary_panel.png"), fig_summary)
    else
        @warn "Not enough variables for summary panel. Summary panel creation skipped."
    end

    @info "Generating CSV-based plots"

    # 10. cluster_size_barplot.png
    cluster_csv_path = joinpath(output_dir, "vulnerability_regime_cluster_summary.csv")
    if isfile(cluster_csv_path)
        cluster_df = CSV.File(cluster_csv_path)
        if :cluster_id in propertynames(cluster_df) && :n_cells in propertynames(cluster_df)
            fig = Figure(size=(800, 600))
            ax = Axis(fig[1, 1], title="Cells per Vulnerability Regime Cluster", xlabel="Cluster ID", ylabel="Number of Cells")

            x_vals = Float32.(cluster_df.cluster_id)
            y_vals = Float32.(cluster_df.n_cells)

            barplot!(ax, x_vals, y_vals, color=:steelblue)

            if :suggested_regime_label in propertynames(cluster_df)
                ax.xticks = (x_vals, string.(cluster_df.suggested_regime_label))
                ax.xticklabelrotation = pi/4
            else
                ax.xticks = (x_vals, string.(cluster_df.cluster_id))
            end
            save(joinpath(figures_dir, "cluster_size_barplot.png"), fig)
        else
            @warn "Expected columns missing from $cluster_csv_path"
        end
    else
        @warn "File not found: $cluster_csv_path. Skipping cluster_size_barplot.png"
    end

    # 9. compound_behavior_profiles_barplot.png
    compound_csv_path = joinpath(output_dir, "compound_behavior_profiles.csv")
    if isfile(compound_csv_path)
        compounds_df = CSV.File(compound_csv_path)
        if :behavior_profile in propertynames(compounds_df)
            # Count profiles
            profiles = string.(compounds_df.behavior_profile)
            unique_profiles = unique(profiles)
            counts = [count(x -> x == p, profiles) for p in unique_profiles]

            fig = Figure(size=(800, 600))
            ax = Axis(fig[1, 1], title="Compounds per Behavior Profile", xlabel="Behavior Profile", ylabel="Number of Compounds")

            x_vals = Float32.(1:length(unique_profiles))
            y_vals = Float32.(counts)

            barplot!(ax, x_vals, y_vals, color=:darkorange)

            ax.xticks = (x_vals, unique_profiles)
            ax.xticklabelrotation = pi/4
            save(joinpath(figures_dir, "compound_behavior_profiles_barplot.png"), fig)
        else
            @warn "Expected column 'behavior_profile' missing from $compound_csv_path"
        end
    else
        @warn "File not found: $compound_csv_path. Skipping compound_behavior_profiles_barplot.png"
    end

    @info "Plotting complete. Figures saved to $figures_dir"
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_archetype_compound_memory_10yr_grid_plots()
end