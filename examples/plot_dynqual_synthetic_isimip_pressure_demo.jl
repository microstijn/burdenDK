# plot_dynqual_synthetic_isimip_pressure_demo.jl

using Statistics
using Dates
using JSON
using CSV
using DataFrames
using NCDatasets

function make_dynqual_raw_climatology_figure(fig_dir, lons, lats, raw_clim_maps)
    try
        @eval using CairoMakie
    catch
        error("CairoMakie is required for plotting but is not available.")
    end

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
    try
        @eval using CairoMakie
    catch
        error("CairoMakie is required for plotting but is not available.")
    end

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
    try
        @eval using CairoMakie
    catch
        error("CairoMakie is required for plotting but is not available.")
    end

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
    try
        @eval using CairoMakie
    catch
        error("CairoMakie is required for plotting but is not available.")
    end

    println("  -> Figure 4: Regime explanation heatmap")
    fig4 = Figure(size = (1000, 800))
    Label(fig4[0, :], "What distinguishes the vulnerability regimes?", font=:bold, fontsize=20)

    ax_heat = Axis(fig4[1, 1],
        xticks = (1:k_clusters, ["Cluster $i" for i in 1:k_clusters]),
        yticks = (1:length(baseline_kept_names), baseline_kept_names),
        xticklabelrotation = pi/4
    )

    hm_heat = heatmap!(ax_heat, 1:k_clusters, 1:length(baseline_kept_names), baseline_centroids, colormap=:coolwarm)
    Colorbar(fig4[1, 2], hm_heat, label="Standardized Feature Value")

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

    cluster_summary_csv = joinpath(output_dir, "dynqual_vulnerability_cluster_summary.csv")
    if !isfile(cluster_summary_csv)
        @warn "Cluster summary CSV missing: dynqual_vulnerability_cluster_summary.csv not found, proceeding without it."
    end

    mkpath(figures_dir)
    println("Reading from cache: $cache_path")

    # Read NetCDF
    local lons, lats, nx, ny, actual_n_tranches, k_clusters
    local mean_raw_bod, mean_raw_fc, mean_raw_tds, mean_raw_bodload
    local mean_org_p, mean_pat_p, mean_ion_p, mean_was_p, mean_low_p, mean_com_p
    local cluster_id_all, baseline_centroids_all

    NCDataset(cache_path, "r") do ds
        lons = ds["lon"][:]
        lats = ds["lat"][:]
        nx = length(lons)
        ny = length(lats)
        actual_n_tranches = length(ds["tranche"])
        k_clusters = length(ds["cluster"])

        mean_raw_bod = ds["raw_BOD_climatology"][:,:]
        mean_raw_fc = ds["raw_pathogen_climatology"][:,:]
        mean_raw_tds = ds["raw_TDSload_climatology"][:,:]
        mean_raw_bodload = ds["raw_BODload_climatology"][:,:]

        mean_org_p = ds["organic_pressure_climatology"][:,:]
        mean_pat_p = ds["pathogen_pressure_climatology"][:,:]
        mean_ion_p = ds["ionic_pressure_climatology"][:,:]
        mean_was_p = ds["wastewater_source_pressure_climatology"][:,:]
        mean_low_p = ds["low_flow_concentration_pressure_climatology"][:,:]
        mean_com_p = ds["combined_wastewater_pressure_climatology"][:,:]

        cluster_id_all = ds["cluster_id"][:, :, :]
        baseline_centroids_all = ds["baseline_centroids"][:, :]
    end

    # Read feature metadata
    df_meta = CSV.read(feature_metadata_csv, DataFrame)

    # Keep only those used for clustering
    df_clust = df_meta[df_meta.used_for_clustering .== true, :]
    baseline_kept_names = df_clust.feature_descriptor # use descriptors for better labels

    raw_clim_maps = [
        (title="BOD / organic", data=mean_raw_bod),
        (title="Pathogen / FC proxy", data=mean_raw_fc),
        (title="TDS load", data=mean_raw_tds),
        (title="BOD load", data=mean_raw_bodload)
    ]
    make_dynqual_raw_climatology_figure(figures_dir, lons, lats, raw_clim_maps)

    derived_clim_maps = [
        (title="Organic O2 Demand Proxy", data=mean_org_p),
        (title="Pathogen Exposure Proxy", data=mean_pat_p),
        (title="Ionic / Salinity Proxy", data=mean_ion_p),
        (title="Wastewater Source Proxy", data=mean_was_p),
        (title="Low-flow Concentration Proxy", data=mean_low_p),
        (title="Combined Wastewater Proxy", data=mean_com_p)
    ]
    make_dynqual_derived_pressure_layers_figure(figures_dir, lons, lats, derived_clim_maps)

    if actual_n_tranches >= 2
        make_dynqual_vulnerability_regime_maps_figure(figures_dir, lons, lats, nx, ny, k_clusters, cluster_id_all[:, :, 1], cluster_id_all[:, :, actual_n_tranches])
    end

    make_dynqual_regime_explanation_heatmap_figure(figures_dir, k_clusters, baseline_kept_names, baseline_centroids_all)

    println("--- Lightweight Plotting Complete ---")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_plot_dynqual_synthetic_isimip_pressure_demo()
end
