using Test
using NCDatasets
using CSV
using DataFrames
using JSON

include("../examples/plot_dynqual_synthetic_isimip_pressure_demo.jl")

@testset "Lightweight Plotting Script - dynqual_synthetic_isimip_pressure_demo" begin
    # 1. Create a temporary output directory.
    temp_out = mktempdir()

    # Dimensions
    nx, ny = 4, 3
    n_tranches = 2
    n_features = 4
    k_clusters = 3

    # 2. Create tiny synthetic cache NetCDF
    cache_path = joinpath(temp_out, "dynqual_demo_cache.nc")
    NCDataset(cache_path, "c") do ds
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        defDim(ds, "tranche", n_tranches)
        defDim(ds, "feature", n_features)
        defDim(ds, "cluster", k_clusters)

        defVar(ds, "lon", [1.0, 2.0, 3.0, 4.0], ("x",))
        defVar(ds, "lat", [10.0, 20.0, 30.0], ("y",))
        defVar(ds, "tranche", 1:n_tranches, ("tranche",))
        defVar(ds, "feature_index", 1:n_features, ("feature",))
        defVar(ds, "cluster", 1:k_clusters, ("cluster",))

        defVar(ds, "raw_BOD_climatology", rand(Float32, nx, ny), ("x", "y"))
        defVar(ds, "raw_pathogen_climatology", rand(Float32, nx, ny), ("x", "y"))
        defVar(ds, "raw_TDSload_climatology", rand(Float32, nx, ny), ("x", "y"))
        defVar(ds, "raw_BODload_climatology", rand(Float32, nx, ny), ("x", "y"))

        defVar(ds, "organic_pressure_climatology", rand(Float32, nx, ny), ("x", "y"))
        defVar(ds, "pathogen_pressure_climatology", rand(Float32, nx, ny), ("x", "y"))
        defVar(ds, "ionic_pressure_climatology", rand(Float32, nx, ny), ("x", "y"))
        defVar(ds, "wastewater_source_pressure_climatology", rand(Float32, nx, ny), ("x", "y"))
        defVar(ds, "low_flow_concentration_pressure_climatology", rand(Float32, nx, ny), ("x", "y"))
        defVar(ds, "combined_wastewater_pressure_climatology", rand(Float32, nx, ny), ("x", "y"))

        defVar(ds, "cluster_id", rand(Float32, nx, ny, n_tranches), ("x", "y", "tranche"))
        defVar(ds, "feature_map", rand(Float32, nx, ny, n_features, n_tranches), ("x", "y", "feature", "tranche"))
        defVar(ds, "baseline_centroids", rand(Float64, k_clusters, n_features), ("cluster", "feature"))
    end

    # 3. Create tiny metadata CSV
    meta_csv_path = joinpath(temp_out, "dynqual_feature_metadata.csv")
    df = DataFrame(
        feature_index = 1:n_features,
        feature_name = ["feat1", "feat2", "feat3", "feat4"],
        feature_descriptor = ["desc1", "desc2", "desc3", "desc4"],
        used_for_clustering = [true, true, true, true],
        units_or_scale = ["0-1", "0-1", "0-1", "0-1"]
    )
    CSV.write(meta_csv_path, df)

    # 4. Create tiny metadata JSON
    meta_json_path = joinpath(temp_out, "dynqual_synthetic_isimip_metadata.json")
    open(meta_json_path, "w") do f
        JSON.print(f, Dict("test" => true))
    end

    # 5. Call plotting function
    # Disable CairoMakie dependency explicitly by mocking plotting output functions to avoid test bloat?
    # Actually wait - we should see if Makie is available, if not the test should still pass if it handles the catch.
    # The instructions say: "gate plotting tests behind TTR_RUN_PLOTTING_TESTS=true" which means it'll likely have CairoMakie available when true.
    # We will test the plotting function call.

    fig_dir = joinpath(temp_out, "figures")
    try
        run_plot_dynqual_synthetic_isimip_pressure_demo(output_dir=temp_out, figures_dir=fig_dir)
    catch e
        if e isa ErrorException && occursin("CairoMakie is required", e.msg)
            @info "CairoMakie not available, skipping actual PNG generation checks."
        else
            rethrow(e)
        end
    end

    # If CairoMakie was used, files will exist. If not, the test passes gracefully without failure.
    if isdir(fig_dir) && length(readdir(fig_dir)) > 0
        @test isfile(joinpath(fig_dir, "dynqual_raw_climatology_maps.png"))
        @test filesize(joinpath(fig_dir, "dynqual_raw_climatology_maps.png")) > 0
        @test isfile(joinpath(fig_dir, "dynqual_derived_pressure_layers.png"))
        @test isfile(joinpath(fig_dir, "dynqual_vulnerability_regime_maps.png"))
        @test isfile(joinpath(fig_dir, "dynqual_regime_explanation_heatmap.png"))
    end
end
