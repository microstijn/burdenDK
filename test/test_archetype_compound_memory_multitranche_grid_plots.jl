using Test
using NCDatasets
using CSV
using DataFrames

# Note: we use Base.include to include it within a module scope to avoid UndefVarError
module TestPlottingWrapper
    import NCDatasets
    import CSV
    import DataFrames
    # The script itself calls `using CairoMakie` etc., so they need to be available.
    Base.include(TestPlottingWrapper, joinpath(@__DIR__, "..", "examples", "plot_archetype_compound_memory_multitranche_grid_demo.jl"))
end

if get(ENV, "TTR_RUN_PLOTTING_TESTS", "false") == "true"
    @info "Running TTR_RUN_PLOTTING_TESTS tests for multi-tranche plots..."

    @testset "Multi-tranche Plotting Script Tests" begin
        # 1. Setup minimal output dir structure
        tmp_dir = mktempdir()
        out_dir = joinpath(tmp_dir, "output", "archetype_compound_memory_multitranche_grid_demo")
        mkpath(out_dir)

        # 2. Setup mock NetCDF
        nc_path = joinpath(out_dir, "vulnerability_regime_multitranche_outputs.nc")
        NCDataset(nc_path, "c") do ds
            defDim(ds, "x", 4)
            defDim(ds, "y", 3)
            defDim(ds, "tranche", 3)

            v_clust = defVar(ds, "cluster_id", Int32, ("x", "y", "tranche"))
            v_clust[:, :, :] = rand(1:5, 4, 3, 3)

            v_p95f = defVar(ds, "p95_F_grouped", Float64, ("x", "y", "tranche"))
            v_p95f[:, :, :] = rand(4, 3, 3)

            v_p95q = defVar(ds, "p95_Q_grouped", Float64, ("x", "y", "tranche"))
            v_p95q[:, :, :] = rand(4, 3, 3)

            v_ent = defVar(ds, "axis_entropy", Float64, ("x", "y", "tranche"))
            v_ent[:, :, :] = rand(4, 3, 3)

            v_ri = defVar(ds, "regime_intensity_delta_from_baseline", Float64, ("x", "y", "tranche"))
            v_ri[:, :, :] = rand(4, 3, 3)
        end

        # 3. Setup mock CSVs
        CSV.write(joinpath(out_dir, "tranche_feature_change_summary.csv"), DataFrame(
            tranche_from=[1, 1], tranche_to=[2, 3], feature_name=["p95_F_grouped", "p95_Q_grouped"], median_delta=[0.1, -0.2]
        ))

        CSV.write(joinpath(out_dir, "tranche_cluster_area_fractions.csv"), DataFrame(
            tranche=[1, 2, 3], cluster_id=[1, 2, 3], n_cells=[10, 15, 5], fraction_cells=[0.5, 0.4, 0.1]
        ))

        CSV.write(joinpath(out_dir, "tranche_cluster_transition_matrix.csv"), DataFrame(
            tranche_from=[1, 1], tranche_to=[3, 3], cluster_from=[1, 2], cluster_to=[2, 3], fraction_of_from_cluster=[0.8, 0.5]
        ))

        CSV.write(joinpath(out_dir, "tranche_cluster_area_change.csv"), DataFrame(
            tranche_from=[1, 1], tranche_to=[2, 3], cluster_id=[1, 2], delta_fraction=[-0.1, 0.2]
        ))

        CSV.write(joinpath(out_dir, "tranche_cluster_distribution_distances.csv"), DataFrame(
            tranche_from=[1, 2], tranche_to=[2, 3], total_variation_distance=[0.15, 0.25]
        ))

        ENV["TTR_MULTITRANCHE_DEMO_OUTPUT_DIR"] = out_dir
        TestPlottingWrapper.run_archetype_compound_memory_multitranche_grid_plots()

        fig_dir = joinpath(out_dir, "figures")

        @test isdir(fig_dir)
        @test isfile(joinpath(fig_dir, "cluster_maps_by_tranche.png"))
        @test filesize(joinpath(fig_dir, "cluster_maps_by_tranche.png")) > 0

        @test isfile(joinpath(fig_dir, "cluster_transition_heatmap_T1_to_final.png"))
        @test filesize(joinpath(fig_dir, "cluster_transition_heatmap_T1_to_final.png")) > 0

        @test isfile(joinpath(fig_dir, "feature_change_heatmap_from_baseline.png"))
        @test filesize(joinpath(fig_dir, "feature_change_heatmap_from_baseline.png")) > 0

        @test isfile(joinpath(fig_dir, "cluster_area_fraction_heatmap.png"))
        @test filesize(joinpath(fig_dir, "cluster_area_fraction_heatmap.png")) > 0

        @test isfile(joinpath(fig_dir, "cluster_area_delta_heatmap.png"))
        @test filesize(joinpath(fig_dir, "cluster_area_delta_heatmap.png")) > 0

        @test isfile(joinpath(fig_dir, "tranche_distance_heatmap.png"))
        @test filesize(joinpath(fig_dir, "tranche_distance_heatmap.png")) > 0

        @test isfile(joinpath(fig_dir, "regime_intensity_delta_map_final.png"))
        @test filesize(joinpath(fig_dir, "regime_intensity_delta_map_final.png")) > 0

        @test isfile(joinpath(fig_dir, "p95_F_delta_map_final.png"))
        @test filesize(joinpath(fig_dir, "p95_F_delta_map_final.png")) > 0

        @test isfile(joinpath(fig_dir, "p95_Q_delta_map_final.png"))
        @test filesize(joinpath(fig_dir, "p95_Q_delta_map_final.png")) > 0

        delete!(ENV, "TTR_MULTITRANCHE_DEMO_OUTPUT_DIR")
    end
else
    @info "Skipping multi-tranche plotting tests. Set TTR_RUN_PLOTTING_TESTS=true to run."
end
