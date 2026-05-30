using Test
using TwoTimescaleResilience
using Dates
using JSON
using CSV
using NCDatasets

@testset "Archetype Compound Memory Multi-Tranche Grid Demo" begin
    # Only run if explicitly requested
    if get(ENV, "TTR_RUN_EXTENDED_TESTS", "false") != "true" && get(ENV, "TTR_RUN_EXAMPLE_TESTS", "false") != "true"
        @info "Skipping archetype compound memory multitranche grid demo test. Set TTR_RUN_EXTENDED_TESTS=true to run."
        return
    end

    # Check for archetype database and skip clearly if missing unless fallback is enabled
    allow_fallback = get(ENV, "TTR_ALLOW_NON_ARCHETYPE_FALLBACK", "false") == "true"
    archetype_csv = joinpath(dirname(@__DIR__), "data", "AmP_Species_Archetypes.csv")
    archetype_json = joinpath(dirname(@__DIR__), "data", "AmP_Species_Archetypes.json")

    if !isfile(archetype_csv) && !isfile(archetype_json) && !allow_fallback
        @info "Skipping archetype compound memory multitranche grid demo test: Archetype database missing and fallback not explicitly enabled."
        return
    end

    # Set up deterministic environment for test execution
    ENV["TTR_GRID_NX"] = "20"
    ENV["TTR_GRID_NY"] = "15"
    ENV["TTR_N_TRANCHES"] = "3"
    ENV["TTR_TRANCHE_LENGTH_YEARS"] = "2"
    ENV["TTR_N_ARCHETYPE_SPECIES"] = "8"
    ENV["TTR_N_COMPOUNDS"] = "8"
    ENV["TTR_VULN_CLUSTER_K"] = "3"
    ENV["TTR_MAKE_EXAMPLE_PLOTS"] = "false"

    # We will modify the demo function to accept an output dir, so we can test it cleanly in a temp dir
    # To do this safely without modifying the script file during test, we'll include it inside a module

    mktempdir() do temp_dir
        # Wrap the include in a module to avoid namespace pollution
        test_mod = Module(:TestMultitrancheDemo)
        Core.eval(test_mod, :(
            using TwoTimescaleResilience, CSV, JSON, Dates, NCDatasets, Statistics;
            Base.include(TestMultitrancheDemo, joinpath($(dirname(@__DIR__)), "examples", "archetype_compound_memory_multitranche_grid_demo.jl"))
        ))

        # 1. Example runs without error
        @test_nowarn Core.eval(test_mod, :(run_archetype_compound_memory_multitranche_grid_demo(output_dir=$temp_dir)))

        # 2. Output directory exists
        @test isdir(temp_dir)

        # 3. Selected archetype species
        species_csv = joinpath(temp_dir, "selected_archetype_species.csv")
        @test isfile(species_csv)
        species_df = CSV.File(species_csv)
        @test length(species_df) >= 8

        # 4. Selected memory compounds
        compounds_csv = joinpath(temp_dir, "selected_memory_compounds.csv")
        @test isfile(compounds_csv)
        compounds_df = CSV.File(compounds_csv)
        @test length(compounds_df) >= 8

        # 5. Tranche definitions
        tranches_csv = joinpath(temp_dir, "tranche_definitions.csv")
        @test isfile(tranches_csv)
        tranches_df = CSV.File(tranches_csv)
        @test length(tranches_df) == 3

        # 6. Compound tranche trajectories
        traj_csv = joinpath(temp_dir, "compound_tranche_trajectories.csv")
        @test isfile(traj_csv)
        traj_df = CSV.File(traj_csv)
        @test "tranche_1_multiplier" in string.(propertynames(traj_df))

        # 7. Analytical warmup summary
        warmup_csv = joinpath(temp_dir, "analytical_warmup_summary.csv")
        @test isfile(warmup_csv)
        @test length(CSV.File(warmup_csv)) > 0

        # 8. Tranche memory carryover
        carryover_csv = joinpath(temp_dir, "tranche_memory_carryover_summary.csv")
        @test isfile(carryover_csv)
        @test length(CSV.File(carryover_csv)) > 0

        # 9. Feature change summary
        feature_change_csv = joinpath(temp_dir, "tranche_feature_change_summary.csv")
        @test isfile(feature_change_csv)
        @test length(CSV.File(feature_change_csv)) > 0

        # 10. Cluster transition matrix
        trans_mat_csv = joinpath(temp_dir, "tranche_cluster_transition_matrix.csv")
        @test isfile(trans_mat_csv)
        @test length(CSV.File(trans_mat_csv)) > 0

        # 11. Distribution distances
        dist_csv = joinpath(temp_dir, "tranche_cluster_distribution_distances.csv")
        @test isfile(dist_csv)
        @test length(CSV.File(dist_csv)) > 0

        # 12 & 13. NetCDF exists and contains tranche dimension
        nc_path = joinpath(temp_dir, "vulnerability_regime_multitranche_outputs.nc")
        @test isfile(nc_path)
        NCDataset(nc_path, "r") do ds
            @test haskey(ds.dim, "tranche")
            @test ds.dim["tranche"] == 3
            @test haskey(ds, "cluster_id")
            cluster_var = ds["cluster_id"]
            @test size(cluster_var) == (20, 15, 3)
        end

        # 14. Simulation metadata verification
        meta_json = joinpath(temp_dir, "simulation_metadata.json")
        @test isfile(meta_json)
        metadata = JSON.parsefile(meta_json)

        @test get(metadata, "n_tranches", -1) == 3
        @test get(metadata, "tranche_length_years", -1) == 2
        @test get(metadata, "memory_carryover", "") == "continuous_across_tranches"
        @test get(metadata, "standardisation_reference", "") == "tranche_1"
        @test get(metadata, "clustering_reference", "") == "tranche_1"
        @test get(metadata, "tranche_change_model", "") == "compound_specific_stepwise_tranche_multipliers"

        # 15. No output feature names contain forbidden patterns
        feature_names = String[]
        if haskey(metadata, "standardized_feature_names")
            append!(feature_names, metadata["standardized_feature_names"])
        end
        # Also check column names in feature change summary
        fc_df = CSV.File(feature_change_csv)
        append!(feature_names, string.(fc_df.feature_name))

        forbidden = ["_gt_", "_lt_", "threshold", "exceedance", "above", "below"]
        for fn in feature_names
            for f in forbidden
                @test !occursin(f, fn)
            end
        end

    end
end
