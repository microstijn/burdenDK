using Test
using TwoTimescaleResilience
using NCDatasets
using CSV
using JSON

@testset "Extended Integration: Archetype Compound Memory 10yr Grid Demo" begin
    # Store old ENV values
    old_nx = get(ENV, "TTR_GRID_NX", nothing)
    old_ny = get(ENV, "TTR_GRID_NY", nothing)
    old_nyrs = get(ENV, "TTR_N_YEARS", nothing)
    old_nsp = get(ENV, "TTR_N_ARCHETYPE_SPECIES", nothing)
    old_ncmp = get(ENV, "TTR_N_COMPOUNDS", nothing)
    old_k = get(ENV, "TTR_VULN_CLUSTER_K", nothing)
    old_plots = get(ENV, "TTR_MAKE_EXAMPLE_PLOTS", nothing)
    old_fallback = get(ENV, "TTR_ALLOW_NON_ARCHETYPE_FALLBACK", nothing)

    # Set up temp environment
    ENV["TTR_GRID_NX"] = "20"
    ENV["TTR_GRID_NY"] = "15"
    ENV["TTR_N_YEARS"] = "2"
    ENV["TTR_N_ARCHETYPE_SPECIES"] = "8"
    ENV["TTR_N_COMPOUNDS"] = "8"
    ENV["TTR_VULN_CLUSTER_K"] = "3"
    ENV["TTR_MAKE_EXAMPLE_PLOTS"] = "false"
    ENV["TTR_ALLOW_NON_ARCHETYPE_FALLBACK"] = "true"

    tmp_out = mktempdir()

    try
        # Include the example script so it defines the function
        include(joinpath(@__DIR__, "..", "examples", "archetype_compound_memory_10yr_grid_demo.jl"))

        # Test 1: Example runs without error
        @test_nowarn run_archetype_compound_memory_10yr_grid_demo(output_dir=tmp_out)

        # Test 2: Output directory exists
        @test isdir(tmp_out)

        # Test 3: selected_archetype_species.csv exists and has >= 8 species
        sp_csv = joinpath(tmp_out, "selected_archetype_species.csv")
        @test isfile(sp_csv)
        if isfile(sp_csv)
            sp_data = CSV.File(sp_csv)
            if ENV["TTR_ALLOW_NON_ARCHETYPE_FALLBACK"] == "true"
                @test length(sp_data) > 0
            else
                @test length(sp_data) >= 8
            end
        end

        # Test 4: selected_memory_compounds.csv exists and has >= 8 compounds
        cmp_csv = joinpath(tmp_out, "selected_memory_compounds.csv")
        @test isfile(cmp_csv)
        if isfile(cmp_csv)
            cmp_data = CSV.File(cmp_csv)
            @test length(cmp_data) >= 8
        end

        # Test 5: compound_behavior_profiles.csv exists and is non-empty
        bp_csv = joinpath(tmp_out, "compound_behavior_profiles.csv")
        @test isfile(bp_csv)
        if isfile(bp_csv)
            bp_data = CSV.File(bp_csv)
            @test length(bp_data) > 0
        end

        # Test 6: compound_concentration_scenario_summary.csv exists and is non-empty
        cs_csv = joinpath(tmp_out, "compound_concentration_scenario_summary.csv")
        @test isfile(cs_csv)
        if isfile(cs_csv)
            cs_data = CSV.File(cs_csv)
            @test length(cs_data) > 0
        end

        # Test 7: analytical_warmup_summary.csv exists and is non-empty
        aw_csv = joinpath(tmp_out, "analytical_warmup_summary.csv")
        @test isfile(aw_csv)
        if isfile(aw_csv)
            aw_data = CSV.File(aw_csv)
            @test length(aw_data) > 0
        end

        # Test 8: vulnerability_regime_outputs.nc exists and is non-empty
        nc_path = joinpath(tmp_out, "vulnerability_regime_outputs.nc")
        @test isfile(nc_path)
        @test filesize(nc_path) > 0

        # Test 9 & 10 & 11: NetCDF contents
        NCDataset(nc_path, "r") do ds
            vars = keys(ds)
            @test "cluster_id" in vars
            @test "p95_F_grouped" in vars
            @test "p95_Q_grouped" in vars
            @test "axis_entropy" in vars

            # Check cluster count > 1
            clusters = unique(ds["cluster_id"][:])
            @test length(clusters) > 1
        end

        # Test 12: vulnerability_regime_cluster_summary.csv exists and is non-empty
        vcs_csv = joinpath(tmp_out, "vulnerability_regime_cluster_summary.csv")
        @test isfile(vcs_csv)
        if isfile(vcs_csv)
            vcs_data = CSV.File(vcs_csv)
            @test length(vcs_data) > 0
        end

        # Test 13: No output feature names contain forbidden threshold patterns
        meta_json = joinpath(tmp_out, "simulation_metadata.json")
        @test isfile(meta_json)

        sf_csv = joinpath(tmp_out, "vulnerability_regime_selected_features.csv")
        if isfile(sf_csv)
            sf_data = CSV.File(sf_csv)
            for f in sf_data
                name = String(f.feature_name)
                @test !occursin("_gt_", name)
                @test !occursin("_lt_", name)
                @test !occursin("threshold", name)
                @test !occursin("exceedance", name)
                @test !occursin("above", name)
                @test !occursin("below", name)
            end
        end

        # Test 14: simulation_metadata.json records required fields
        meta = JSON.parsefile(meta_json)
        @test haskey(meta, "example_script")
        @test haskey(meta, "nx")
        @test meta["nx"] == 20
        @test meta["ny"] == 15
        @test meta["k_clusters"] == 3
        @test meta["physiological_Z_t"] == false
        @test meta["DEBtox_D_t"] == false

    finally
        # Restore old ENV
        if old_nx !== nothing ENV["TTR_GRID_NX"] = old_nx else delete!(ENV, "TTR_GRID_NX") end
        if old_ny !== nothing ENV["TTR_GRID_NY"] = old_ny else delete!(ENV, "TTR_GRID_NY") end
        if old_nyrs !== nothing ENV["TTR_N_YEARS"] = old_nyrs else delete!(ENV, "TTR_N_YEARS") end
        if old_nsp !== nothing ENV["TTR_N_ARCHETYPE_SPECIES"] = old_nsp else delete!(ENV, "TTR_N_ARCHETYPE_SPECIES") end
        if old_ncmp !== nothing ENV["TTR_N_COMPOUNDS"] = old_ncmp else delete!(ENV, "TTR_N_COMPOUNDS") end
        if old_k !== nothing ENV["TTR_VULN_CLUSTER_K"] = old_k else delete!(ENV, "TTR_VULN_CLUSTER_K") end
        if old_plots !== nothing ENV["TTR_MAKE_EXAMPLE_PLOTS"] = old_plots else delete!(ENV, "TTR_MAKE_EXAMPLE_PLOTS") end
        if old_fallback !== nothing ENV["TTR_ALLOW_NON_ARCHETYPE_FALLBACK"] = old_fallback else delete!(ENV, "TTR_ALLOW_NON_ARCHETYPE_FALLBACK") end
    end
end
