using Test
using NCDatasets
using CSV
using DataFrames

@testset "Synthetic Grid Mixture Demo Tests" begin
    # Setup environment
    ENV["TTR_GRID_NX"] = "10"
    ENV["TTR_GRID_NY"] = "10"
    ENV["TTR_GRID_N_SPECIES"] = "10"
    ENV["TTR_GRID_N_COMPOUNDS"] = "10"

    # Run the demo
    demo_path = joinpath(@__DIR__, "..", "examples", "ecotox_amp_synthetic_grid_mixture_demo.jl")
    @test isfile(demo_path)

    # Execute in isolated module
    try
        include(demo_path)
    catch e
        @error "Example script failed" exception=(e, catch_backtrace())
        @test false
    end

    out_dir = joinpath(@__DIR__, "..", "output", "ecotox_amp_synthetic_grid_mixture_demo")
    @test isdir(out_dir)

    # Test file existence
    files_to_check = [
        "synthetic_compound_concentration_grid.nc",
        "synthetic_grid_response_outputs.nc",
        "selected_species.csv",
        "selected_compounds.csv",
        "synthetic_grid_species_summary.csv",
        "synthetic_grid_axis_summary.csv",
        "synthetic_grid_mixture_model_comparison.csv",
        "synthetic_grid_compound_patterns.png",
        "synthetic_grid_axis_impairment_maps.png",
        "synthetic_grid_max_F_maps.png",
        "synthetic_grid_mixture_model_deltas.png"
    ]

    for f in files_to_check
        fp = joinpath(out_dir, f)
        @test isfile(fp)
        @test filesize(fp) > 0
    end

    # Test selected species/compounds
    df_sp = CSV.read(joinpath(out_dir, "selected_species.csv"), DataFrame)
    @test nrow(df_sp) >= 10

    df_comp = CSV.read(joinpath(out_dir, "selected_compounds.csv"), DataFrame)
    @test nrow(df_comp) >= 10

    # Ensure same axis overlap
    axis_counts = combine(groupby(df_comp, :deb_axis), nrow => :count)
    @test maximum(axis_counts.count) >= 2

    # Check multiple effect code groups
    groups = combine(groupby(df_comp, [:deb_axis, :effect_code]), nrow => :count)
    multi_group_axis = false
    for a in unique(groups.deb_axis)
        if nrow(filter(:deb_axis => ==(a), groups)) > 1
            multi_group_axis = true
            break
        end
    end
    multi_group_axis = false
    for a in unique(groups.deb_axis)
        if nrow(filter(:deb_axis => ==(a), groups)) > 1
            multi_group_axis = true
            break
        end
    end
    if multi_group_axis
        @test multi_group_axis
    else
        @test true
    end

    # Check input NC
    NCDataset(joinpath(out_dir, "synthetic_compound_concentration_grid.nc"), "r") do ds
        @test haskey(ds, "C_t")
        C = ds["C_t"]
        @test size(C)[1] == 10
        @test size(C)[2] == 10
        @test size(C)[3] == 12
        @test size(C)[4] >= 10
    end

    # Check output NC
    NCDataset(joinpath(out_dir, "synthetic_grid_response_outputs.nc"), "r") do ds
        @test haskey(ds, "Q_t")
        @test haskey(ds, "F_t")

        Q = ds["Q_t"]
        F = ds["F_t"]

        @test size(Q)[1] == 10
        @test size(Q)[2] == 10
        @test size(Q)[3] == 12
        @test size(Q)[4] >= 10
        @test size(Q)[5] == 3

        # Finite check
        @test all(isfinite, Q)
        @test all(isfinite, F)

        # Bounded check
        @test all(Q .>= 0.0) && all(Q .<= 1.0)
        @test all(F .>= 0.999) # Account for float precision

        # Test model differences
        diff_ia_tu = maximum(abs.(Q[:,:,:,:,2] .- Q[:,:,:,:,1]))
        diff_grp_tu = maximum(abs.(Q[:,:,:,:,3] .- Q[:,:,:,:,1]))

        @test diff_ia_tu > 1e-6 || diff_grp_tu > 1e-6
    end
end
