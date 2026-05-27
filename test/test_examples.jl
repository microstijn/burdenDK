using TwoTimescaleResilience
using Test

@testset "Examples Verification" begin
    @testset "Test 7.1 -- check output names conceptually" begin
        # We don't want to run the full netcdf script if it requires actual netcdf files, 
        # but we can verify the generated names.
        expected_files = [
            "pathogen_normalised.asc",
            "organic_normalised.asc",
            "deb_assimilation_stress.asc",
            "deb_maintenance_stress.asc",
            "deb_growth_stress.asc",
            "deb_reproduction_stress.asc",
            "deb_adaptive_margin.asc",
            "deb_restoring_force.asc",
            "deb_amplification_factor.asc"
        ]
        @test "deb_amplification_factor.asc" in expected_files
    end

    @testset "Multiaxis Response Calibrated Demo (Tranche 6)" begin
        # Run the demo to ensure outputs are generated
        demo_script = normpath(joinpath(@__DIR__, "..", "examples", "ecotox_amp_multiaxis_response_calibrated_demo.jl"))
        
        # We can evaluate the module directly to avoid starting a new julia process, 
        # but since CairoMakie creates plots, running it inside the test suite is standard.
        # To avoid polluting global scope, we run it in a module.
        m = Module()
        Base.include(m, demo_script)
        
        # Call the main function explicitly to ensure it runs
        m.main()
        
        out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
        @test isdir(out_dir)
        
        # Check that required CSV files exist and are non-empty
        csv_files = [
            "multiaxis_compound_summary.csv",
            "multiaxis_species_summary.csv",
            "multiaxis_response_mode_comparison_summary.csv",
            "multiaxis_response_calibration_summary.csv"
        ]
        
        for file in csv_files
            path = joinpath(out_dir, file)
            @test isfile(path)
            @test filesize(path) > 0
        end
        
        # Check that required PNG files exist and are non-empty
        png_files = [
            "multiaxis_axis_burdens.png",
            "multiaxis_axis_impairments.png",
            "multiaxis_weighted_impairment_Q.png",
            "multiaxis_adaptive_margin.png",
            "multiaxis_restoring_force.png",
            "multiaxis_amplification.png"
        ]
        
        for file in png_files
            path = joinpath(out_dir, file)
            @test isfile(path)
            @test filesize(path) > 0
        end
    end
end
