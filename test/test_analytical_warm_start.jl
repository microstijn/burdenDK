using Test
using TwoTimescaleResilience

# Include the example script to access the helper functions.
# (Since the script has an abspath/PROGRAM_FILE check, it will not run `main()`)
include(joinpath(@__DIR__, "..", "examples", "ecotox_amp_multispecies_multicompound_monthly_memory_demo.jl"))

@testset "Analytical Warm-Start Helpers" begin
    @testset "Tranche 1: analytical_initial_burden" begin
        # explicit recurrence test
        rho = 0.9
        K = 10.0
        C_bg = 0.1
        n = 24
        B0 = 0.0

        B_recurrence = B0
        for i in 1:n
            B_recurrence = rho * B_recurrence + (1 - rho) * K * C_bg
        end

        B_analytical = analytical_initial_burden(rho, K, C_bg, n; B0 = B0)
        @test isapprox(B_recurrence, B_analytical, atol=1e-9)

        # non-zero B0 test
        B0_nz = 0.5
        B_recurrence_nz = B0_nz
        for i in 1:n
            B_recurrence_nz = rho * B_recurrence_nz + (1 - rho) * K * C_bg
        end
        B_analytical_nz = analytical_initial_burden(rho, K, C_bg, n; B0 = B0_nz)
        @test isapprox(B_recurrence_nz, B_analytical_nz, atol=1e-9)

        # C_bg = 0, B0 = 0
        @test isapprox(analytical_initial_burden(0.9, 10.0, 0.0, 24; B0 = 0.0), 0.0, atol=1e-9)

        # rho = 0, K = 1
        @test isapprox(analytical_initial_burden(0.0, 1.0, 0.5, 24; B0 = 0.0), 0.5, atol=1e-9)

        # Edge cases and ArgumentErrors
        @test_throws ArgumentError analytical_initial_burden(1.1, 10.0, 0.1, 24) # rho >= 1
        @test_throws ArgumentError analytical_initial_burden(-0.1, 10.0, 0.1, 24) # rho < 0
        @test_throws ArgumentError analytical_initial_burden(0.9, -1.0, 0.1, 24) # K <= 0 (assuming strictly > 0 for bioacc)
        @test_throws ArgumentError analytical_initial_burden(0.9, 10.0, -0.1, 24) # C_bg < 0
        @test_throws ArgumentError analytical_initial_burden(0.9, 10.0, 0.1, -1) # n < 0
    end

    @testset "Tranche 2: background_for_target_burden" begin
        rho = 0.9
        K = 10.0
        n = 24
        B0 = 0.0

        # target_B = 0.5 * NOEC test
        NOEC_dummy = 10.0
        target_B = 0.5 * NOEC_dummy

        # compute required background
        C_bg_computed = background_for_target_burden(target_B, rho, K, n; B0 = B0)

        # compute forward again using analytical function to verify recovery
        B_initial = analytical_initial_burden(rho, K, C_bg_computed, n; B0 = B0)

        @test isapprox(B_initial, target_B, atol=1e-9)
        @test C_bg_computed < (0.25 * NOEC_dummy)  # Should be lower because of bioaccumulation K

        # Verify edge cases and ArgumentErrors
        @test_throws ArgumentError background_for_target_burden(-1.0, rho, K, n) # target_B < 0
        @test_throws ArgumentError background_for_target_burden(target_B, -0.1, K, n) # rho < 0
        @test_throws ArgumentError background_for_target_burden(target_B, rho, -1.0, n) # K <= 0

        # denominator zero logic test (when rho = 1, which throws from rho validation first)
        @test_throws ArgumentError background_for_target_burden(target_B, 1.0, K, n)

        # negative C_bg expected from large initial B0
        B0_large = 1000.0
        @test_throws ArgumentError background_for_target_burden(target_B, rho, K, n; B0 = B0_large)
    end

    @testset "Tranche 3: Add zero_start and analytical_warm_start scenarios" begin
        # To test scenario execution, we can call the example script's main and capture output
        # For a clean test, we'll invoke the example workflow and check the output CSVs
        output_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multispecies_multicompound_monthly_memory_demo"))

        # run the demo script directly inside the test framework to generate outputs
        main()

        compound_csv = joinpath(output_dir, "monthly_compound_summary.csv")
        species_csv = joinpath(output_dir, "monthly_species_summary.csv")

        @test isfile(compound_csv)
        @test isfile(species_csv)

        df_comp = CSV.read(compound_csv, DataFrame)
        df_spec = CSV.read(species_csv, DataFrame)

        # Tranche 4 Row count checks (which proves scenario loop ran correctly)
        @test nrow(df_comp) == 432
        @test nrow(df_spec) == 144

        # Verify exact scenarios exist
        scenarios = unique(df_comp.scenario)
        @test sort(scenarios) == ["analytical_warm_start", "zero_start"]

        df_comp_zero = filter(r -> r.scenario == "zero_start", df_comp)
        df_comp_warm = filter(r -> r.scenario == "analytical_warm_start", df_comp)

        df_spec_zero = filter(r -> r.scenario == "zero_start", df_spec)
        df_spec_warm = filter(r -> r.scenario == "analytical_warm_start", df_spec)

        @test nrow(df_comp_zero) == 216
        @test nrow(df_comp_warm) == 216
        @test nrow(df_spec_zero) == 72
        @test nrow(df_spec_warm) == 72

        # Check months are strictly 1:12
        @test sort(unique(df_comp_zero.month)) == collect(1:12)
        @test sort(unique(df_comp_warm.month)) == collect(1:12)
        @test sort(unique(df_spec_zero.month)) == collect(1:12)
        @test sort(unique(df_spec_warm.month)) == collect(1:12)

        # Check initial values for month 1
        m1_zero = filter(r -> r.month == 1, df_comp_zero)
        m1_warm = filter(r -> r.month == 1, df_comp_warm)

        @test all(m1_zero.initial_B_t_at_reported_month_1 .== 0.0)

        # Check warm start values at month 1
        for row in eachrow(m1_warm)
            if row.chemical_name == "Sodium chloride"
                @test row.initial_B_t_at_reported_month_1 == 0.0
            else
                # For Cd and Hg, there should be a positive initial B_t since they have NOEC/EC50
                @test row.initial_B_t_at_reported_month_1 > 0.0
            end
        end

        # Check that spinup_used and method fields are populated correctly
        @test all(df_comp_zero.spinup_used .== false)
        @test all(df_comp_zero.spinup_months .== 0)
        @test all(df_comp_zero.spinup_method .== "none")

        @test all(df_comp_warm.spinup_used .== true)
        @test all(df_comp_warm.spinup_months .== 24)
        @test all(df_comp_warm.spinup_method .== "analytical_closed_form")

        # F_t approximately 1.0 for zero_start month 1 assuming C_t[1] == 0.0
        m1_spec_zero = filter(r -> r.month == 1, df_spec_zero)
        @test all(isapprox.(m1_spec_zero.F_t, 1.0, atol=1e-5))

        # Ensure all columns are finite
        @test all(isfinite.(df_comp.B_t))
        @test all(isfinite.(df_comp.x_t))
        @test all(isfinite.(df_spec.A_t))
        @test all(isfinite.(df_spec.lambda_t))
        @test all(isfinite.(df_spec.F_t))
    end
end
