using Test
using TwoTimescaleResilience

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

        # zero spinup returns B0
        @test isapprox(analytical_initial_burden(0.5, 2.0, 5.0, 0; B0 = 3.0), 3.0, atol=1e-9)

        # rho = 0, gives K*C_bg for spinup > 0
        @test isapprox(analytical_initial_burden(0.0, 2.0, 5.0, 1; B0 = 0.0), 10.0, atol=1e-9)
        @test isapprox(analytical_initial_burden(0.0, 2.0, 5.0, 12; B0 = 0.0), 10.0, atol=1e-9)

        # finite spinup examples
        @test isapprox(analytical_initial_burden(0.5, 2.0, 5.0, 1; B0 = 0.0), 5.0, atol=1e-9)
        @test isapprox(analytical_initial_burden(0.5, 2.0, 5.0, 2; B0 = 0.0), 7.5, atol=1e-9)

        # long spinup approaches K*C_bg
        @test isapprox(analytical_initial_burden(0.9, 10.0, 2.0, 10000; B0 = 0.0), 20.0, atol=1e-9)

        # Edge cases and ArgumentErrors
        @test_throws ArgumentError analytical_initial_burden(1.1, 10.0, 0.1, 24) # rho >= 1
        @test_throws ArgumentError analytical_initial_burden(-0.1, 10.0, 0.1, 24) # rho < 0
        @test_throws ArgumentError analytical_initial_burden(0.9, -1.0, 0.1, 24) # K <= 0 (assuming strictly > 0 for bioacc)
        @test_throws ArgumentError analytical_initial_burden(0.9, 10.0, -0.1, 24) # C_bg < 0
        @test_throws ArgumentError analytical_initial_burden(0.9, 10.0, 0.1, -1) # n < 0
        @test_throws ArgumentError analytical_initial_burden(NaN, 10.0, 0.1, 24) # rho NaN
        @test_throws ArgumentError analytical_initial_burden(0.9, Inf, 0.1, 24) # K Inf
        @test_throws ArgumentError analytical_initial_burden(0.9, 10.0, NaN, 24) # C_bg NaN
        @test_throws ArgumentError analytical_initial_burden(0.9, 10.0, 0.1, 24, B0 = -1.0) # B0 < 0
        @test_throws ArgumentError analytical_initial_burden(0.9, 10.0, 0.1, 24, B0 = Inf) # B0 Inf
    end

    @testset "Tranche 2: background_for_target_burden" begin
        rho = 0.9
        K = 10.0
        n = 24
        B0 = 0.0

        # target_B = 0.5 * NOEC test
        NOEC_dummy = 10.0
        target_B = 0.5 * NOEC_dummy

        # roundtrip with analytical_initial_burden
        C_bg_computed = background_for_target_burden(rho, K, target_B, n; B0 = B0)
        B_initial = analytical_initial_burden(rho, K, C_bg_computed, n; B0 = B0)
        @test isapprox(B_initial, target_B, atol=1e-9)
        @test C_bg_computed < (0.25 * NOEC_dummy)  # Should be lower because of bioaccumulation K

        # rho = 0 gives B_target / K
        @test isapprox(background_for_target_burden(0.0, 2.0, 10.0, n; B0 = B0), 5.0, atol=1e-9)

        # Verify edge cases and ArgumentErrors
        @test_throws ArgumentError background_for_target_burden(rho, K, -1.0, n) # target_B < 0
        @test_throws ArgumentError background_for_target_burden(-0.1, K, target_B, n) # rho < 0
        @test_throws ArgumentError background_for_target_burden(rho, -1.0, target_B, n) # K <= 0

        # spinup_months == 0 throws ArgumentError
        @test_throws ArgumentError background_for_target_burden(rho, K, target_B, 0)

        # unreachable target requiring negative C_bg throws ArgumentError
        B0_large = 1000.0
        @test_throws ArgumentError background_for_target_burden(rho, K, target_B, n; B0 = B0_large)
    end

    @testset "Tranche 3: analytical_periodic_initial_burden" begin
        # constant cycle equals K*C
        C_cycle_const = fill(5.0, 12)
        @test isapprox(analytical_periodic_initial_burden(0.9, 2.0, C_cycle_const), 10.0, atol=1e-9)

        # rho = 0 returns K*C_cycle[end]
        C_cycle_rho_zero = [1.0, 2.0, 4.0]
        @test isapprox(analytical_periodic_initial_burden(0.0, 3.0, C_cycle_rho_zero), 12.0, atol=1e-9)

        # one-cycle periodic consistency
        rho = 0.8
        K = 3.0
        C_cycle_periodic = [1.0, 2.0, 4.0, 8.0]
        B0_periodic = analytical_periodic_initial_burden(rho, K, C_cycle_periodic)

        B_curr = B0_periodic
        for c in C_cycle_periodic
            B_curr = rho * B_curr + (1.0 - rho) * K * c
        end
        @test isapprox(B_curr, B0_periodic, atol=1e-9)

        # invalid inputs throw ArgumentError
        @test_throws ArgumentError analytical_periodic_initial_burden(-0.1, 2.0, C_cycle_const) # rho < 0
        @test_throws ArgumentError analytical_periodic_initial_burden(1.1, 2.0, C_cycle_const) # rho >= 1
        @test_throws ArgumentError analytical_periodic_initial_burden(NaN, 2.0, C_cycle_const) # rho NaN
        @test_throws ArgumentError analytical_periodic_initial_burden(0.9, -1.0, C_cycle_const) # K <= 0
        @test_throws ArgumentError analytical_periodic_initial_burden(0.9, Inf, C_cycle_const) # K Inf
        @test_throws ArgumentError analytical_periodic_initial_burden(0.9, 2.0, Float64[]) # empty C_cycle
        @test_throws ArgumentError analytical_periodic_initial_burden(0.9, 2.0, [1.0, -2.0]) # C_cycle containing negative value
        @test_throws ArgumentError analytical_periodic_initial_burden(0.9, 2.0, [1.0, NaN]) # C_cycle containing NaN
        @test_throws ArgumentError analytical_periodic_initial_burden(0.9, 2.0, [1.0, Inf]) # C_cycle containing Inf
    end
end
