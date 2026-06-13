using Test
using TwoTimescaleResilience

@testset "Movement exposure (occupancy-weighted)" begin
    # two regions, three compounds
    r1 = [1.0, 10.0, 100.0]
    r2 = [3.0, 0.0, 0.0]

    @testset "Weighted average matches hand computation" begin
        c = occupancy_weighted_exposure([r1, r2], [0.25, 0.75])
        @test c ≈ [0.25 * 1.0 + 0.75 * 3.0, 0.25 * 10.0, 0.25 * 100.0]
    end

    @testset "Single occupied region reduces to that region (resident case)" begin
        @test occupancy_weighted_exposure([r1, r2], [1.0, 0.0]) ≈ r1
        @test occupancy_weighted_exposure([r1, r2], [0.0, 1.0]) ≈ r2
    end

    @testset "Occupancy is normalised (residence-time shares, not raw weights)" begin
        # un-normalised weights give the same result as their normalised form
        @test occupancy_weighted_exposure([r1, r2], [1.0, 3.0]) ≈
              occupancy_weighted_exposure([r1, r2], [0.25, 0.75])
        # uniform occupancy = plain mean
        @test occupancy_weighted_exposure([r1, r2], [1.0, 1.0]) ≈ (r1 .+ r2) ./ 2
    end

    @testset "Matrix overload" begin
        M = [1.0 10.0 100.0; 3.0 0.0 0.0]   # regions × compounds
        @test occupancy_weighted_exposure(M, [0.25, 0.75]) ≈
              occupancy_weighted_exposure([r1, r2], [0.25, 0.75])
    end

    @testset "Guards" begin
        @test_throws ArgumentError occupancy_weighted_exposure([r1, r2], [0.0, 0.0])      # all-zero
        @test_throws ArgumentError occupancy_weighted_exposure([r1, r2], [1.0])           # length mismatch
        @test_throws ArgumentError occupancy_weighted_exposure([r1, r2], [-1.0, 2.0])     # negative weight
        @test_throws ArgumentError occupancy_weighted_exposure([r1, [1.0, 2.0]], [0.5, 0.5])  # ragged regions
        @test_throws ArgumentError occupancy_weighted_exposure([[1.0, -1.0]], [1.0])      # negative conc
    end
end

@testset "Surface:volume retention (aquatic toxicokinetics)" begin
    Lref = 26.0

    @testset "Anchored at the reference length" begin
        @test surface_volume_retention(0.7, Lref, Lref) ≈ 0.7
    end

    @testset "Smaller stage equilibrates faster; monotone in L" begin
        @test surface_volume_retention(0.7, 3.0, Lref) < 0.7        # parr: low retention, fast
        @test surface_volume_retention(0.7, 50.0, Lref) > 0.7       # larger than ref: lags more
        @test surface_volume_retention(0.7, 3.0, Lref) <
              surface_volume_retention(0.7, 10.0, Lref) <
              surface_volume_retention(0.7, Lref, Lref)
    end

    @testset "Formula rho_ref^(L_ref/L)" begin
        @test surface_volume_retention(0.7, 13.0, Lref) ≈ 0.7^(Lref / 13.0)
    end

    @testset "Gating by exposure route" begin
        @test waterborne_stage_retention(0.7, 3.0, Lref; waterborne = true) ≈
              surface_volume_retention(0.7, 3.0, Lref)
        @test waterborne_stage_retention(0.7, 3.0, Lref; waterborne = false) ≈ 0.7   # inert for non-aquatic
    end

    @testset "Guards" begin
        @test surface_volume_retention(0.0, 3.0, Lref) == 0.0
        @test_throws ArgumentError surface_volume_retention(1.0, 3.0, Lref)   # rho_ref >= 1
        @test_throws ArgumentError surface_volume_retention(-0.1, 3.0, Lref)  # rho_ref < 0
        @test_throws ArgumentError surface_volume_retention(0.7, 0.0, Lref)   # L <= 0
        @test_throws ArgumentError surface_volume_retention(0.7, 3.0, 0.0)    # L_ref <= 0
    end
end
