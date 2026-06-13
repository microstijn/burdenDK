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
