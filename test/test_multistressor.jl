using Test
using TwoTimescaleResilience

@testset "Tranche 2: Multi-Stressor Background Burden" begin
    # Test 2.1 -- additive burden
    stressors = [
        BackgroundStressor("a", 0.2, 0.5),
        BackgroundStressor("b", 0.4, 0.25),
        BackgroundStressor("c", 0.8, 0.25)
    ]
    # Expected: B = 0.5(0.2) + 0.25(0.4) + 0.25(0.8) = 0.1 + 0.1 + 0.2 = 0.4
    @test isapprox(background_index(stressors), 0.4; atol=1e-10)

    # Test 2.2 -- pairwise interaction
    x = [2.0, 3.0, 4.0]
    M = zeros(3, 3)
    M[1, 2] = 10.0
    M[2, 3] = 5.0
    # Expected: 10(2)(3) + 5(3)(4) = 60 + 60 = 120
    @test isapprox(TwoTimescaleResilience.pairwise_sum(M, x), 120.0; atol=1e-10)

    # Test 2.3 -- background burden with interaction
    interaction = zeros(3, 3)
    interaction[1, 2] = 1.0
    # Expected: B = 0.4 + 1.0(0.2)(0.4) = 0.4 + 0.08 = 0.48
    @test isapprox(background_index(stressors; interaction=interaction), 0.48; atol=1e-10)

    # Test 2.4 -- non-negativity
    stressors_neg = [
        BackgroundStressor("a", 1.0, -1.0),
        BackgroundStressor("b", 1.0, -2.0)
    ]
    @test background_index(stressors_neg) == 0.0
end
