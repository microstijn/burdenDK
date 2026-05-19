using Test
using TwoTimescaleResilience

@testset "Tranche 1: Background Scalar Functions" begin
    # Test 1.1 -- positive part
    @test positive_part(-1.0) == 0.0
    @test positive_part(0.0) == 0.0
    @test positive_part(2.5) == 2.5

    params = BackgroundParams()

    # Test 1.2 -- zero background gives baseline margin
    @test phi_background(0.0, params) == 0.0
    @test adaptive_margin(0.0, params) == params.A0

    # Test 1.3 -- restoring force boundedness
    Bvals = range(0.0, 3.0, length=100)
    for B in Bvals
        l_B = restoring_force(B, params)
        @test l_B >= params.lambda_min
        @test l_B <= params.lambda_max
    end

    # Test 1.4 -- amplification factor baseline
    @test amplification_factor(0.0, params) == 1.0

    # Test 1.5 -- monotonicity
    # For increasing B, adaptive margin should be non-increasing
    # Restoring force should be non-increasing
    # Amplification factor should be non-decreasing
    for i in 1:length(Bvals)-1
        B1 = Bvals[i]
        B2 = Bvals[i+1]

        # We know B1 < B2
        @test adaptive_margin(B1, params) >= adaptive_margin(B2, params) - 1e-12
        @test restoring_force(B1, params) >= restoring_force(B2, params) - 1e-12
        @test amplification_factor(B1, params) <= amplification_factor(B2, params) + 1e-12
    end
end
