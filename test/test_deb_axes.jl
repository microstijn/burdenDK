using TwoTimescaleResilience
using Test

@testset "DEB Axes Scalar Mapping" begin
    @testset "Test 1.1 -- matrix mapping without interactions" begin
        values = [0.5, 0.2]
        W = [
            1.0 0.0;
            0.0 2.0;
            0.5 0.5;
            0.0 1.0
        ]
        mapping = DEBAxisMapping(W=W)
        axes = deb_axes(values, mapping)

        @test isapprox(axes.assimilation, 0.5; atol=1e-8)
        @test isapprox(axes.maintenance, 0.4; atol=1e-8)
        @test isapprox(axes.growth, 0.35; atol=1e-8)
        @test isapprox(axes.reproduction, 0.2; atol=1e-8)
    end

    @testset "Test 1.2 -- axis-specific interactions" begin
        values = [0.5, 0.2]
        W = zeros(4, 2)
        Gamma_A = zeros(2, 2)
        Gamma_M = zeros(2, 2)
        Gamma_G = zeros(2, 2)
        Gamma_R = zeros(2, 2)
        Gamma_M[1, 2] = 3.0

        interactions = [Gamma_A, Gamma_M, Gamma_G, Gamma_R]
        mapping = DEBAxisMapping(W=W, interactions=interactions)

        axes = deb_axes(values, mapping)

        @test isapprox(axes.assimilation, 0.0; atol=1e-8)
        @test isapprox(axes.maintenance, 0.3; atol=1e-8) # 3.0 * 0.5 * 0.2 = 0.3
        @test isapprox(axes.growth, 0.0; atol=1e-8)
        @test isapprox(axes.reproduction, 0.0; atol=1e-8)
    end

    @testset "Test 1.3 -- non-negativity clamp" begin
        values = [0.5, 0.2]
        W = [
            -1.0 0.0;
            0.0 -2.0;
            0.5 0.5;
            0.0 1.0
        ]
        mapping_clamp = DEBAxisMapping(W=W, clamp_nonnegative=true)
        axes_clamp = deb_axes(values, mapping_clamp)
        @test isapprox(axes_clamp.assimilation, 0.0; atol=1e-8)
        @test isapprox(axes_clamp.maintenance, 0.0; atol=1e-8)
        @test isapprox(axes_clamp.growth, 0.35; atol=1e-8)

        mapping_noclamp = DEBAxisMapping(W=W, clamp_nonnegative=false)
        axes_noclamp = deb_axes(values, mapping_noclamp)
        @test isapprox(axes_noclamp.assimilation, -0.5; atol=1e-8)
        @test isapprox(axes_noclamp.maintenance, -0.4; atol=1e-8)
    end

    @testset "Test 1.4 -- optional unit clamp" begin
        values = [0.8, 0.8]
        W = [
            2.0 0.0;
            0.0 2.0;
            0.5 0.5;
            0.0 1.0
        ]
        mapping_clamp = DEBAxisMapping(W=W, clamp_unit=true)
        axes = deb_axes(values, mapping_clamp)

        @test isapprox(axes.assimilation, 1.0; atol=1e-8)
        @test isapprox(axes.maintenance, 1.0; atol=1e-8)
        @test isapprox(axes.growth, 0.8; atol=1e-8)
    end

    @testset "Test 1.5 -- adaptive margin" begin
        s = [0.1, 0.2, 0.3, 0.4]
        params = DEBAxisParams(A0=10.0, alpha_axes=(1.0, 2.0, 3.0, 4.0))

        A = deb_adaptive_margin(s, params)
        # Expected: 10 - (1*0.1 + 2*0.2 + 3*0.3 + 4*0.4) = 10 - (0.1 + 0.4 + 0.9 + 1.6) = 10 - 3.0 = 7.0
        @test isapprox(A, 7.0; atol=1e-8)

        # Test NamedTuple interface
        axes_nt = (assimilation=0.1, maintenance=0.2, growth=0.3, reproduction=0.4)
        A_nt = deb_adaptive_margin(axes_nt, params)
        @test isapprox(A_nt, 7.0; atol=1e-8)
    end

    @testset "Test 1.6 -- restoring force boundedness" begin
        params = DEBAxisParams(lambda_min=0.04, lambda_max=1.0, KA=0.30)
        for A in range(-2.0, 2.0, length=50)
            lambda_val = restoring_force_from_margin(A, params)
            @test lambda_val >= params.lambda_min
            @test lambda_val <= params.lambda_max
        end
    end

    @testset "Test 1.7 -- amplification baseline" begin
        params = DEBAxisParams()
        amp = amplification_from_margin(params.A0, params)
        @test isapprox(amp, 1.0; atol=1e-8)
    end
end
