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
        params = DEBAxisParams(lambda_min=0.04, lambda_max=1.0)
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

@testset "Response from Precomputed Axis Impairments: Tranche 2" begin
    params = TwoTimescaleResilience.DEBAxisParams(
        A0 = 1000.0,
        alpha_axes = (1.0, 1.0, 1.0, 1.0),
        lambda_max = 0.05
    )

    # 1. Zero impairments -> Q_t=0, A_t=A0, F_t=1
    zeros_imp = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
    res_zero = compute_adaptive_margin_response_from_impairment(zeros_imp, params)
    @test res_zero.Q_t == 0.0
    @test res_zero.A_t == params.A0
    @test res_zero.F_t == 1.0
    @test isnan(res_zero.X_assimilation) # default X_axis=nothing

    # 2. Maintenance = 0.5
    half_m_imp = (assimilation=0.0, maintenance=0.5, growth=0.0, reproduction=0.0)
    res_half = compute_adaptive_margin_response_from_impairment(half_m_imp, params)
    @test res_half.Q_t > 0.0
    @test res_half.A_t < params.A0
    @test res_half.F_t >= 1.0

    # 3. Check Q_t = sum(w * E)
    @test res_half.Q_t ≈ (res_half.w_maintenance * 0.5)

    # 4. A_t = A0 * max(...)
    @test res_half.A_t ≈ params.A0 * max(1e-6, 1.0 - res_half.Q_t)

    # 5. lambda_t & 6. F_t
    lam_t = TwoTimescaleResilience.restoring_force_from_margin(res_half.A_t, params)
    @test res_half.lambda_t ≈ lam_t
    @test res_half.F_t ≈ (TwoTimescaleResilience.restoring_force_from_margin(params.A0, params) / lam_t)

    # 7. Invalid E < 0
    @test_throws ArgumentError compute_adaptive_margin_response_from_impairment((assimilation=-0.1, maintenance=0.0, growth=0.0, reproduction=0.0), params)

    # 8. Invalid E > 1
    @test_throws ArgumentError compute_adaptive_margin_response_from_impairment((assimilation=1.1, maintenance=0.0, growth=0.0, reproduction=0.0), params)

    # 9. NaN E
    @test_throws ArgumentError compute_adaptive_margin_response_from_impairment((assimilation=NaN, maintenance=0.0, growth=0.0, reproduction=0.0), params)

    # 10. Inf E
    @test_throws ArgumentError compute_adaptive_margin_response_from_impairment((assimilation=Inf, maintenance=0.0, growth=0.0, reproduction=0.0), params)

    # 11. Invalid A_floor_fraction
    @test_throws ArgumentError compute_adaptive_margin_response_from_impairment(zeros_imp, params, A_floor_fraction=-0.1)
    @test_throws ArgumentError compute_adaptive_margin_response_from_impairment(zeros_imp, params, A_floor_fraction=1.5)

    # 12. Matches EC50 response path exactly
    X_vals = (assimilation=1.0, maintenance=2.0, growth=0.5, reproduction=0.0)
    # Get standard response
    res_std = compute_adaptive_margin_response(X_vals, params; response_mode="ec50_anchored_fractional_impairment")
    
    # Calculate E manually 
    E_vals = (
        assimilation = X_vals.assimilation / (1 + X_vals.assimilation),
        maintenance = X_vals.maintenance / (1 + X_vals.maintenance),
        growth = X_vals.growth / (1 + X_vals.growth),
        reproduction = X_vals.reproduction / (1 + X_vals.reproduction)
    )
    
    res_pre = compute_adaptive_margin_response_from_impairment(E_vals, params; X_axis=X_vals)
    @test res_pre.Q_t ≈ res_std.Q_t
    @test res_pre.A_t ≈ res_std.A_t
    @test res_pre.F_t ≈ res_std.F_t
    @test res_pre.X_maintenance == res_std.X_maintenance

end
