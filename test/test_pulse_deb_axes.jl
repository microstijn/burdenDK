using TwoTimescaleResilience
using Test

@testset "DEB-Like Pulse Mapping" begin
    @testset "Test 5.1 -- pulse axes from burden matrix" begin
        # ammonia, pesticide
        D = [
            0.5 0.0;
            0.0 0.8;
            0.5 0.8;
            0.0 0.0
        ]
        
        W = [
            0.0 0.5;
            0.6 0.0;
            0.0 0.3;
            0.0 0.1
        ]
        mapping = DEBAxisMapping(W=W)
        pulse_axes = pulse_deb_axes_timeseries(D, mapping)
        
        @test isapprox(pulse_axes.assimilation[1], 0.0; atol=1e-8)
        @test isapprox(pulse_axes.maintenance[1], 0.6 * 0.5; atol=1e-8)
        
        @test isapprox(pulse_axes.assimilation[2], 0.5 * 0.8; atol=1e-8)
        @test isapprox(pulse_axes.maintenance[2], 0.0; atol=1e-8)
        
        @test isapprox(pulse_axes.assimilation[3], 0.5 * 0.8; atol=1e-8)
        @test isapprox(pulse_axes.maintenance[3], 0.6 * 0.5; atol=1e-8)
        
        @test isapprox(pulse_axes.assimilation[4], 0.0; atol=1e-8)
        @test isapprox(pulse_axes.maintenance[4], 0.0; atol=1e-8)
    end
    
    @testset "Test 5.2 -- pulse interaction" begin
        D = [
            0.5 0.8;
            0.0 0.8
        ]
        W = [
            0.0 0.5;
            0.6 0.0;
            0.0 0.3;
            0.0 0.1
        ]
        Gamma_A = zeros(2, 2)
        Gamma_M = zeros(2, 2)
        Gamma_G = zeros(2, 2)
        Gamma_R = zeros(2, 2)
        
        eta = 0.2
        Gamma_M[1, 2] = eta
        
        mapping = DEBAxisMapping(W=W, interactions=[Gamma_A, Gamma_M, Gamma_G, Gamma_R])
        pulse_axes = pulse_deb_axes_timeseries(D, mapping)
        
        @test isapprox(pulse_axes.assimilation[1], 0.5 * 0.8; atol=1e-8)
        @test isapprox(pulse_axes.maintenance[1], 0.6 * 0.5 + eta * 0.5 * 0.8; atol=1e-8)
        @test isapprox(pulse_axes.growth[1], 0.3 * 0.8; atol=1e-8)
        @test isapprox(pulse_axes.reproduction[1], 0.1 * 0.8; atol=1e-8)
        
        # for second time step, ammonia is 0, so interaction is 0
        @test isapprox(pulse_axes.maintenance[2], 0.0; atol=1e-8)
    end
    
    @testset "Test 5.3 -- total margin decreases under pulse stress" begin
        D = [0.5 0.8; 0.2 0.1]
        W = [0.0 0.5; 0.6 0.0; 0.0 0.3; 0.0 0.1]
        mapping = DEBAxisMapping(W=W)
        pulse_axes = pulse_deb_axes_timeseries(D, mapping)
        
        params = DEBAxisParams()
        A_bg = 0.8
        A_total = total_deb_margin_timeseries(A_bg, pulse_axes, params)
        
        @test all(A_total .<= A_bg)
    end
    
    @testset "Test 5.4 -- zero pulse axes leave margin unchanged" begin
        D = zeros(5, 2)
        W = [0.0 0.5; 0.6 0.0; 0.0 0.3; 0.0 0.1]
        mapping = DEBAxisMapping(W=W)
        pulse_axes = pulse_deb_axes_timeseries(D, mapping)
        
        params = DEBAxisParams()
        A_bg = 0.8
        A_total = total_deb_margin_timeseries(A_bg, pulse_axes, params)
        
        @test all(A_total .== A_bg)
    end
end
