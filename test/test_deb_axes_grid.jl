using TwoTimescaleResilience
using Test

@testset "DEB Axes Grid Mapping" begin
    @testset "Test 2.1 -- grid mapping matches scalar mapping" begin
        L1 = fill(0.5, 2, 2)
        L2 = fill(0.2, 2, 2)
        layers = [L1, L2]
        W = [
            1.0 0.0;
            0.0 2.0;
            0.5 0.5;
            0.0 1.0
        ]
        mapping = DEBAxisMapping(W=W)
        axes_grid = deb_axes_grid(layers, mapping)
        
        @test all(axes_grid.assimilation .≈ 0.5)
        @test all(axes_grid.maintenance .≈ 0.4)
        @test all(axes_grid.growth .≈ 0.35)
        @test all(axes_grid.reproduction .≈ 0.2)
    end
    
    @testset "Test 2.2 -- shape preservation" begin
        L1 = rand(7, 11)
        L2 = rand(7, 11)
        layers = [L1, L2]
        W = zeros(4, 2)
        mapping = DEBAxisMapping(W=W)
        axes_grid = deb_axes_grid(layers, mapping)
        
        @test size(axes_grid.assimilation) == (7, 11)
        @test size(axes_grid.maintenance) == (7, 11)
        @test size(axes_grid.growth) == (7, 11)
        @test size(axes_grid.reproduction) == (7, 11)
    end
    
    @testset "Test 2.3 -- NaN propagation" begin
        L1 = fill(0.5, 3, 3)
        L2 = fill(0.2, 3, 3)
        L1[2, 2] = NaN
        layers = [L1, L2]
        W = [
            1.0 0.0;
            0.0 2.0;
            0.5 0.5;
            0.0 1.0
        ]
        mapping = DEBAxisMapping(W=W)
        axes_grid = deb_axes_grid(layers, mapping)
        
        @test isnan(axes_grid.assimilation[2, 2])
        @test isnan(axes_grid.maintenance[2, 2])
        @test isnan(axes_grid.growth[2, 2])
        @test isnan(axes_grid.reproduction[2, 2])
        
        @test !isnan(axes_grid.assimilation[1, 1])
    end
    
    @testset "Test 2.4 -- DEB adaptive margin grid" begin
        L1 = fill(0.1, 2, 2)
        L2 = fill(0.2, 2, 2)
        layers = [L1, L2]
        W = [
            1.0 0.0;
            0.0 1.0;
            1.5 0.0;
            0.0 2.0
        ]
        # sA=0.1, sM=0.2, sG=0.15, sR=0.4
        mapping = DEBAxisMapping(W=W)
        axes_grid = deb_axes_grid(layers, mapping)
        
        params = DEBAxisParams(A0=10.0, alpha_axes=(1.0, 2.0, 3.0, 4.0))
        A_grid = deb_adaptive_margin_grid(axes_grid, params)
        
        scalar_A = deb_adaptive_margin([0.1, 0.2, 0.15, 0.4], params)
        
        @test all(A_grid .≈ scalar_A)
    end
    
    @testset "Test 2.5 -- amplification lower bound" begin
        L1 = fill(0.5, 2, 2)
        L2 = fill(0.5, 2, 2)
        layers = [L1, L2]
        W = [
            1.0 0.0;
            0.0 1.0;
            0.0 0.0;
            0.0 0.0
        ]
        mapping = DEBAxisMapping(W=W)
        axes_grid = deb_axes_grid(layers, mapping)
        
        params = DEBAxisParams(A0=1.0)
        A_grid = deb_adaptive_margin_grid(axes_grid, params)
        
        # A <= A0 everywhere because alpha_axes > 0 and sA > 0, sM > 0
        @test all(A_grid .<= params.A0)
        
        F_grid = amplification_from_margin_grid(A_grid, params)
        
        @test all(F_grid .>= 1.0)
    end
end
