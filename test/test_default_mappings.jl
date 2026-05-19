using TwoTimescaleResilience
using Test

@testset "Default DEB Mapping" begin
    @testset "Test 3.1 -- default mapping no interaction" begin
        pathogen = 0.2
        organic = 0.4
        values = [pathogen, organic]
        
        mapping = default_pathogen_organic_deb_mapping(interaction_strength=0.0)
        axes = deb_axes(values, mapping)
        
        @test isapprox(axes.assimilation, 0.2; atol=1e-8)
        @test isapprox(axes.maintenance, 0.23; atol=1e-8)
        @test isapprox(axes.growth, 0.1; atol=1e-8)
        @test isapprox(axes.reproduction, 0.11; atol=1e-8)
    end
    
    @testset "Test 3.2 -- maintenance interaction" begin
        pathogen = 0.2
        organic = 0.4
        values = [pathogen, organic]
        
        mapping = default_pathogen_organic_deb_mapping(interaction_strength=0.25)
        axes = deb_axes(values, mapping)
        
        @test isapprox(axes.assimilation, 0.2; atol=1e-8)
        @test isapprox(axes.maintenance, 0.25; atol=1e-8) # 0.23 + 0.25 * 0.2 * 0.4 = 0.25
        @test isapprox(axes.growth, 0.1; atol=1e-8)
        @test isapprox(axes.reproduction, 0.11; atol=1e-8)
    end
    
    @testset "Test 3.3 -- grid default mapping" begin
        pathogen_grid = fill(0.2, 3, 3)
        organic_grid = fill(0.4, 3, 3)
        layers = [pathogen_grid, organic_grid]
        
        mapping = default_pathogen_organic_deb_mapping(interaction_strength=0.25)
        axes_grid = deb_axes_grid(layers, mapping)
        
        @test all(axes_grid.assimilation .≈ 0.2)
        @test all(axes_grid.maintenance .≈ 0.25)
        @test all(axes_grid.growth .≈ 0.1)
        @test all(axes_grid.reproduction .≈ 0.11)
    end
end
