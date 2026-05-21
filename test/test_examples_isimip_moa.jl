using Test
using TwoTimescaleResilience

@testset "Examples Verification - ISIMIP MoA" begin
    # Smoke test the examples
    # These should run without error and define their output variables in Main or local scope

    # 1. 3x3 Demo
    try
        include(joinpath(@__DIR__, "..", "examples", "isimip_moa_deb_3x3_demo.jl"))
        @test true
    catch e
        @error "Example failed" exception=(e, catch_backtrace())
        @test false
    end

    # 2. Species Comparison
    try
        include(joinpath(@__DIR__, "..", "examples", "species_comparison_3x3_demo.jl"))
        @test true
    catch e
        @error "Example failed" exception=(e, catch_backtrace())
        @test false
    end
end
