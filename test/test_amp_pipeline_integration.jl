using Test
using TwoTimescaleResilience

@testset "AmP-derived SpeciesProfile Pipeline Integration - Tranche 3" begin
    library_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
    library = load_amp_species_library(library_path)

    profile = amp_species_profile(
        library,
        "Abatus cordatus";
        exposure_filter = aquatic_exposure_filter(7),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping()
    )

    @testset "Single-cell zero stress smoke test" begin
        values = zeros(7)

        result = isimip_deb_pipeline(
            values,
            profile.exposure_filter,
            profile.moa_mapping,
            profile.moa_deb_mapping,
            profile.deb_params
        )

        @test result isa NamedTuple
        @test haskey(result, :effective_values)
        @test haskey(result, :modes)
        @test haskey(result, :axes)
        @test haskey(result, :A)
        @test haskey(result, :lambda)
        @test haskey(result, :amplification)

        @test result.amplification ≈ 1.0 atol=1e-8
        @test result.A ≈ profile.deb_params.A0 atol=1e-8

        lambda0 = restoring_force_from_margin(profile.deb_params.A0, profile.deb_params)
        @test result.lambda ≈ lambda0 atol=1e-8
    end

    @testset "Single-cell monotonic stress smoke test" begin
        values_low = fill(0.1, 7)
        values_high = fill(0.5, 7)

        result_low = isimip_deb_pipeline(
            values_low,
            profile.exposure_filter,
            profile.moa_mapping,
            profile.moa_deb_mapping,
            profile.deb_params
        )

        result_high = isimip_deb_pipeline(
            values_high,
            profile.exposure_filter,
            profile.moa_mapping,
            profile.moa_deb_mapping,
            profile.deb_params
        )

        @test result_low.amplification >= 1.0
        @test result_high.amplification >= 1.0
        @test result_high.A <= result_low.A + eps()
        @test result_high.lambda <= result_low.lambda + eps()
        @test result_high.amplification >= result_low.amplification - eps()
    end

    @testset "Grid pipeline zero stress smoke test" begin
        layers_zero = [zeros(2, 2) for _ in 1:7]

        grid_result = isimip_deb_pipeline_grid(
            layers_zero,
            profile.exposure_filter,
            profile.moa_mapping,
            profile.moa_deb_mapping,
            profile.deb_params
        )

        @test grid_result isa NamedTuple
        @test haskey(grid_result, :effective_layers)
        @test haskey(grid_result, :modes)
        @test haskey(grid_result, :axes)
        @test haskey(grid_result, :A)
        @test haskey(grid_result, :lambda)
        @test haskey(grid_result, :amplification)

        @test size(grid_result.A) == (2, 2)
        @test size(grid_result.lambda) == (2, 2)
        @test size(grid_result.amplification) == (2, 2)

        @test all(grid_result.amplification .≈ 1.0)
        @test all(grid_result.A .≈ profile.deb_params.A0)

        lambda0 = restoring_force_from_margin(profile.deb_params.A0, profile.deb_params)
        @test all(grid_result.lambda .≈ lambda0)
    end
end
