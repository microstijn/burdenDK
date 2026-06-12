using Test
using TwoTimescaleResilience

@testset "AmP SpeciesProfile Construction - Tranche 2" begin
    library_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
    library = load_amp_species_library(library_path)

    @testset "Default profile construction" begin
        profile = amp_species_profile(
            library,
            "Abatus cordatus";
            exposure_filter = aquatic_exposure_filter(7),
            moa_mapping = default_isimip_moa_mapping(),
            moa_deb_mapping = default_moa_to_deb_mapping()
        )

        @test profile isa SpeciesProfile
        @test profile.name == "Abatus cordatus"
        @test profile.exposure_filter isa ExposureFilter
        @test profile.moa_mapping isa ModeOfActionMapping
        @test profile.moa_deb_mapping isa MoAToDEBMapping
        @test profile.deb_params isa DEBAxisParams
        @test profile.buffer_params === nothing
        @test occursin("AmP", profile.description)

        # Check AmP parameters
        @test profile.deb_params.A0 ≈ 1539.9871108901957
        @test profile.deb_params.lambda_min ≈ 0.005783592166791565   # = k_M
        @test profile.deb_params.lambda_max ≈ 0.011568702452292915
    end

    @testset "Custom profile attributes" begin
        profile_custom = amp_species_profile(
            library,
            "Abatus_cordatus";
            exposure_filter = aquatic_exposure_filter(7),
            moa_mapping = default_isimip_moa_mapping(),
            moa_deb_mapping = default_moa_to_deb_mapping(),
            name = "Custom Abatus",
            description = "custom description"
        )

        @test profile_custom.name == "Custom Abatus"
        @test profile_custom.description == "custom description"
    end

    @testset "Missing species handling" begin
        @test_throws KeyError amp_species_profile(
            library,
            "Not a real species";
            exposure_filter = aquatic_exposure_filter(7),
            moa_mapping = default_isimip_moa_mapping(),
            moa_deb_mapping = default_moa_to_deb_mapping()
        )
    end
end
