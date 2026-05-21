using Test
using TwoTimescaleResilience

@testset "SpeciesProfiles" begin
    profiles = [
        default_species_profile(),
        aquatic_invertebrate_profile(),
        fish_profile(),
        bivalve_profile(),
        amphibian_profile(),
        bird_profile(),
        small_mammal_profile(),
        human_profile()
    ]

    # 1. Each profile constructs successfully (handled by above array initialization).
    # 2. Each profile has exposure filter length 7.
    for p in profiles
        @test length(p.exposure_filter.multipliers) == 7
    end

    # 3. Each profile can run scalar isimip_deb_pipeline on a 7-element vector.
    vals = [1.0, 0.5, 0.1, 0.2, 0.8, 0.4, 0.3]
    for p in profiles
        res = isimip_deb_pipeline(vals, p.exposure_filter, p.moa_mapping, p.moa_deb_mapping, p.deb_params; buffer_params=p.buffer_params)
        @test res.amplification >= 1.0
    end

    # 4. Each profile can run grid pipeline on seven 2x2 layers.
    layers = [fill(0.1, 2, 2) for _ in 1:7]
    for p in profiles
        res_grid = isimip_deb_pipeline_grid(layers, p.exposure_filter, p.moa_mapping, p.moa_deb_mapping, p.deb_params; buffer_params=p.buffer_params)
        @test size(res_grid.amplification) == (2, 2)
        @test res_grid.amplification[1, 1] >= 1.0
    end

    # 5. Human exposure filter attenuates at least one environmental variable relative to aquatic exposure.
    hum_p = human_profile()
    aq_p = aquatic_invertebrate_profile()

    val = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    hum_eff = apply_exposure_filter(val, hum_p.exposure_filter)
    aq_eff = apply_exposure_filter(val, aq_p.exposure_filter)

    @test any(hum_eff .< aq_eff)
end
