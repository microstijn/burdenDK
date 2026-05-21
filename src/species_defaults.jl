export SpeciesProfile
export aquatic_invertebrate_profile, fish_profile, bivalve_profile, amphibian_profile
export bird_profile, small_mammal_profile, human_profile, default_species_profile

Base.@kwdef struct SpeciesProfile
    name::String
    exposure_filter::ExposureFilter
    moa_mapping::ModeOfActionMapping
    moa_deb_mapping::MoAToDEBMapping
    deb_params::DEBAxisParams
    buffer_params::Union{Nothing, ConditionBufferParams} = nothing
    description::String = ""
end

function default_species_profile()
    return SpeciesProfile(
        name = "default",
        exposure_filter = default_exposure_filter(7),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping(),
        deb_params = DEBAxisParams(),
        description = "Default uncalibrated profile."
    )
end

function aquatic_invertebrate_profile()
    return SpeciesProfile(
        name = "aquatic_invertebrate",
        exposure_filter = aquatic_exposure_filter(7),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping(),
        # higher maintenance sensitivity hypothesis
        deb_params = DEBAxisParams(alpha_axes = (0.25, 0.45, 0.15, 0.15)),
        description = "Hypothetical aquatic invertebrate sensitivity test profile."
    )
end

function fish_profile()
    return SpeciesProfile(
        name = "fish",
        exposure_filter = aquatic_exposure_filter(7),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping(),
        # higher assimilation sensitivity hypothesis
        deb_params = DEBAxisParams(alpha_axes = (0.40, 0.30, 0.15, 0.15)),
        description = "Hypothetical fish sensitivity test profile."
    )
end

function bivalve_profile()
    return SpeciesProfile(
        name = "bivalve",
        exposure_filter = aquatic_exposure_filter(7),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping(),
        # higher growth/feeding sensitivity hypothesis
        deb_params = DEBAxisParams(alpha_axes = (0.25, 0.30, 0.30, 0.15)),
        description = "Hypothetical bivalve sensitivity test profile."
    )
end

function amphibian_profile()
    return SpeciesProfile(
        name = "amphibian",
        exposure_filter = aquatic_exposure_filter(7),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping(),
        # high maintenance and reproduction sensitivity
        deb_params = DEBAxisParams(alpha_axes = (0.20, 0.40, 0.10, 0.30)),
        description = "Hypothetical amphibian sensitivity test profile."
    )
end

function bird_profile()
    return SpeciesProfile(
        name = "bird",
        exposure_filter = ExposureFilter(multipliers=[0.05, 0.1, 0.1, 0.2, 0.1, 0.5, 0.5], name="bird_exposure"),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping(),
        deb_params = DEBAxisParams(alpha_axes = (0.25, 0.30, 0.15, 0.30)),
        description = "Hypothetical bird sensitivity test profile."
    )
end

function small_mammal_profile()
    return SpeciesProfile(
        name = "small_mammal",
        exposure_filter = ExposureFilter(multipliers=[0.05, 0.1, 0.1, 0.2, 0.1, 0.5, 0.2], name="small_mammal_exposure"),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping(),
        deb_params = DEBAxisParams(alpha_axes = (0.30, 0.40, 0.15, 0.15)),
        description = "Hypothetical small mammal sensitivity test profile."
    )
end

function human_profile()
    return SpeciesProfile(
        name = "human",
        exposure_filter = human_exposure_filter(),
        moa_mapping = default_isimip_moa_mapping(),
        moa_deb_mapping = default_moa_to_deb_mapping(),
        # Emphasize maintenance/immune/repair/detoxification
        deb_params = DEBAxisParams(alpha_axes = (0.20, 0.50, 0.15, 0.15)),
        description = "Hypothetical human profile. Output is a vulnerability multiplier, not disease risk."
    )
end
