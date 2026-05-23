using JSON

"""
    amp_species_key(name::AbstractString)::String

Normalize a user-facing species name to the JSON key convention.
Trims leading/trailing whitespace and replaces any run of whitespace with a single underscore.
"""
function amp_species_key(name::AbstractString)::String
    return replace(strip(name), r"\s+" => "_")
end

"""
    load_amp_species_library(path::AbstractString = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json"))

Load the AmP JSON library.
"""
function load_amp_species_library(path::AbstractString = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json"))
    return JSON.parsefile(path)
end

"""
    validate_amp_record(record)::Bool

Validate that one AmP species record contains the fields required to construct DEBAxisParams.
Throws ArgumentError for invalid records, returns true for valid ones.
"""
function validate_amp_record(record)::Bool
    if !haskey(record, "A0")
        throw(ArgumentError("Missing required field 'A0'"))
    end
    A0 = Float64(record["A0"])
    if !isfinite(A0) || A0 <= 0
        throw(ArgumentError("'A0' must be finite and > 0"))
    end

    if !haskey(record, "alpha_axes")
        throw(ArgumentError("Missing required field 'alpha_axes'"))
    end
    alpha_axes = record["alpha_axes"]
    if length(alpha_axes) != 4
        throw(ArgumentError("'alpha_axes' must have exactly length 4"))
    end
    for a in alpha_axes
        val = Float64(a)
        if !isfinite(val) || val < 0
            throw(ArgumentError("'alpha_axes' values must be finite and >= 0"))
        end
    end

    if !haskey(record, "lambda_bounds")
        throw(ArgumentError("Missing required field 'lambda_bounds'"))
    end
    lambda_bounds = record["lambda_bounds"]

    if !haskey(lambda_bounds, "KA") || !haskey(lambda_bounds, "lambda_min") || !haskey(lambda_bounds, "lambda_max")
        throw(ArgumentError("'lambda_bounds' must contain 'KA', 'lambda_min', and 'lambda_max'"))
    end

    KA = Float64(lambda_bounds["KA"])
    if !isfinite(KA) || KA <= 0
        throw(ArgumentError("'KA' must be finite and > 0"))
    end

    lambda_min = Float64(lambda_bounds["lambda_min"])
    if !isfinite(lambda_min) || lambda_min <= 0
        throw(ArgumentError("'lambda_min' must be finite and > 0"))
    end

    lambda_max = Float64(lambda_bounds["lambda_max"])
    if !isfinite(lambda_max) || lambda_max < lambda_min
        throw(ArgumentError("'lambda_max' must be finite and >= lambda_min"))
    end

    return true
end

"""
    amp_record_to_deb_params(record)::DEBAxisParams

Convert one validated AmP JSON record into DEBAxisParams.
"""
function amp_record_to_deb_params(record)::DEBAxisParams
    validate_amp_record(record)

    return DEBAxisParams(
        A0 = Float64(record["A0"]),
        alpha_axes = NTuple{4, Float64}((
            Float64(record["alpha_axes"][1]),
            Float64(record["alpha_axes"][2]),
            Float64(record["alpha_axes"][3]),
            Float64(record["alpha_axes"][4])
        )),
        lambda_min = Float64(record["lambda_bounds"]["lambda_min"]),
        lambda_max = Float64(record["lambda_bounds"]["lambda_max"]),
        KA = Float64(record["lambda_bounds"]["KA"])
    )
end

"""
    amp_species_deb_params(library, species_name::AbstractString)::DEBAxisParams

Look up a species in the parsed AmP library and return DEBAxisParams.
"""
function amp_species_deb_params(library, species_name::AbstractString)::DEBAxisParams
    if haskey(library, species_name)
        return amp_record_to_deb_params(library[species_name])
    end

    normalized_key = amp_species_key(species_name)
    if haskey(library, normalized_key)
        return amp_record_to_deb_params(library[normalized_key])
    end

    throw(KeyError(normalized_key))
end

"""
    amp_species_profile(library, species_name::AbstractString; kwargs...)

Construct a full SpeciesProfile from an AmP species record and user-supplied runtime mapping objects.
"""
function amp_species_profile(
    library,
    species_name::AbstractString;
    exposure_filter::ExposureFilter,
    moa_mapping::ModeOfActionMapping,
    moa_deb_mapping::MoAToDEBMapping,
    name::String = replace(amp_species_key(species_name), "_" => " "),
    description::String = "AmP-derived species profile"
)::SpeciesProfile
    deb_params = amp_species_deb_params(library, species_name)

    return SpeciesProfile(
        name = name,
        exposure_filter = exposure_filter,
        moa_mapping = moa_mapping,
        moa_deb_mapping = moa_deb_mapping,
        deb_params = deb_params,
        buffer_params = nothing,
        description = description
    )
end
