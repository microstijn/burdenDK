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

    # NOTE: the recovery curve is now linear in margin (no half-saturation
    # constant), so a "KA" field is no longer required or read. Any "KA" present
    # in a legacy JSON record is ignored. See src/deb_axes.jl.
    if !haskey(lambda_bounds, "lambda_min") || !haskey(lambda_bounds, "lambda_max")
        throw(ArgumentError("'lambda_bounds' must contain 'lambda_min' and 'lambda_max'"))
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
        lambda_max = Float64(record["lambda_bounds"]["lambda_max"])
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

# --- Stage-resolved capacity (life-stage integration) --------------------------------------
# Record-facing wrappers over the pure DEB math in deb_axes.jl. These are ADDITIVE: the
# whole-organism path (amp_record_to_deb_params) is unchanged, and a species whose JSON record
# carries no "ontogeny" block simply has no stage API (a clear error is thrown).

"""
    amp_species_record(library, species_name)

Return the raw AmP JSON record for a species (with name normalization), or throw `KeyError`.
"""
function amp_species_record(library, species_name::AbstractString)
    if haskey(library, species_name)
        return library[species_name]
    end
    normalized_key = amp_species_key(species_name)
    if haskey(library, normalized_key)
        return library[normalized_key]
    end
    throw(KeyError(normalized_key))
end

"""
    has_ontogeny(record)::Bool

True if `record` carries an ontogeny block with the fields needed for stage-resolved capacity.
"""
function has_ontogeny(record)::Bool
    haskey(record, "ontogeny") || return false
    ont = record["ontogeny"]
    for k in ("s_M", "L_b", "L_j", "L_p", "L_i")
        haskey(ont, k) || return false
    end
    return true
end

"""
    validate_ontogeny_record(record)::Bool

Validate the ontogeny block required for stage-resolved capacity. Throws `ArgumentError` if the
block is missing or inconsistent; returns `true` otherwise.
"""
function validate_ontogeny_record(record)::Bool
    if !haskey(record, "ontogeny")
        throw(ArgumentError("Record has no 'ontogeny' block; stage-resolved capacity is unavailable for this species"))
    end
    ont = record["ontogeny"]
    for k in ("s_M", "L_b", "L_j", "L_p", "L_i")
        if !haskey(ont, k)
            throw(ArgumentError("'ontogeny' is missing required field '$k'"))
        end
    end
    s_M = Float64(ont["s_M"]); L_b = Float64(ont["L_b"]); L_j = Float64(ont["L_j"])
    L_p = Float64(ont["L_p"]); L_i = Float64(ont["L_i"])
    for (nm, val) in (("s_M", s_M), ("L_b", L_b), ("L_j", L_j), ("L_p", L_p), ("L_i", L_i))
        if !isfinite(val) || val <= 0
            throw(ArgumentError("'ontogeny.$nm' must be finite and > 0"))
        end
    end
    if s_M < 1.0
        throw(ArgumentError("'ontogeny.s_M' must be >= 1"))
    end
    if !(L_b <= L_p <= L_i) || !(L_b <= L_j <= L_i)
        throw(ArgumentError("'ontogeny' lengths must satisfy 0 < L_b <= L_p <= L_i and L_b <= L_j <= L_i"))
    end
    return true
end

# Extract the DEB quantities the stage API needs from a record. The conductance is recovered as
# v = p_Am / A0, because the offline translator sets A0 = E_m = p_Am / v.
function _amp_stage_inputs(record)
    validate_amp_record(record)
    validate_ontogeny_record(record)
    A0 = Float64(record["A0"])
    aux = record["auxiliary_metrics"]
    p_Am = Float64(aux["p_Am"])
    k_M = Float64(aux["k_M"])
    v = p_Am / A0
    ont = record["ontogeny"]
    return (
        base = amp_record_to_deb_params(record),
        v = v, k_M = k_M,
        L_b = Float64(ont["L_b"]), L_j = Float64(ont["L_j"]),
        L_p = Float64(ont["L_p"]), L_i = Float64(ont["L_i"]), s_M = Float64(ont["s_M"])
    )
end

"""
    deb_params_at_length(record, L) -> DEBStageProfile
    deb_params_at_length(library, species_name, L) -> DEBStageProfile

Stage-resolved capacity at structural length `L`, with `lambda_max = v_eff(L)/L`. The
whole-organism `A0`, `alpha_axes`, and `lambda_min = k_M` are preserved; only the recovery
ceiling becomes length-dependent (younger/smaller -> faster recovery; abj acceleration via s_M).
"""
function deb_params_at_length(record::AbstractDict, L::Real)::DEBStageProfile
    f = _amp_stage_inputs(record)
    return deb_stage_profile(f.base, L; v=f.v, lambda_min=f.k_M,
                             L_b=f.L_b, L_j=f.L_j, L_p=f.L_p, L_i=f.L_i, s_M=f.s_M)
end

function deb_params_at_length(library, species_name::AbstractString, L::Real)::DEBStageProfile
    return deb_params_at_length(amp_species_record(library, species_name), L)
end

"""
    deb_params_for_stage(record, stage) -> DEBStageProfile
    deb_params_for_stage(library, species_name, stage) -> DEBStageProfile

Discrete-stage convenience wrapper over [`deb_params_at_length`](@ref) at a representative
structural length: `:juvenile -> (L_b + L_p)/2` and `:adult -> L_i`. `:adult` reproduces the
whole-organism profile (since `L_i = s_M * L_m`, the recovery ceiling `v_eff(L_i)/L_i = v/L_m`).
The embryo stage (`L < L_b`) has no environmental feeding, so it is out of scope for the discrete
API; query `deb_params_at_length(record, L)` with `L < L_b` directly if genuinely needed.
"""
function deb_params_for_stage(record::AbstractDict, stage::Symbol)::DEBStageProfile
    f = _amp_stage_inputs(record)
    L = if stage === :juvenile
        (f.L_b + f.L_p) / 2
    elseif stage === :adult
        f.L_i
    elseif stage === :embryo
        throw(ArgumentError("Stage :embryo has no environmental feeding (assimilation axis not applicable) and is out of scope for the discrete stage API; use deb_params_at_length(record, L) with L < L_b if an embryo profile is genuinely needed"))
    else
        throw(ArgumentError("Unknown stage $(repr(stage)); expected :juvenile or :adult"))
    end
    return deb_params_at_length(record, L)
end

function deb_params_for_stage(library, species_name::AbstractString, stage::Symbol)::DEBStageProfile
    return deb_params_for_stage(amp_species_record(library, species_name), stage)
end
