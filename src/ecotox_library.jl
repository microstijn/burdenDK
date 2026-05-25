import JSON
import CSV

function load_ecotox_library(path::AbstractString = normpath(joinpath(@__DIR__, "..", "data", "ECOTOX_Toxicity_Library.json")))
    if !isfile(path)
        throw(ArgumentError("ECOTOX library file not found at path: $path"))
    end
    return JSON.parsefile(path)
end

function _record_has_key(record, key)
    return haskey(record, key)
end

function _as_float(x, fieldname)
    if ismissing(x) || isnothing(x)
        throw(ArgumentError("Field $fieldname cannot be missing or nothing"))
    end
    val = try
        Float64(x)
    catch
        throw(ArgumentError("Field $fieldname must be numeric"))
    end
    if !isfinite(val)
        throw(ArgumentError("Field $fieldname must be finite"))
    end
    return val
end

function _as_nonnegative_int(x, fieldname)
    if ismissing(x) || isnothing(x)
        throw(ArgumentError("Field $fieldname cannot be missing or nothing"))
    end
    val = try
        Int(x)
    catch
        throw(ArgumentError("Field $fieldname must be integer-like"))
    end
    if val < 0
        throw(ArgumentError("Field $fieldname must be >= 0"))
    end
    return val
end

function validate_ecotox_record(record)::Bool
    required_keys = ["cas", "cas_norm", "cas_hyphenated", "taxon_class", "effect_code", "NOEC_median", "EC50_median", "n_NOEC", "n_EC50"]
    for key in required_keys
        if !_record_has_key(record, key)
            throw(ArgumentError("Record is missing required field: $key"))
        end
    end

    cas_norm = record["cas_norm"]
    if !isa(cas_norm, AbstractString) || isempty(strip(cas_norm))
        throw(ArgumentError("cas_norm must be a non-empty string"))
    end

    cas_hyphenated = record["cas_hyphenated"]
    if !isa(cas_hyphenated, AbstractString) || isempty(strip(cas_hyphenated))
        throw(ArgumentError("cas_hyphenated must be a non-empty string"))
    end

    taxon_class = record["taxon_class"]
    if !isa(taxon_class, AbstractString) || isempty(strip(taxon_class))
        throw(ArgumentError("taxon_class must be a non-empty string"))
    end

    effect_code = record["effect_code"]
    if !isa(effect_code, AbstractString) || isempty(strip(effect_code))
        throw(ArgumentError("effect_code must be a non-empty string"))
    end

    NOEC_median = _as_float(record["NOEC_median"], "NOEC_median")
    if NOEC_median < 0
        throw(ArgumentError("NOEC_median must be >= 0"))
    end

    EC50_median = _as_float(record["EC50_median"], "EC50_median")

    if EC50_median <= NOEC_median
        throw(ArgumentError("EC50_median must be > NOEC_median"))
    end

    _as_nonnegative_int(record["n_NOEC"], "n_NOEC")
    _as_nonnegative_int(record["n_EC50"], "n_EC50")

    return true
end

function compound_retention(cas; memory_library=nothing)::Float64
    cas_norm = _normalize_cas_runtime(cas)
    if isempty(cas_norm)
        return 0.0
    end

    lib = memory_library === nothing ? load_compound_memory_library() : memory_library

    for record in lib
        if string(record["cas_norm"]) == cas_norm
            validate_compound_memory_record(record)
            return Float64(record["retention_rho_monthly"])
        end
    end

    return 0.0
end

ecotox_default_retention(cas; memory_library=nothing) = compound_retention(cas; memory_library=memory_library)

function compound_bioaccumulation_factor(cas; memory_library=nothing)::Float64
    cas_norm = _normalize_cas_runtime(cas)
    if isempty(cas_norm)
        return 1.0
    end

    lib = memory_library === nothing ? load_compound_memory_library() : memory_library

    for record in lib
        if string(record["cas_norm"]) == cas_norm
            validate_compound_memory_record(record)
            return Float64(record["bioaccumulation_factor"])
        end
    end

    return 1.0
end


Base.@kwdef mutable struct EcotoxExposureState
    internal_burdens::Dict{String, Float64} = Dict{String, Float64}()
end

function get_internal_burden(state::EcotoxExposureState, cas)::Float64
    cas_norm = _normalize_cas_runtime(cas)
    if isempty(cas_norm)
        return 0.0
    end
    return get(state.internal_burdens, cas_norm, 0.0)
end

function set_internal_burden!(state::EcotoxExposureState, cas, value)::EcotoxExposureState
    cas_norm = _normalize_cas_runtime(cas)
    if isempty(cas_norm)
        return state
    end

    val = Float64(value)
    if !isfinite(val) || val < 0.0
        throw(ArgumentError("Internal burden value must be finite and >= 0"))
    end

    state.internal_burdens[cas_norm] = val
    return state
end

function reset_internal_burdens!(state::EcotoxExposureState)::EcotoxExposureState
    empty!(state.internal_burdens)
    return state
end

function update_internal_burden!(
    state::EcotoxExposureState,
    cas,
    concentration;
    retention = nothing,
    bioaccumulation_factor = nothing,
    memory_library = nothing
)::Float64
    cas_norm = _normalize_cas_runtime(cas)
    if isempty(cas_norm)
        return 0.0
    end

    C = Float64(concentration)
    if !isfinite(C) || C < 0.0
        throw(ArgumentError("concentration must be finite and >= 0"))
    end

    rho = retention === nothing ? compound_retention(cas_norm; memory_library=memory_library) : Float64(retention)

    if !isfinite(rho) || rho < 0.0 || rho >= 1.0
        throw(ArgumentError("retention must be finite and satisfy 0.0 <= rho < 1.0"))
    end

    K = bioaccumulation_factor === nothing ? compound_bioaccumulation_factor(cas_norm; memory_library=memory_library) : Float64(bioaccumulation_factor)

    if !isfinite(K) || K < 0.0
        throw(ArgumentError("bioaccumulation_factor must be finite and >= 0.0"))
    end

    B_old = get_internal_burden(state, cas_norm)
    B_new = rho * B_old + (1.0 - rho) * K * C

    set_internal_burden!(state, cas_norm, B_new)
    return B_new
end

function ecotox_active_stress(concentration, NOEC, EC50)::Float64
    conc_val = _as_float(concentration, "concentration")
    noec_val = _as_float(NOEC, "NOEC")
    ec50_val = _as_float(EC50, "EC50")

    if noec_val < 0
        throw(ArgumentError("NOEC must be >= 0"))
    end

    if ec50_val <= noec_val
        throw(ArgumentError("EC50 must be > NOEC"))
    end

    return max(0.0, conc_val - noec_val) / (ec50_val - noec_val)
end

function ecotox_effect_to_deb_axis(effect_code)::Symbol
    if !isa(effect_code, AbstractString) || isempty(strip(effect_code))
        throw(ArgumentError("effect_code must be a non-empty string"))
    end

    code = uppercase(strip(effect_code))

    if code in ["MOR", "ITX", "MPH", "BEH", "POP", "ENZ", "FDB"]
        return :maintenance
    elseif code in ["GRO", "DVP"]
        return :growth
    elseif code in ["REP", "FEC"]
        return :reproduction
    elseif code in ["BCM", "FEED", "FED", "ING", "FOOD"]
        return :assimilation
    else
        throw(ArgumentError("Unknown ECOTOX effect code: $code"))
    end
end

function deb_axis_index(axis::Symbol)::Int
    if axis === :assimilation
        return 1
    elseif axis === :maintenance
        return 2
    elseif axis === :growth
        return 3
    elseif axis === :reproduction
        return 4
    else
        throw(ArgumentError("Unknown DEB axis: $axis"))
    end
end

function ecotox_record_to_deb_burden(concentration, record)
    validate_ecotox_record(record)

    x = ecotox_active_stress(
        concentration,
        record["NOEC_median"],
        record["EC50_median"]
    )

    axis = ecotox_effect_to_deb_axis(record["effect_code"])

    return (
        assimilation = axis === :assimilation ? x : 0.0,
        maintenance = axis === :maintenance ? x : 0.0,
        growth = axis === :growth ? x : 0.0,
        reproduction = axis === :reproduction ? x : 0.0
    )
end

function _normalize_cas_runtime(cas)::String
    if ismissing(cas) || isnothing(cas)
        return ""
    end
    s = string(cas)
    return replace(s, r"[^\d]" => "")
end

function ecotox_records_to_deb_burden(concentrations, records)
    norm_concentrations = Dict{String, Float64}()
    for (k, v) in concentrations
        k_norm = _normalize_cas_runtime(k)
        if !isempty(k_norm)
            norm_concentrations[k_norm] = Float64(v)
        end
    end

    total_A = 0.0
    total_M = 0.0
    total_G = 0.0
    total_R = 0.0

    for record in records
        validate_ecotox_record(record)

        cas_norm = record["cas_norm"]
        if !haskey(norm_concentrations, cas_norm)
            continue
        end

        C = norm_concentrations[cas_norm]
        burden = ecotox_record_to_deb_burden(C, record)

        total_A += burden.assimilation
        total_M += burden.maintenance
        total_G += burden.growth
        total_R += burden.reproduction
    end

    return (
        assimilation = total_A,
        maintenance = total_M,
        growth = total_G,
        reproduction = total_R
    )
end

function ecotox_records_to_deb_burden_stateful!(
    state::EcotoxExposureState,
    concentrations,
    records;
    retention = nothing,
    bioaccumulation_factor = nothing,
    memory_library = nothing
)
    norm_concentrations = Dict{String, Float64}()
    for (k, v) in concentrations
        k_norm = _normalize_cas_runtime(k)
        if !isempty(k_norm)
            norm_concentrations[k_norm] = Float64(v)
        end
    end

    total_A = 0.0
    total_M = 0.0
    total_G = 0.0
    total_R = 0.0

    for record in records
        validate_ecotox_record(record)

        cas_norm = record["cas_norm"]
        if !haskey(norm_concentrations, cas_norm)
            continue
        end

        C = norm_concentrations[cas_norm]

        rho = nothing
        if retention isa Number
            rho = Float64(retention)
        elseif retention isa Dict
            norm_retention = Dict(string(_normalize_cas_runtime(k)) => v for (k,v) in retention)
            rho = get(norm_retention, cas_norm, nothing)
        end

        K = nothing
        if bioaccumulation_factor isa Number
            K = Float64(bioaccumulation_factor)
        elseif bioaccumulation_factor isa Dict
            norm_k = Dict(string(_normalize_cas_runtime(k)) => v for (k,v) in bioaccumulation_factor)
            K = get(norm_k, cas_norm, nothing)
        end

        B = update_internal_burden!(state, cas_norm, C; retention=rho, bioaccumulation_factor=K, memory_library=memory_library)
        burden = ecotox_record_to_deb_burden(B, record)

        total_A += burden.assimilation
        total_M += burden.maintenance
        total_G += burden.growth
        total_R += burden.reproduction
    end

    return (
        assimilation = total_A,
        maintenance = total_M,
        growth = total_G,
        reproduction = total_R
    )
end

function ecotox_burden_to_response(burden, params::DEBAxisParams)
    A = deb_adaptive_margin(burden, params)
    lambda = restoring_force_from_margin(A, params)
    F = amplification_from_margin(A, params)

    return (
        A = A,
        lambda = lambda,
        amplification = F
    )
end

function ecotox_filter_records(
    library;
    cas = nothing,
    taxon_class = nothing,
    effect_code = nothing
)
    result = Vector{Dict{String, Any}}()

    # Pre-process filters
    cas_list = cas === nothing ? nothing : (cas isa Vector ? _normalize_cas_runtime.(cas) : [_normalize_cas_runtime(cas)])
    if cas_list !== nothing
        filter!(x -> !isempty(x), cas_list)
    end

    tax_list = taxon_class === nothing ? nothing : (taxon_class isa Vector ? lowercase.(strip.(taxon_class)) : [lowercase(strip(taxon_class))])
    eff_list = effect_code === nothing ? nothing : (effect_code isa Vector ? uppercase.(strip.(effect_code)) : [uppercase(strip(effect_code))])

    for record in library
        match = true

        if cas_list !== nothing && !isempty(cas_list)
            rec_cas = get(record, "cas_norm", "")
            if !(rec_cas in cas_list)
                match = false
            end
        end

        if match && tax_list !== nothing && !isempty(tax_list)
            rec_tax = lowercase(strip(string(get(record, "taxon_class", ""))))
            if !(rec_tax in tax_list)
                match = false
            end
        end

        if match && eff_list !== nothing && !isempty(eff_list)
            rec_eff = uppercase(strip(string(get(record, "effect_code", ""))))
            if !(rec_eff in eff_list)
                match = false
            end
        end

        if match
            push!(result, record)
        end
    end

    return result
end

function ecotox_records_for_taxon(library, taxon_class; cas=nothing, effect_code=nothing)
    return ecotox_filter_records(library; taxon_class=taxon_class, cas=cas, effect_code=effect_code)
end

function load_compound_memory_library(
    path::AbstractString = normpath(joinpath(@__DIR__, "..", "data", "Compound_Memory_Library.csv"))
)
    if !isfile(path)
        throw(ArgumentError("Compound memory library file not found at path: $path"))
    end
    df = CSV.File(path)
    records = []
    for row in df
        record = Dict{String, Any}()
        for k in propertynames(row)
            record[string(k)] = row[k]
        end
        push!(records, record)
    end
    return records
end

function validate_compound_memory_record(record)::Bool
    required_keys = ["cas_norm", "cas_hyphenated", "chemical_name", "memory_class", "retention_rho_monthly", "bioaccumulation_factor", "basis", "confidence", "notes"]
    for key in required_keys
        if !haskey(record, key)
            throw(ArgumentError("Record is missing required field: $key"))
        end
    end

    cas_norm = string(record["cas_norm"])
    if isempty(strip(cas_norm)) || !occursin(r"^\d+$", strip(cas_norm))
        throw(ArgumentError("cas_norm must be a non-empty digits-only string"))
    end

    for key in ["cas_hyphenated", "chemical_name", "memory_class", "basis", "confidence"]
        val = record[key]
        if ismissing(val) || isnothing(val) || (isa(val, AbstractString) && isempty(strip(string(val))))
            throw(ArgumentError("$key must be a non-empty string"))
        end
    end

    rho_raw = record["retention_rho_monthly"]
    if ismissing(rho_raw) || isnothing(rho_raw)
        throw(ArgumentError("retention_rho_monthly cannot be missing"))
    end
    rho = try
        Float64(rho_raw)
    catch
        throw(ArgumentError("retention_rho_monthly must be numeric"))
    end

    if !isfinite(rho) || rho < 0.0 || rho >= 1.0
        throw(ArgumentError("retention_rho_monthly must be finite and satisfy 0.0 <= rho < 1.0"))
    end

    k_raw = record["bioaccumulation_factor"]
    if ismissing(k_raw) || isnothing(k_raw)
        throw(ArgumentError("bioaccumulation_factor cannot be missing"))
    end
    k_val = try
        Float64(k_raw)
    catch
        throw(ArgumentError("bioaccumulation_factor must be numeric"))
    end

    if !isfinite(k_val) || k_val < 0.0
        throw(ArgumentError("bioaccumulation_factor must be finite and >= 0.0"))
    end

    return true
end
