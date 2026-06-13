# pmoa_stressor_routing.jl — load the pMoA stressor-routing table and turn a stressor
# class into per-axis weights for the margin pipeline. Mirrors load_compound_memory_library.
#
# The table (data/pMoA_Stressor_Routing.csv) routes AGGREGATE water-quality stressors
# (microplastics, BOD, salinity, …) onto the four DEB axes by physiological mode of action,
# replacing the prototype's inline w_A/w_M/w_G/w_R tuning knobs with a declared, documented
# mechanistic assignment. Basis: docs/notes/pmoa_stressor_routing_basis.md.
#
# `axis_split` is "1.0" (primary only), "p/s" (e.g. "0.7/0.3"), or "NA".
# `primary_axis == NONE` flags a rate modifier (temperature) — NOT a pressure axis.

const _PMOA_AXES = ("assimilation", "maintenance", "growth", "reproduction")

function load_pmoa_stressor_routing(
    path::AbstractString = normpath(joinpath(@__DIR__, "..", "data", "pMoA_Stressor_Routing.csv"))
)
    if !isfile(path)
        throw(ArgumentError("pMoA stressor-routing file not found at path: $path"))
    end
    records = Dict{String, Any}[]
    for row in CSV.File(path)
        record = Dict{String, Any}()
        for k in propertynames(row)
            record[string(k)] = row[k]
        end
        push!(records, record)
    end
    return records
end

# "0.7/0.3" -> (0.7, 0.3); "1.0" -> (1.0, 0.0); "NA"/missing/"" -> (1.0, 0.0)
function _parse_axis_split(s)::Tuple{Float64, Float64}
    (s === missing || s === nothing) && return (1.0, 0.0)
    str = strip(string(s))
    (isempty(str) || uppercase(str) == "NA") && return (1.0, 0.0)
    if occursin("/", str)
        parts = split(str, "/")
        length(parts) == 2 || throw(ArgumentError("axis_split must be 'p/s' or a single number; got '$str'"))
        return (parse(Float64, strip(parts[1])), parse(Float64, strip(parts[2])))
    end
    return (parse(Float64, str), 0.0)
end

_blank(v) = v === missing || v === nothing || (isa(v, AbstractString) && isempty(strip(string(v))))

function validate_pmoa_stressor_record(record)::Bool
    required = ["stressor_class", "primary_axis", "secondary_axis", "axis_split",
                "tissue_retention_rho", "ambient_condition", "pmoa_basis", "confidence", "mechanism_notes"]
    for k in required
        haskey(record, k) || throw(ArgumentError("pMoA record missing required field: $k"))
    end
    prim = lowercase(strip(string(record["primary_axis"])))
    if prim != "none" && !(prim in _PMOA_AXES)
        throw(ArgumentError("primary_axis must be one of $(_PMOA_AXES) or NONE; got '$(record["primary_axis"])'"))
    end
    if !_blank(record["secondary_axis"])
        sec = lowercase(strip(string(record["secondary_axis"])))
        sec in _PMOA_AXES || throw(ArgumentError("secondary_axis must be one of $(_PMOA_AXES) or blank; got '$(record["secondary_axis"])'"))
    end
    if prim != "none"
        wp, ws = _parse_axis_split(record["axis_split"])
        (wp < 0 || ws < 0) && throw(ArgumentError("axis_split weights must be >= 0"))
        if !_blank(record["secondary_axis"]) && !isapprox(wp + ws, 1.0; atol = 1e-6)
            throw(ArgumentError("axis_split must sum to 1 when a secondary axis is given (got $wp + $ws) for '$(record["stressor_class"])'"))
        end
        rho = record["tissue_retention_rho"]
        if !_blank(rho) && uppercase(strip(string(rho))) != "NA"
            r = parse(Float64, strip(string(rho)))
            (isfinite(r) && 0.0 <= r < 1.0) || throw(ArgumentError("tissue_retention_rho must satisfy 0 <= rho < 1; got $r"))
        end
    end
    return true
end

is_rate_modifier(record)::Bool = lowercase(strip(string(record["primary_axis"]))) == "none"

# stressor class -> (assimilation, maintenance, growth, reproduction) weights summing to 1
# (all-zero for a rate modifier such as temperature).
function stressor_axis_weights(record::AbstractDict)
    w = Dict(a => 0.0 for a in _PMOA_AXES)
    is_rate_modifier(record) && return (assimilation = 0.0, maintenance = 0.0, growth = 0.0, reproduction = 0.0)
    wp, ws = _parse_axis_split(record["axis_split"])
    w[lowercase(strip(string(record["primary_axis"])))] = wp
    if !_blank(record["secondary_axis"])
        w[lowercase(strip(string(record["secondary_axis"])))] += ws
    end
    return (assimilation = w["assimilation"], maintenance = w["maintenance"],
            growth = w["growth"], reproduction = w["reproduction"])
end

function _lookup_stressor(stressor_class::AbstractString, routing)
    table = routing === nothing ? load_pmoa_stressor_routing() : routing
    for rec in table
        string(rec["stressor_class"]) == stressor_class && return rec
    end
    throw(KeyError(stressor_class))
end

stressor_axis_weights(stressor_class::AbstractString; routing = nothing) =
    stressor_axis_weights(_lookup_stressor(stressor_class, routing))

# monthly tissue-memory rho for a stressor class (0.0 for ambient conditions; NaN for rate modifiers)
function stressor_tissue_rho(stressor_class::AbstractString; routing = nothing)::Float64
    rec = _lookup_stressor(stressor_class, routing)
    is_rate_modifier(rec) && return NaN
    rho = rec["tissue_retention_rho"]
    (_blank(rho) || uppercase(strip(string(rho))) == "NA") && return 0.0
    return parse(Float64, strip(string(rho)))
end
