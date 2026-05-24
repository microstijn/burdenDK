module ECOTOXParser

using CSV, DataFrames, Statistics, JSON

export parse_ecotox_data,
       summarize_ecotox_endpoints,
       write_ecotox_library_json,
       build_ecotox_toxicity_library

function parse_ecotox_data(results_path::String, tests_path::String, species_path::String, target_cas::String)
    # Load files
    results = CSV.read(results_path, DataFrame, delim='|', stringtype=String, silencewarnings=true, strict=false)
    tests = CSV.read(tests_path, DataFrame, delim='|', stringtype=String, silencewarnings=true, strict=false)
    species = CSV.read(species_path, DataFrame, delim='|', stringtype=String, silencewarnings=true, strict=false)

    # Inner join on test_id
    df_tests_results = innerjoin(results, tests, on = :test_id, makeunique=true)

    # Inner join on species_number
    df = innerjoin(df_tests_results, species, on = :species_number, makeunique=true)

    # Filter for test_cas matching target_cas (stripping whitespace)
    # Handling missing or Nothing gracefully
    target_cas_stripped = strip(target_cas)
    filter_func(x) = ismissing(x) ? false : strip(string(x)) == target_cas_stripped

    return filter(:test_cas => filter_func, df)
end

function _safe_float(x)::Union{Nothing, Float64}
    if ismissing(x) || x === nothing
        return nothing
    end
    if isa(x, Number)
        return Float64(x)
    end

    val = strip(string(x))
    if val == ""
        return nothing
    end

    parsed = tryparse(Float64, val)
    return parsed
end

function _first_existing_column(df::DataFrame, candidates::Vector{Symbol})::Union{Nothing, Symbol}
    df_names = propertynames(df)
    for c in candidates
        if c in df_names
            return c
        end
    end
    return nothing
end

end
