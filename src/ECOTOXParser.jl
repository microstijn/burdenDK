module ECOTOXParser

using CSV, DataFrames, Statistics, JSON

export parse_ecotox_data,
       summarize_ecotox_endpoints,
       write_ecotox_library_json,
       build_ecotox_toxicity_library,
       build_ecotox_toxicity_library_multi,
       normalize_cas,
       hyphenate_cas

function normalize_cas(cas)::String
    if ismissing(cas) || cas === nothing
        return ""
    end
    # strip whitespace and keep only digits
    s = strip(string(cas))
    return filter(isdigit, s)
end

function hyphenate_cas(cas)::String
    norm = normalize_cas(cas)
    if length(norm) < 4
        return norm
    end
    prefix = norm[1:end-3]
    middle = norm[end-2:end-1]
    check = norm[end:end]
    return "$(prefix)-$(middle)-$(check)"
end

function parse_ecotox_data(results_path::String, tests_path::String, species_path::String, target_cas::String)
    target_cas_norm = normalize_cas(target_cas)
    if target_cas_norm == ""
        throw(ArgumentError("target_cas normalizes to an empty string"))
    end

    # Load files
    results = CSV.read(results_path, DataFrame, delim='|', stringtype=String, silencewarnings=true, strict=false)
    tests = CSV.read(tests_path, DataFrame, delim='|', stringtype=String, silencewarnings=true, strict=false)
    species = CSV.read(species_path, DataFrame, delim='|', stringtype=String, silencewarnings=true, strict=false)

    # Inner join on test_id
    df_tests_results = innerjoin(results, tests, on = :test_id, makeunique=true)

    # Inner join on species_number
    df = innerjoin(df_tests_results, species, on = :species_number, makeunique=true)

    # Filter for test_cas matching normalized target_cas
    # Handling missing or Nothing gracefully
    filter_func(x) = normalize_cas(x) == target_cas_norm

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

function summarize_ecotox_endpoints(df::DataFrame; cas::AbstractString)
    cas_norm = normalize_cas(cas)
    if cas_norm == ""
        throw(ArgumentError("cas normalizes to an empty string"))
    end
    cas_hyphenated = hyphenate_cas(cas_norm)

    taxon_col = _first_existing_column(df, [:class, :taxonomic_class, :species_class, :taxon_class])
    if taxon_col === nothing
        throw(ArgumentError("No taxonomic class column could be detected."))
    end

    # We will assume effect code is in :effect and endpoint is in :endpoint
    if !(:effect in propertynames(df)) || !(:endpoint in propertynames(df))
        throw(ArgumentError("Missing required :effect or :endpoint column."))
    end

    conc_col = _first_existing_column(df, [:conc1_mean_num, :conc1_mean, :concentration, :result_value, :effect_conc, :endpoint_value])
    if conc_col === nothing
        throw(ArgumentError("No concentration column could be detected."))
    end

    # Keep only required columns and add normalized versions
    df_clean = DataFrame()
    df_clean.taxon_class = [ismissing(x) ? missing : strip(string(x)) for x in df[!, taxon_col]]
    df_clean.effect_code = [ismissing(x) ? missing : uppercase(strip(string(x))) for x in df[!, :effect]]
    df_clean.endpoint_type = [ismissing(x) ? missing : uppercase(strip(string(x))) for x in df[!, :endpoint]]
    df_clean.conc = [_safe_float(x) for x in df[!, conc_col]]

    # Drop rows with missing taxonomic class, effect code, or endpoint type, or nothing in conc
    df_filtered = filter(row -> !ismissing(row.taxon_class) && !ismissing(row.effect_code) && !ismissing(row.endpoint_type) && row.conc !== nothing, df_clean)

    df_filtered = filter(row -> row.endpoint_type == "NOEC" || row.endpoint_type == "EC50", df_filtered)

    if nrow(df_filtered) == 0
        return DataFrame(
            cas = String[],
            cas_norm = String[],
            cas_hyphenated = String[],
            taxon_class = String[],
            effect_code = String[],
            NOEC_median = Union{Missing, Float64}[],
            EC50_median = Union{Missing, Float64}[],
            n_NOEC = Int[],
            n_EC50 = Int[]
        )
    end

    # Group by taxon_class and effect_code
    summary_df = combine(groupby(df_filtered, [:taxon_class, :effect_code])) do subdf
        noecs = subdf[subdf.endpoint_type .== "NOEC", :conc]
        ec50s = subdf[subdf.endpoint_type .== "EC50", :conc]

        n_NOEC = length(noecs)
        n_EC50 = length(ec50s)

        NOEC_median = n_NOEC > 0 ? median(noecs) : missing
        EC50_median = n_EC50 > 0 ? median(ec50s) : missing

        (NOEC_median = NOEC_median, EC50_median = EC50_median, n_NOEC = n_NOEC, n_EC50 = n_EC50)
    end

    if nrow(summary_df) > 0
        summary_df.cas = fill(cas_norm, nrow(summary_df))
        summary_df.cas_norm = fill(cas_norm, nrow(summary_df))
        summary_df.cas_hyphenated = fill(cas_hyphenated, nrow(summary_df))
    else
        summary_df.cas = String[]
        summary_df.cas_norm = String[]
        summary_df.cas_hyphenated = String[]
    end

    # Reorder columns to match schema
    select!(summary_df, :cas, :cas_norm, :cas_hyphenated, :taxon_class, :effect_code, :NOEC_median, :EC50_median, :n_NOEC, :n_EC50)

    return summary_df
end

function write_ecotox_library_json(summary_df::DataFrame, output_path::AbstractString)
    # Convert DataFrame to array of dicts to correctly handle missing => null via JSON.jl
    records = []
    for row in eachrow(summary_df)
        d = Dict{String, Any}()
        for name in propertynames(summary_df)
            val = row[name]
            d[string(name)] = ismissing(val) ? nothing : val
        end
        push!(records, d)
    end

    # Create parent directories if they don't exist
    dir = dirname(output_path)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end

    open(output_path, "w") do io
        JSON.print(io, records, 2)
    end
    return output_path
end

function build_ecotox_toxicity_library(
    results_path::String,
    tests_path::String,
    species_path::String,
    target_cas::String;
    output_path::Union{Nothing, AbstractString}=nothing
)
    df = parse_ecotox_data(results_path, tests_path, species_path, target_cas)
    summary_df = summarize_ecotox_endpoints(df; cas=target_cas)

    if output_path !== nothing
        write_ecotox_library_json(summary_df, output_path)
    end

    return summary_df
end

function build_ecotox_toxicity_library_multi(
    results_path::String,
    tests_path::String,
    species_path::String,
    target_cas_list::AbstractVector{<:AbstractString};
    output_path::Union{Nothing, AbstractString}=nothing,
    skip_empty::Bool=true
)
    if isempty(target_cas_list)
        throw(ArgumentError("target_cas_list cannot be empty"))
    end

    # Deduplicate CAS based on normalized values, preserving first-seen order
    seen = Set{String}()
    deduped_cas_list = String[]
    for cas in target_cas_list
        norm = normalize_cas(cas)
        if norm == ""
            throw(ArgumentError("Invalid or empty CAS provided: '$cas'"))
        end
        if !(norm in seen)
            push!(seen, norm)
            push!(deduped_cas_list, cas)
        end
    end

    summaries = DataFrame[]

    for cas in deduped_cas_list
        df = parse_ecotox_data(results_path, tests_path, species_path, cas)
        if nrow(df) == 0
            if skip_empty
                @warn "No valid ECOTOX records found for CAS '$cas'. Skipping."
                continue
            else
                @warn "No valid ECOTOX records found for CAS '$cas'. Including empty summary."
            end
        end

        summary_df = summarize_ecotox_endpoints(df; cas=cas)
        push!(summaries, summary_df)
    end

    if isempty(summaries)
        # Create an empty DataFrame with the correct schema
        combined_df = DataFrame(
            cas = String[],
            cas_norm = String[],
            cas_hyphenated = String[],
            taxon_class = String[],
            effect_code = String[],
            NOEC_median = Union{Missing, Float64}[],
            EC50_median = Union{Missing, Float64}[],
            n_NOEC = Int[],
            n_EC50 = Int[]
        )
    else
        combined_df = vcat(summaries...; cols=:union)
    end

    if output_path !== nothing
        write_ecotox_library_json(combined_df, output_path)
    end

    return combined_df
end

end
