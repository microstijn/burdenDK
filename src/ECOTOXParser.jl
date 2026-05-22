module ECOTOXParser

using CSV, DataFrames

export parse_ecotox_data

function parse_ecotox_data(results_path::String, tests_path::String, target_cas::String)
    # Load files
    results = CSV.read(results_path, DataFrame, delim='|', stringtype=String, silencewarnings=true, strict=false)
    tests = CSV.read(tests_path, DataFrame, delim='|', stringtype=String, silencewarnings=true, strict=false)

    # Inner join on test_id
    df = innerjoin(results, tests, on = :test_id, makeunique=true)

    # Filter for test_cas matching target_cas (stripping whitespace)
    # Handling missing or Nothing gracefully
    target_cas_stripped = strip(target_cas)
    filter_func(x) = ismissing(x) ? false : strip(string(x)) == target_cas_stripped

    return filter(:test_cas => filter_func, df)
end

end
