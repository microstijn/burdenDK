using Statistics
using DataFrames

include("../src/ECOTOXParser.jl")
using .ECOTOXParser

function main()
    target_cas = "7664417"
    results_path = "test_data/sample_results.txt"
    tests_path = "test_data/sample_tests.txt"

    df = parse_ecotox_data(results_path, tests_path, target_cas)

    if isempty(df)
        println("Sample too small for Ammonia targets, but pipeline executed successfully.")
        return
    end

    # Clean Data
    # 1. Filter for endpoint NOEC or EC50 (ignoring missing and stripping whitespace)
    function is_target_endpoint(x)
        if ismissing(x)
            return false
        end
        val = strip(string(x))
        return val == "NOEC" || val == "EC50"
    end

    df_filtered = filter(:endpoint => is_target_endpoint, df)

    if isempty(df_filtered)
        println("Sample too small for Ammonia targets, but pipeline executed successfully.")
        return
    end

    # 2. Filter out missing, "", "NR", "NC" from conc1_mean
    function is_valid_conc(x)
        if ismissing(x) || x === nothing
            return false
        end
        val = strip(string(x))
        return val != "" && val != "NR" && val != "NC"
    end

    df_valid_conc = filter(:conc1_mean => is_valid_conc, df_filtered)

    if isempty(df_valid_conc)
        println("Sample too small for Ammonia targets, but pipeline executed successfully.")
        return
    end

    # 3. Parse conc1_mean to Float64
    df_valid_conc.conc1_mean_num = parse.(Float64, strip.(string.(df_valid_conc.conc1_mean)))

    # 4. Calculate and print Medians
    df_noec = filter(:endpoint => x -> strip(string(x)) == "NOEC", df_valid_conc)
    df_ec50 = filter(:endpoint => x -> strip(string(x)) == "EC50", df_valid_conc)

    if !isempty(df_noec)
        median_noec = median(df_noec.conc1_mean_num)
        println("Median NOEC: $median_noec")
    end

    if !isempty(df_ec50)
        median_ec50 = median(df_ec50.conc1_mean_num)
        println("Median EC50: $median_ec50")
    end
end

main()
