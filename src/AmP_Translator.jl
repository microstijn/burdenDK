using MAT
using JSON

function process_species_data()
    # Load MAT file
    mat_file = "data/allStat.mat"
    if !isfile(mat_file)
        println("Error: $mat_file not found.")
        return
    end

    file = matopen(mat_file)
    if !haskey(file, "allStat")
        println("Error: 'allStat' key not found in $mat_file")
        close(file)
        return
    end

    allStat = read(file, "allStat")
    close(file)

    species_names = keys(allStat)

    results = Dict{String, Any}()
    skipped_count = 0
    processed_count = 0

    for species in species_names
        species_data = allStat[species]

        try
            p_Am = species_data["p_Am"]
            p_M = species_data["p_M"]
            kap = species_data["kap"]
            v = species_data["v"]

            # Check validity
            if !(p_Am isa Float64) || !(p_M isa Float64) || !(kap isa Float64) || !(v isa Float64) ||
               isnan(p_Am) || isnan(p_M) || isnan(kap) || isnan(v) ||
               p_Am <= 0 || p_M <= 0 || kap <= 0 || v <= 0
                skipped_count += 1
                continue
            end

            # Math Translation
            E_m = p_Am / v
            A_0 = E_m

            alpha_M = p_M / (kap * p_Am)
            alpha_A = 1.0 / E_m
            alpha_R = 1.0 - kap
            alpha_G = kap

            L_m = kap * p_Am / p_M

            lambda_max = v * 0.1
            
            # Ensure the slow recovery floor is strictly smaller (e.g., 5% of max)
            # but apply an absolute minimum floor of 1e-5 to prevent division-by-zero.
            lambda_min = max(lambda_max * 0.05, 1e-5)
            
            # Safety catch: just in case v was incredibly close to 0, ensure max > min
            lambda_max = max(lambda_max, lambda_min * 2.0)
            
            # KA remains the same
            K_A = A_0 * 0.3

            # Tranche 3: JSON Structure formatting
            results[species] = Dict(
                "A0" => A_0,
                "alpha_axes" => [alpha_A, alpha_M, alpha_G, alpha_R],
                "lambda_bounds" => Dict(
                    "lambda_min" => lambda_min,
                    "lambda_max" => lambda_max,
                    "KA" => K_A
                ),
                "auxiliary_metrics" => Dict(
                    "L_m" => L_m,
                    "p_Am" => p_Am,
                    "p_M" => p_M
                )
            )
            processed_count += 1

        catch e
            skipped_count += 1
        end
    end

    println("Processed $processed_count species. Skipped $skipped_count.")

    if haskey(results, "Donax_trunculus")
        dt = results["Donax_trunculus"]
        println("\nDonax_trunculus stats:")
        println("  A0: ", dt["A0"])
        println("  alpha_M: ", dt["alpha_axes"][2])
        println("  L_m: ", dt["auxiliary_metrics"]["L_m"])
        println("  p_Am: ", dt["auxiliary_metrics"]["p_Am"])
        println("  p_M: ", dt["auxiliary_metrics"]["p_M"])
    end

    # Write JSON
    output_file = "data/AmP_Species_Library.json"
    open(output_file, "w") do f
        JSON.print(f, results, 2)
    end
    println("Results successfully written to $output_file")

    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    process_species_data()
end
