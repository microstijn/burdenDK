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
            k_M = species_data["k_M"]   # somatic maintenance rate constant [p_M]/[E_G]
            E_G = species_data["E_G"]   # volume-specific cost of structure
            g_ratio = species_data["g"] # energy investment ratio [E_G]/(kap [E_m])

            # Check validity
            if !(p_Am isa Float64) || !(p_M isa Float64) || !(kap isa Float64) || !(v isa Float64) ||
               !(k_M isa Float64) || !(E_G isa Float64) || !(g_ratio isa Float64) ||
               isnan(p_Am) || isnan(p_M) || isnan(kap) || isnan(v) ||
               isnan(k_M) || isnan(E_G) || isnan(g_ratio) ||
               p_Am <= 0 || p_M <= 0 || kap <= 0 || v <= 0 ||
               k_M <= 0 || E_G <= 0 || g_ratio <= 0
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

            # Fast recovery rate (1/time): reserve mobilization, conductance over max length.
            lambda_max = v / L_m

            # Slow recovery floor (1/time): the DEB somatic maintenance rate constant
            #   k_M = [p_M] / [E_G].
            # This replaces the previous lambda_min = p_M / A_0 = [p_M]/[E_m], which normalized
            # maintenance by reserve density instead of the cost of structure and therefore forced
            # lambda_max / lambda_min == 1/kappa (the "kappa-collapse"). With k_M as the floor,
            #   lambda_max / lambda_min == (v/L_m)/k_M == g  (the energy investment ratio),
            # a primary DEB parameter that varies across species by reserve/structure economy.
            #
            # Clamp so lambda_min <= lambda_max: reserve-rich species (g < 1, k_M > v/L_m) get no
            # timescale separation -> lambda is constant -> F ~ 1 (resilient by construction).
            lambda_min = min(k_M, lambda_max)

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
                    "p_M" => p_M,
                    "k_M" => k_M,
                    "E_G" => E_G,
                    "g" => g_ratio
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

process_species_data()

