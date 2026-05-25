using Pkg
try
    Pkg.activate(joinpath(@__DIR__, ".."))
catch err
    @warn "Could not activate project" exception=(err, catch_backtrace())
end

using TwoTimescaleResilience

function main()
    println("============================================================")
    println("TRANCHE 1: MULTI-SPECIES × MULTI-COMPOUND SCALAR DEMO")
    println("============================================================")

    # Load libraries
    amp_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
    ecotox_path = joinpath(@__DIR__, "..", "data", "ECOTOX_Toxicity_Library.json")

    amp = load_amp_species_library(amp_path)
    ecotox = load_ecotox_library(ecotox_path)

    println("This demo combines selected ECOTOX taxon-level toxicity records with AmP species-level defense parameters to demonstrate the runtime pipeline. Taxon matching will be refined in later workflows.")

    # ------------------------------------------------------------
    # Species selection
    # ------------------------------------------------------------
    preferred_species = ["Abatus cordatus", "Podarcis muralis", "Thalia democratica"]
    selected_species_names = String[]
    species_params = Dict{String, Any}()

    for sp_name in preferred_species
        if haskey(amp, sp_name)
            push!(selected_species_names, sp_name)
            species_params[sp_name] = amp_species_deb_params(amp, sp_name)
        else
            @warn "Preferred species $sp_name missing from AmP library."
        end
    end

    # Fallback to get 3 species
    if length(selected_species_names) < 3
        for (sp_name, data) in amp
            if !(sp_name in selected_species_names)
                try
                    params = amp_species_deb_params(amp, sp_name)
                    push!(selected_species_names, sp_name)
                    species_params[sp_name] = params
                    if length(selected_species_names) == 3
                        break
                    end
                catch
                    # Ignore invalid params
                end
            end
        end
    end

    println("\n--- Selected Species ---")
    for sp_name in selected_species_names
        params = species_params[sp_name]
        lambda0 = restoring_force_from_margin(params.A0, params)
        println("Species: $sp_name")
        println("  A0: $(params.A0)")
        println("  alpha_axes: $(params.alpha_axes)")
        println("  lambda_min: $(params.lambda_min)")
        println("  lambda_max: $(params.lambda_max)")
        println("  KA: $(params.KA)")
        println("  baseline lambda0: $lambda0\n")
    end

    # ------------------------------------------------------------
    # Compound selection
    # ------------------------------------------------------------
    preferred_compounds = ["7647-14-5", "7440-43-9", "7440-50-8"]
    selected_records = []

    # Get a list of all available CAS from the library just in case
    available_cas = unique([get(rec, "cas_hyphenated", get(rec, "cas_norm", get(rec, "cas", ""))) for rec in ecotox])

    # Helper to find a record
    function find_valid_record(target_cas)
        recs = TwoTimescaleResilience.ecotox_filter_records(ecotox; cas=target_cas)
        # effect code priority
        for eff in ["MOR", "GRO", "REP"]
            for rec in recs
                if get(rec, "effect_code", "") == eff && validate_ecotox_record(rec)
                    return rec
                end
            end
        end
        # any valid
        for rec in recs
            if validate_ecotox_record(rec)
                return rec
            end
        end
        return nothing
    end

    selected_cas_list = String[]
    for cas in preferred_compounds
        rec = find_valid_record(cas)
        if rec !== nothing
            push!(selected_records, rec)
            push!(selected_cas_list, cas)
        else
            @warn "Preferred compound CAS $cas has no valid records. Will select fallback."
        end
    end

    # Fallback to get 3 compounds
    if length(selected_records) < 3
        for cas in available_cas
            if !(cas in selected_cas_list) && cas != ""
                rec = find_valid_record(cas)
                if rec !== nothing
                    push!(selected_records, rec)
                    push!(selected_cas_list, cas)
                    if length(selected_records) == 3
                        break
                    end
                end
            end
        end
    end

    println("\n--- Selected Compound Records ---")
    for rec in selected_records
        cas_hyph = get(rec, "cas_hyphenated", "")
        cas_norm = get(rec, "cas_norm", "")
        taxon = get(rec, "taxon_class", "")
        eff = get(rec, "effect_code", "")
        routed = ecotox_effect_to_deb_axis(eff)
        n_m = get(rec, "NOEC_median", NaN)
        e_m = get(rec, "EC50_median", NaN)
        n_noec = get(rec, "n_NOEC", 0)
        n_ec50 = get(rec, "n_EC50", 0)
        println("CAS Hyphenated: $cas_hyph")
        println("CAS Norm:       $cas_norm")
        println("Taxon Class:    $taxon")
        println("Effect Code:    $eff")
        println("Routed DEB:     $routed")
        println("NOEC_median:    $n_m")
        println("EC50_median:    $e_m")
        println("n_NOEC:         $n_noec")
        println("n_EC50:         $n_ec50\n")
    end

    # ------------------------------------------------------------
    # Define concentration scenarios
    # ------------------------------------------------------------
    concentrations_baseline = Dict{String, Float64}()
    concentrations_moderate = Dict{String, Float64}()
    concentrations_high = Dict{String, Float64}()

    for rec in selected_records
        cas = get(rec, "cas_hyphenated", "")
        noec = Float64(get(rec, "NOEC_median", 0.0))
        ec50 = Float64(get(rec, "EC50_median", 0.0))
        concentrations_baseline[cas] = noec
        concentrations_moderate[cas] = ec50
        concentrations_high[cas] = 3.0 * ec50
    end

    scenarios = [
        ("Baseline", concentrations_baseline),
        ("Moderate", concentrations_moderate),
        ("High", concentrations_high)
    ]

    # ------------------------------------------------------------
    # Run scalar responses
    # ------------------------------------------------------------

    # Table 1: Compound stress summary per scenario
    println("\n--- Table 1: Compound stress summary per scenario ---")
    println(rpad("Scenario", 12), rpad("Compound CAS", 15), rpad("Effect", 8), rpad("Axis", 15), rpad("Conc", 12), rpad("NOEC", 12), rpad("EC50", 12), "Active Stress")

    # Keep track for sanity checks
    active_stresses = Dict()
    scenario_burdens = Dict()
    scenario_species_results = Dict()

    for (scen_name, concs) in scenarios
        active_stresses[scen_name] = Dict()
        for rec in selected_records
            cas = get(rec, "cas_hyphenated", "")
            eff = get(rec, "effect_code", "")
            axis = string(ecotox_effect_to_deb_axis(eff))
            noec = Float64(get(rec, "NOEC_median", 0.0))
            ec50 = Float64(get(rec, "EC50_median", 0.0))
            c = concs[cas]
            astress = ecotox_active_stress(c, noec, ec50)
            active_stresses[scen_name][cas] = astress
            println(rpad(scen_name, 12), rpad(cas, 15), rpad(eff, 8), rpad(axis, 15), rpad(round(c, digits=4), 12), rpad(round(noec, digits=4), 12), rpad(round(ec50, digits=4), 12), round(astress, digits=4))
        end
    end

    # Table 2: Aggregated DEB burden per scenario
    println("\n--- Table 2: Aggregated DEB burden per scenario ---")
    println(rpad("Scenario", 12), rpad("Assimilation", 15), rpad("Maintenance", 15), rpad("Growth", 15), "Reproduction")

    for (scen_name, concs) in scenarios
        burden = ecotox_records_to_deb_burden(concs, selected_records)
        scenario_burdens[scen_name] = burden
        println(rpad(scen_name, 12), rpad(round(burden.assimilation, digits=4), 15), rpad(round(burden.maintenance, digits=4), 15), rpad(round(burden.growth, digits=4), 15), round(burden.reproduction, digits=4))
    end

    # Table 3: Species response per scenario
    println("\n--- Table 3: Species response per scenario ---")
    println(rpad("Scenario", 12), rpad("Species", 20), rpad("A", 12), rpad("lambda", 12), "Amplification")

    for (scen_name, concs) in scenarios
        burden = scenario_burdens[scen_name]
        scenario_species_results[scen_name] = Dict()
        for sp_name in selected_species_names
            params = species_params[sp_name]
            response = ecotox_burden_to_response(burden, params)
            scenario_species_results[scen_name][sp_name] = response
            println(rpad(scen_name, 12), rpad(sp_name, 20), rpad(round(response.A, digits=4), 12), rpad(round(response.lambda, digits=4), 12), round(response.amplification, digits=4))
        end
    end

    println("\n--- Expected scientific behavior ---")
    println("- Baseline scenario uses concentration = NOEC, so active stress should be zero or near zero.")
    println("- Moderate scenario uses concentration = EC50, so each selected compound should have active stress near 1.")
    println("- High scenario uses concentration = 3 * EC50, so active stress should exceed 1.")
    println("- The same aggregated ECOTOX attack burden is applied to each species.")
    println("- Different species produce different A, lambda, and amplification because their AmP-derived defense parameters differ.")
    println("- Amplification should generally increase from baseline to moderate to high.")

    # ------------------------------------------------------------
    # Sanity checks
    # ------------------------------------------------------------
    println("\n--- Sanity Checks ---")
    checks_passed = true

    # 1. Baseline active stress approx zero
    for (cas, st) in active_stresses["Baseline"]
        if abs(st) > 1e-4
            @warn "Baseline active stress for $cas is $st (expected ~0)"
            checks_passed = false
        end
    end

    # 2. Moderate active stress approx 1
    for (cas, st) in active_stresses["Moderate"]
        if abs(st - 1.0) > 1e-4
            @warn "Moderate active stress for $cas is $st (expected ~1)"
            checks_passed = false
        end
    end

    # 3. High active stress >= Moderate
    for cas in keys(active_stresses["High"])
        shigh = active_stresses["High"][cas]
        smod = active_stresses["Moderate"][cas]
        if shigh < smod
            @warn "High active stress ($shigh) < Moderate active stress ($smod) for $cas"
            checks_passed = false
        end
    end

    # 4. Species checks
    for sp_name in selected_species_names
        r_base = scenario_species_results["Baseline"][sp_name]
        r_high = scenario_species_results["High"][sp_name]
        if r_high.amplification < r_base.amplification - 1e-4
            @warn "High F ($(r_high.amplification)) < Baseline F ($(r_base.amplification)) for $sp_name"
            checks_passed = false
        end
        if r_high.A > r_base.A + 1e-4
            @warn "High A ($(r_high.A)) > Baseline A ($(r_base.A)) for $sp_name"
            checks_passed = false
        end
    end

    if checks_passed
        println("Scalar multi-species multi-compound sanity checks passed.")
    else
        println("Some sanity checks failed. Check warnings.")
    end

    println("\nDone.")
end

main()
