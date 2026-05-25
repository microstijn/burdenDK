using TwoTimescaleResilience
using Printf

function main()
    println("Running Monthly Memory Demo...")

    amp = load_amp_species_library()
    params = amp_species_deb_params(amp, "Abatus cordatus")
    ecotox = load_ecotox_library()
    memory = load_compound_memory_library()

    cadmium_cas = "7440-43-9"

    # Try finding MOR, then GRO, then REP, then fallback
    records = ecotox_filter_records(ecotox; cas=cadmium_cas, effect_code="MOR")
    if isempty(records)
        records = ecotox_filter_records(ecotox; cas=cadmium_cas, effect_code="GRO")
    end
    if isempty(records)
        records = ecotox_filter_records(ecotox; cas=cadmium_cas, effect_code="REP")
    end
    if isempty(records)
        records = ecotox_filter_records(ecotox; cas=cadmium_cas)
    end

    if isempty(records)
        error("No ECOTOX records found for Cadmium (7440-43-9). Cannot proceed with demo.")
    end

    record = first(records)

    println("Selected ECOTOX Record:")
    println("  CAS: ", record["cas_hyphenated"])
    println("  Effect: ", record["effect_code"])
    println("  NOEC: ", record["NOEC_median"])
    println("  EC50: ", record["EC50_median"])
    println()

    EC50 = record["EC50_median"]
    low = 0.0
    pulse = 3.0 * EC50

    C_t = [
        low, low,
        pulse, pulse, pulse,
        low, low, low,
        pulse, pulse,
        low, low
    ]

    state = EcotoxExposureState()

    # Output formatting
    @printf("%-6s | %-8s | %-12s | %-15s | %-14s | %-12s | %-15s\n",
            "month", "C", "B_stateful", "F_stateless", "F_stateful", "A_stateful", "lambda_stateful")
    println("-"^97)

    for (m, C) in enumerate(C_t)
        # Stateless scenario
        stateless_concs = Dict(cadmium_cas => C)
        stateless_burden = ecotox_records_to_deb_burden(stateless_concs, [record])
        stateless_response = ecotox_burden_to_response(stateless_burden, params)
        F_stateless = stateless_response.amplification

        # Stateful scenario
        stateful_concs = Dict(cadmium_cas => C)
        stateful_burden = ecotox_records_to_deb_burden_stateful!(
            state, stateful_concs, [record]; memory_library=memory
        )
        stateful_response = ecotox_burden_to_response(stateful_burden, params)

        B_stateful = get_internal_burden(state, cadmium_cas)
        F_stateful = stateful_response.amplification
        A_stateful = stateful_response.A
        lambda_stateful = stateful_response.lambda

        @printf("%-6d | %-8.3f | %-12.3f | %-15.3f | %-14.3f | %-12.3f | %-15.3f\n",
                m, C, B_stateful, F_stateless, F_stateful, A_stateful, lambda_stateful)
    end
end

main()
