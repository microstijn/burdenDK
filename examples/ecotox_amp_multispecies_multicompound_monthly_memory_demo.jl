using TwoTimescaleResilience
using Printf
using CSV
using DataFrames
using CairoMakie

"""
3 species × 3 compounds × 12 months monthly memory diagnostic.

This example incorporates an analytical warm-start initialisation for chemical memory (B_t)
to represent embodied prior environmental exposure before reported month 1. This is an analytical
calculation only and does not involve physiological condition memory (Z_t) or formal mixture toxicity.

Species:
- Abatus cordatus
- Podarcis muralis
- Thalia democratica

Compounds:
- Sodium chloride (7647-14-5) - low memory
- Cadmium (7440-43-9) - persistent
- Mercury (7439-97-6) - highly persistent

The diagnostic computes monthly internal burdens B_t, active stress x_t,
DEB-axis burdens, adaptive margin A_t, restoring force lambda_t, and amplification F_t.
Note: Current aggregation is additive/axis-wise; formal mixture toxicity is not implemented.
Physiological condition memory Z_t and raster integration are also not implemented.

This diagnostic is meant to show how ambient concentration C_t (pressure proxy), internal burden B_t (retained memory),
active stress x_t, DEB-axis burdens, adaptive margin A_t (physiological capacity proxy),
restoring force lambda_t, and amplification F_t evolve through time.
"""

function analytical_initial_burden(rho::Float64, K::Float64, C_bg::Float64, spinup_months::Int; B0::Float64 = 0.0)
    if !(0.0 <= rho < 1.0)
        throw(ArgumentError("retention_rho_monthly must be finite and satisfy 0.0 <= rho < 1.0. Got $rho"))
    end
    if K <= 0.0 || !isfinite(K)
        throw(ArgumentError("bioaccumulation factor K must be positive and finite. Got $K"))
    end
    if C_bg < 0.0 || !isfinite(C_bg)
        throw(ArgumentError("background concentration C_bg must be >= 0 and finite. Got $C_bg"))
    end
    if spinup_months < 0
        throw(ArgumentError("spinup_months must be >= 0. Got $spinup_months"))
    end
    if B0 < 0.0 || !isfinite(B0)
        throw(ArgumentError("initial B0 must be >= 0 and finite. Got $B0"))
    end

    if spinup_months == 0
        return B0
    end

    # Explicit calculation: B_n = rho^n * B_0 + K * C_bg * (1 - rho^n)
    return (rho^spinup_months) * B0 + K * C_bg * (1.0 - (rho^spinup_months))
end

function background_for_target_burden(target_B::Float64, rho::Float64, K::Float64, spinup_months::Int; B0::Float64 = 0.0)
    if !(0.0 <= rho < 1.0)
        throw(ArgumentError("retention_rho_monthly must be finite and satisfy 0.0 <= rho < 1.0. Got $rho"))
    end
    if K <= 0.0 || !isfinite(K)
        throw(ArgumentError("bioaccumulation factor K must be positive and finite. Got $K"))
    end
    if target_B < 0.0 || !isfinite(target_B)
        throw(ArgumentError("target_B must be >= 0 and finite. Got $target_B"))
    end
    if spinup_months < 0
        throw(ArgumentError("spinup_months must be >= 0. Got $spinup_months"))
    end
    if B0 < 0.0 || !isfinite(B0)
        throw(ArgumentError("initial B0 must be >= 0 and finite. Got $B0"))
    end

    denom = K * (1.0 - (rho^spinup_months))
    if denom <= 0.0 || !isfinite(denom)
        throw(ArgumentError("Denominator K*(1-rho^spinup_months) must be strictly positive and finite. Got $denom"))
    end

    C_bg = (target_B - (rho^spinup_months) * B0) / denom

    if C_bg < 0.0
        throw(ArgumentError("Computed C_bg is negative ($C_bg) for target $target_B with B0 $B0. This is invalid."))
    end

    return C_bg
end

function generate_scenario_concentrations(EC50::Float64)
    low = 0.0
    pulse = 10.0 * EC50
    mod_pulse = 5.0 * EC50
    return [
        low, low, low,
        pulse, pulse,
        low, low, low,
        mod_pulse, mod_pulse,
        low, low
    ]
end

function get_ecotox_record(ecotox, cas)
    function is_valid_record(r)
        has_noec = haskey(r, "NOEC_median") && r["NOEC_median"] !== nothing
        has_ec50 = haskey(r, "EC50_median") && r["EC50_median"] !== nothing
        if has_noec && has_ec50
            return Float64(r["EC50_median"]) > Float64(r["NOEC_median"])
        end
        return false
    end
    valid_filter = is_valid_record

    records = filter(valid_filter, ecotox_filter_records(ecotox; cas=cas, effect_code="MOR"))
    if isempty(records)
        records = filter(valid_filter, ecotox_filter_records(ecotox; cas=cas, effect_code="GRO"))
    end
    if isempty(records)
        records = filter(valid_filter, ecotox_filter_records(ecotox; cas=cas, effect_code="REP"))
    end
    if isempty(records)
        records = filter(valid_filter, ecotox_filter_records(ecotox; cas=cas))
    end

    if isempty(records)
        error("No ECOTOX records with finite NOEC/EC50 found for CAS \$cas. Cannot proceed.")
    end
    return first(records)
end

function main()
    println("Loading libraries...")
    amp = load_amp_species_library()
    ecotox = load_ecotox_library()
    memory = load_compound_memory_library()

    species_list = ["Abatus cordatus", "Podarcis muralis", "Thalia democratica"]
    params_map = Dict(sp => amp_species_deb_params(amp, sp) for sp in species_list)

    compounds = [
        (cas="7647-14-5", name="Sodium chloride"),
        (cas="7440-43-9", name="Cadmium"),
        (cas="7439-97-6", name="Mercury")
    ]

    records_map = Dict{String, Any}()
    for c in compounds
        rec = get_ecotox_record(ecotox, c.cas)
        # Deepcopy the JSON.Object and explicitly set type
        r2 = Dict{String, Any}()
        for (k, v) in rec
            r2[k] = v
        end
        r2["NOEC_median"] = Float64(rec["NOEC_median"])
        r2["EC50_median"] = Float64(rec["EC50_median"])
        records_map[c.cas] = r2
    end

    # Tranche 3: Compute stateful monthly burden with Scenarios
    compound_results = []
    species_results = []

    scenarios = ["zero_start", "analytical_warm_start"]
    mixture_methods = ["additive_axis_burden", "axis_toxic_unit_sum"]

    ec50_ref = records_map["7440-43-9"]["EC50_median"] # pedogogical scalar based on Cd EC50
    C_t = generate_scenario_concentrations(ec50_ref)

    for scenario in scenarios
        for mixture_method in mixture_methods
        # Initialize an independent state for each species
        species_states = Dict(sp => EcotoxExposureState() for sp in species_list)

        spinup_used = (scenario == "analytical_warm_start")
        spinup_months = spinup_used ? 24 : 0
        spinup_method = spinup_used ? "analytical_closed_form" : "none"

        # Pre-compute the initial analytical burden values mapping cas_norm -> value
        # Note: B_t applies to the internal chemical state, independent of species DEB params.
        # But we do keep one EcotoxExposureState per species just conceptually for the pipeline.
        analytical_initial_map = Dict{String, Float64}()
        analytical_C_bg_map = Dict{String, Float64}()

        for c in compounds
            record = records_map[c.cas]
            cas_norm = record["cas_norm"]

            if !spinup_used
                analytical_initial_map[cas_norm] = 0.0
                analytical_C_bg_map[cas_norm] = 0.0
                continue
            end

            # Analytical Warm Start logic
            if c.name == "Sodium chloride"
                analytical_initial_map[cas_norm] = 0.0
                analytical_C_bg_map[cas_norm] = 0.0
                continue
            end

            rho = TwoTimescaleResilience.compound_retention(cas_norm; memory_library=memory)
            K = TwoTimescaleResilience.compound_bioaccumulation_factor(cas_norm; memory_library=memory)

            has_noec = haskey(record, "NOEC_median") && record["NOEC_median"] !== nothing
            has_ec50 = haskey(record, "EC50_median") && record["EC50_median"] !== nothing

            target_B_initial = 0.0
            if has_noec && Float64(record["NOEC_median"]) > 0.0
                target_B_initial = 0.5 * Float64(record["NOEC_median"])
            elseif has_ec50 && Float64(record["EC50_median"]) > 0.0
                target_B_initial = 0.05 * Float64(record["EC50_median"])
            else
                # No usable NOEC or EC50 - skip warm-start for this compound
                analytical_initial_map[cas_norm] = 0.0
                analytical_C_bg_map[cas_norm] = 0.0
                continue
            end

            # Use inverse helper to find background
            # If compound has rho=0 or K=1 or typical values, catch ArgumentErrors or safe behavior.
            try
                C_bg = background_for_target_burden(target_B_initial, rho, K, spinup_months; B0=0.0)
                B_init = analytical_initial_burden(rho, K, C_bg, spinup_months; B0=0.0)
                analytical_initial_map[cas_norm] = B_init
                analytical_C_bg_map[cas_norm] = C_bg
            catch e
                if e isa ArgumentError
                    # if the memory variables trigger an error (e.g. K <= 0), fallback to zero
                    analytical_initial_map[cas_norm] = 0.0
                    analytical_C_bg_map[cas_norm] = 0.0
                else
                    rethrow(e)
                end
            end
        end

        # Pre-apply analytical burden to states before starting Month 1
        for sp in species_list
            state = species_states[sp]
            for (cas_norm, B_init) in analytical_initial_map
                TwoTimescaleResilience.set_internal_burden!(state, cas_norm, B_init)
            end
        end

        # Iterate over reported months
        for month in 1:12
            for sp in species_list
                state = species_states[sp]
                params = params_map[sp]

                for c in compounds
                    record = records_map[c.cas]
                    cas_norm = record["cas_norm"]
                    conc = C_t[month] # same concentration scenario for all compounds to compare memory

                    # Capture initial burden strictly AT the start of reported month 1 (before applying C_t[month])
                    initial_B_t_at_month_1 = analytical_initial_map[cas_norm]

                    # The stateful burden expects concs mapped by cas_norm, not cas hyphenated
                    concs = Dict(cas_norm => conc)

                    # update state and compute burden for this compound
                    stateful_burden = TwoTimescaleResilience.ecotox_records_to_deb_burden_stateful!(
                        state, concs, [record]; memory_library=memory
                    )

                    B_t = TwoTimescaleResilience.get_internal_burden(state, c.cas)
                    x_t = stateful_burden.maintenance # active stress is predominantly maintenance

                    push!(compound_results, (
                        scenario = scenario,
                        mixture_method = mixture_method,
                        spinup_used = spinup_used,
                        spinup_months = spinup_months,
                        spinup_method = spinup_method,
                        species_key = replace(sp, " " => "_"),
                        species_name = sp,
                        month = month,
                        cas_norm = cas_norm,
                        cas_hyphenated = record["cas_hyphenated"],
                        chemical_name = c.name,
                        C_t = conc,
                        B_t = B_t,
                        x_t = x_t,
                        burden_assimilation = stateful_burden.assimilation,
                        burden_maintenance = stateful_burden.maintenance,
                        burden_growth = stateful_burden.growth,
                        burden_reproduction = stateful_burden.reproduction,
                        spinup_background_C_t = analytical_C_bg_map[cas_norm],
                        initial_B_t_at_reported_month_1 = initial_B_t_at_month_1
                    ))
                end
            end
        end

        # Compute aggregated species response calculation for this scenario and mixture method
        for month in 1:12
            for sp in species_list
                params = params_map[sp]

                # Aggregate burdens across all compounds for this scenario, mixture method, species and month
                compounds_in_month = filter(r -> r.scenario == scenario && r.mixture_method == mixture_method && r.species_name == sp && r.month == month, compound_results)

                # Runtime analytical exact aggregation
                agg_burden = TwoTimescaleResilience.aggregate_deb_axis_burdens(compounds_in_month; mixture_method = mixture_method)

                diag = TwoTimescaleResilience.mixture_contribution_diagnostics(compounds_in_month)

                # For physiological response, map to standard NamedTuple expected by response functions
                mapped_burden = (
                    assimilation = agg_burden.total_burden_assimilation,
                    maintenance = agg_burden.total_burden_maintenance,
                    growth = agg_burden.total_burden_growth,
                    reproduction = agg_burden.total_burden_reproduction
                )

                # Compute physiological response
                resp = TwoTimescaleResilience.ecotox_burden_to_response(mapped_burden, params)

                # For baseline restoring force
                zero_burden = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
                zero_resp = TwoTimescaleResilience.ecotox_burden_to_response(zero_burden, params)

                push!(species_results, (
                    scenario = scenario,
                    mixture_method = mixture_method,
                    spinup_used = spinup_used,
                    spinup_months = spinup_months,
                    spinup_method = spinup_method,
                    species_key = replace(sp, " " => "_"),
                    species_name = sp,
                    month = month,

                    total_burden_assimilation = agg_burden.total_burden_assimilation,
                    total_burden_maintenance = agg_burden.total_burden_maintenance,
                    total_burden_growth = agg_burden.total_burden_growth,
                    total_burden_reproduction = agg_burden.total_burden_reproduction,

                    n_compounds_contributing_assimilation = diag.n_compounds_contributing_assimilation,
                    n_compounds_contributing_maintenance = diag.n_compounds_contributing_maintenance,
                    n_compounds_contributing_growth = diag.n_compounds_contributing_growth,
                    n_compounds_contributing_reproduction = diag.n_compounds_contributing_reproduction,

                    dominant_compound_assimilation = diag.dominant_compound_assimilation,
                    dominant_compound_maintenance = diag.dominant_compound_maintenance,
                    dominant_compound_growth = diag.dominant_compound_growth,
                    dominant_compound_reproduction = diag.dominant_compound_reproduction,

                    max_single_compound_fraction_assimilation = diag.max_single_compound_fraction_assimilation,
                    max_single_compound_fraction_maintenance = diag.max_single_compound_fraction_maintenance,
                    max_single_compound_fraction_growth = diag.max_single_compound_fraction_growth,
                    max_single_compound_fraction_reproduction = diag.max_single_compound_fraction_reproduction,

                    A_t = resp.A,
                    lambda_t = resp.lambda,
                    lambda0 = zero_resp.lambda,
                    F_t = resp.amplification
                ))
            end
        end
        end
    end

    println("Tranche 3 & 4 complete. Compound Rows: ", length(compound_results), " | Species Rows: ", length(species_results))

    # Tranche 5: CSV output
    output_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multispecies_multicompound_monthly_memory_demo"))
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Convert vectors of NamedTuples to DataFrames
    df_comp = DataFrame(compound_results)
    df_spec = DataFrame(species_results)

    compound_csv = joinpath(output_dir, "monthly_compound_summary.csv")
    species_csv = joinpath(output_dir, "monthly_species_summary.csv")

    CSV.write(compound_csv, df_comp)
    CSV.write(species_csv, df_spec)

    println("Tranche 5 complete. Wrote CSVs to ", output_dir)

    # Tranche 6: Diagnostic plots
    # Import CairoMakie inside the main script block conditionally if available or just assume it is used
    eval(:(using CairoMakie))

    # Basic Makie style setup
    theme = Theme(
        fontsize = 16,
        Axis = (xgridvisible = false, ygridvisible = false)
    )
    CairoMakie.set_theme!(theme)

    colors = [:black, :firebrick, :steelblue]
    linestyles = [:solid, :dash, :dot]

    # Plot 1: Concentrations
    # Because methods perfectly overlap, we use additive_axis_burden as the baseline for visualization
    df_comp_base = filter(r -> r.mixture_method == "additive_axis_burden", df_comp)
    df_spec_base = filter(r -> r.mixture_method == "additive_axis_burden", df_spec)

    fig_C = CairoMakie.Figure(size=(800, 600))
    ax_C = CairoMakie.Axis(fig_C[1,1], title="Ambient Concentrations (C_t)", xlabel="Month", ylabel="Concentration")
    for (i, cas) in enumerate(unique(df_comp_base.cas_norm))
        sub_df = filter(r -> r.cas_norm == cas && r.species_name == species_list[1] && r.scenario == "zero_start", df_comp_base) # C_t is same for all scenarios/species
        CairoMakie.lines!(ax_C, sub_df.month, sub_df.C_t, label=sub_df.chemical_name[1], color=colors[i], linewidth=2)
    end
    CairoMakie.axislegend(ax_C)
    CairoMakie.save(joinpath(output_dir, "monthly_concentrations.png"), fig_C)

    # Plot 2: Internal Burdens
    fig_B = CairoMakie.Figure(size=(1200, 400))
    for (s, sp) in enumerate(species_list)
        ax_B = CairoMakie.Axis(fig_B[1,s], title="Internal Burden - $sp", xlabel="Month", ylabel="Burden (B_t)")
        for (i, cas) in enumerate(unique(df_comp_base.cas_norm))
            sub_df_zero = filter(r -> r.cas_norm == cas && r.species_name == sp && r.scenario == "zero_start", df_comp_base)
            sub_df_warm = filter(r -> r.cas_norm == cas && r.species_name == sp && r.scenario == "analytical_warm_start", df_comp_base)
            CairoMakie.lines!(ax_B, sub_df_zero.month, sub_df_zero.B_t, label=sub_df_zero.chemical_name[1] * " (Zero)", color=colors[i], linewidth=2, linestyle=:solid)
            CairoMakie.lines!(ax_B, sub_df_warm.month, sub_df_warm.B_t, label=sub_df_warm.chemical_name[1] * " (Warm)", color=colors[i], linewidth=2, linestyle=:dash)
        end
        if s == 1
            CairoMakie.axislegend(ax_B)
        end
    end
    CairoMakie.save(joinpath(output_dir, "monthly_internal_burdens.png"), fig_B)

    # Plot 3: Active stress
    fig_x = CairoMakie.Figure(size=(1200, 400))
    for (s, sp) in enumerate(species_list)
        ax_x = CairoMakie.Axis(fig_x[1,s], title="Active Stress - $sp", xlabel="Month", ylabel="Stress (x_t)")
        for (i, cas) in enumerate(unique(df_comp_base.cas_norm))
            sub_df_zero = filter(r -> r.cas_norm == cas && r.species_name == sp && r.scenario == "zero_start", df_comp_base)
            sub_df_warm = filter(r -> r.cas_norm == cas && r.species_name == sp && r.scenario == "analytical_warm_start", df_comp_base)
            CairoMakie.lines!(ax_x, sub_df_zero.month, sub_df_zero.x_t, label=sub_df_zero.chemical_name[1] * " (Zero)", color=colors[i], linewidth=2, linestyle=:solid)
            CairoMakie.lines!(ax_x, sub_df_warm.month, sub_df_warm.x_t, label=sub_df_warm.chemical_name[1] * " (Warm)", color=colors[i], linewidth=2, linestyle=:dash)
        end
        if s == 1
            CairoMakie.axislegend(ax_x)
        end
    end
    CairoMakie.save(joinpath(output_dir, "monthly_active_stress.png"), fig_x)

    # Plot 4: Axis Burdens (Aggregated)
    fig_axis = CairoMakie.Figure(size=(1200, 400))
    for (s, sp) in enumerate(species_list)
        ax_axis = CairoMakie.Axis(fig_axis[1,s], title="Axis Burdens - $sp", xlabel="Month", ylabel="Aggregated Burden")

        agg_maint_zero = [sum(r.burden_maintenance for r in compound_results if r.mixture_method == "additive_axis_burden" && r.species_name == sp && r.month == m && r.scenario == "zero_start") for m in 1:12]
        agg_maint_warm = [sum(r.burden_maintenance for r in compound_results if r.mixture_method == "additive_axis_burden" && r.species_name == sp && r.month == m && r.scenario == "analytical_warm_start") for m in 1:12]

        CairoMakie.lines!(ax_axis, 1:12, agg_maint_zero, label="Maintenance (Zero)", color=:red, linewidth=2, linestyle=:solid)
        CairoMakie.lines!(ax_axis, 1:12, agg_maint_warm, label="Maintenance (Warm)", color=:red, linewidth=2, linestyle=:dash)
        if s == 1
            CairoMakie.axislegend(ax_axis)
        end
    end
    CairoMakie.save(joinpath(output_dir, "monthly_axis_burdens.png"), fig_axis)

    # Plot 5: Adaptive Margin
    fig_A = CairoMakie.Figure(size=(800, 600))
    ax_A = CairoMakie.Axis(fig_A[1,1], title="Adaptive Margin (A_t)", xlabel="Month", ylabel="Adaptive Margin")
    for (s, sp) in enumerate(species_list)
        sub_df_zero = filter(r -> r.species_name == sp && r.scenario == "zero_start", df_spec_base)
        sub_df_warm = filter(r -> r.species_name == sp && r.scenario == "analytical_warm_start", df_spec_base)
        CairoMakie.lines!(ax_A, sub_df_zero.month, sub_df_zero.A_t, label=sp * " (Zero)", color=colors[s], linewidth=2, linestyle=:solid)
        CairoMakie.lines!(ax_A, sub_df_warm.month, sub_df_warm.A_t, label=sp * " (Warm)", color=colors[s], linewidth=2, linestyle=:dash)
    end
    CairoMakie.axislegend(ax_A)
    CairoMakie.save(joinpath(output_dir, "monthly_adaptive_margin.png"), fig_A)

    # Plot 6: Restoring Force
    fig_lam = CairoMakie.Figure(size=(800, 600))
    ax_lam = CairoMakie.Axis(fig_lam[1,1], title="Restoring Force (lambda_t)", xlabel="Month", ylabel="Restoring Force")
    for (s, sp) in enumerate(species_list)
        sub_df_zero = filter(r -> r.species_name == sp && r.scenario == "zero_start", df_spec_base)
        sub_df_warm = filter(r -> r.species_name == sp && r.scenario == "analytical_warm_start", df_spec_base)
        CairoMakie.lines!(ax_lam, sub_df_zero.month, sub_df_zero.lambda_t, label=sp * " (Zero)", color=colors[s], linewidth=2, linestyle=:solid)
        CairoMakie.lines!(ax_lam, sub_df_warm.month, sub_df_warm.lambda_t, label=sp * " (Warm)", color=colors[s], linewidth=2, linestyle=:dash)
    end
    CairoMakie.axislegend(ax_lam)
    CairoMakie.save(joinpath(output_dir, "monthly_restoring_force.png"), fig_lam)

    # Plot 7: Amplification
    fig_F = CairoMakie.Figure(size=(800, 600))
    ax_F = CairoMakie.Axis(fig_F[1,1], title="Amplification (F_t)", xlabel="Month", ylabel="Amplification (Baseline=1.0)")
    for (s, sp) in enumerate(species_list)
        sub_df_zero = filter(r -> r.species_name == sp && r.scenario == "zero_start", df_spec_base)
        sub_df_warm = filter(r -> r.species_name == sp && r.scenario == "analytical_warm_start", df_spec_base)
        CairoMakie.lines!(ax_F, sub_df_zero.month, sub_df_zero.F_t, label=sp * " (Zero)", color=colors[s], linewidth=2, linestyle=:solid)
        CairoMakie.lines!(ax_F, sub_df_warm.month, sub_df_warm.F_t, label=sp * " (Warm)", color=colors[s], linewidth=2, linestyle=:dash)
    end
    CairoMakie.axislegend(ax_F)
    CairoMakie.save(joinpath(output_dir, "monthly_amplification.png"), fig_F)

    println("Tranche 6 complete. Wrote PNGs to ", output_dir)
end


main()

