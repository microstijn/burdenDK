using TwoTimescaleResilience
using Printf
using CSV
using DataFrames
using CairoMakie

"""
3 species × 3 compounds × 12 months monthly memory diagnostic.

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

    # Tranche 3: Compute stateful monthly burden
    # Create empty states for each species
    species_states = Dict(sp => EcotoxExposureState() for sp in species_list)

    # Store compound-level results
    compound_results = []

    # Iterate over scenario
    ec50_ref = records_map["7440-43-9"]["EC50_median"] # pedogogical scalar based on Cd EC50
    C_t = generate_scenario_concentrations(ec50_ref)

    for month in 1:12
        for sp in species_list
            state = species_states[sp]
            params = params_map[sp]

            for c in compounds
                record = records_map[c.cas]
                conc = C_t[month] # same concentration scenario for all compounds to compare memory

                # The stateful burden expects concs mapped by cas_norm, not cas hyphenated
                concs = Dict(record["cas_norm"] => conc)

                # update state and compute burden for this compound
                stateful_burden = ecotox_records_to_deb_burden_stateful!(
                    state, concs, [record]; memory_library=memory
                )

                B_t = get_internal_burden(state, c.cas)
                x_t = stateful_burden.maintenance # active stress is predominantly maintenance

                push!(compound_results, (
                    species_name = sp,
                    month = month,
                    cas_norm = record["cas_norm"],
                    chemical_name = c.name,
                    C_t = conc,
                    B_t = B_t,
                    x_t = x_t,
                    burden_assimilation = stateful_burden.assimilation,
                    burden_maintenance = stateful_burden.maintenance,
                    burden_growth = stateful_burden.growth,
                    burden_reproduction = stateful_burden.reproduction
                ))
            end
        end
    end

    println("Tranche 3 complete. Rows: ", length(compound_results))

    # Tranche 4: Compute aggregated species response calculation
    species_results = []

    for month in 1:12
        for sp in species_list
            params = params_map[sp]

            # Aggregate burdens across all compounds for this species and month
            compounds_in_month = filter(r -> r.species_name == sp && r.month == month, compound_results)

            agg_assimilation = sum(r.burden_assimilation for r in compounds_in_month)
            agg_maintenance = sum(r.burden_maintenance for r in compounds_in_month)
            agg_growth = sum(r.burden_growth for r in compounds_in_month)
            agg_reproduction = sum(r.burden_reproduction for r in compounds_in_month)

            agg_burden = (
                assimilation = agg_assimilation,
                maintenance = agg_maintenance,
                growth = agg_growth,
                reproduction = agg_reproduction
            )

            # Compute physiological response
            resp = ecotox_burden_to_response(agg_burden, params)

            # For baseline restoring force
            zero_burden = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
            zero_resp = ecotox_burden_to_response(zero_burden, params)

            push!(species_results, (
                species_key = replace(sp, " " => "_"),
                species_name = sp,
                month = month,
                A_t = resp.A,
                lambda_t = resp.lambda,
                lambda0 = zero_resp.lambda,
                F_t = resp.amplification
            ))
        end
    end

    println("Tranche 4 complete. Species-month Rows: ", length(species_results))

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
    fig_C = CairoMakie.Figure(size=(800, 600))
    ax_C = CairoMakie.Axis(fig_C[1,1], title="Ambient Concentrations (C_t)", xlabel="Month", ylabel="Concentration")
    for (i, cas) in enumerate(unique(df_comp.cas_norm))
        sub_df = filter(r -> r.cas_norm == cas && r.species_name == species_list[1], df_comp) # C_t is same for all species
        CairoMakie.lines!(ax_C, sub_df.month, sub_df.C_t, label=sub_df.chemical_name[1], color=colors[i], linewidth=2)
    end
    CairoMakie.axislegend(ax_C)
    CairoMakie.save(joinpath(output_dir, "monthly_concentrations.png"), fig_C)

    # Plot 2: Internal Burdens
    fig_B = CairoMakie.Figure(size=(1200, 400))
    for (s, sp) in enumerate(species_list)
        ax_B = CairoMakie.Axis(fig_B[1,s], title="Internal Burden - $sp", xlabel="Month", ylabel="Burden (B_t)")
        for (i, cas) in enumerate(unique(df_comp.cas_norm))
            sub_df = filter(r -> r.cas_norm == cas && r.species_name == sp, df_comp)
            CairoMakie.lines!(ax_B, sub_df.month, sub_df.B_t, label=sub_df.chemical_name[1], color=colors[i], linewidth=2)
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
        for (i, cas) in enumerate(unique(df_comp.cas_norm))
            sub_df = filter(r -> r.cas_norm == cas && r.species_name == sp, df_comp)
            CairoMakie.lines!(ax_x, sub_df.month, sub_df.x_t, label=sub_df.chemical_name[1], color=colors[i], linewidth=2)
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
        sub_df = filter(r -> r.species_name == sp, df_spec)
        # Note: we need to pull the axis burdens from compound_results by summing them, but wait, those were aggregated internally. Let's just re-aggregate them here for plotting, or just plot the primary stress driver.
        # Actually, let's just re-sum them here since they aren't in df_spec
        agg_maint = [sum(r.burden_maintenance for r in compound_results if r.species_name == sp && r.month == m) for m in 1:12]
        CairoMakie.lines!(ax_axis, 1:12, agg_maint, label="Maintenance", color=:red, linewidth=2)
        if s == 1
            CairoMakie.axislegend(ax_axis)
        end
    end
    CairoMakie.save(joinpath(output_dir, "monthly_axis_burdens.png"), fig_axis)

    # Plot 5: Adaptive Margin
    fig_A = CairoMakie.Figure(size=(800, 600))
    ax_A = CairoMakie.Axis(fig_A[1,1], title="Adaptive Margin (A_t)", xlabel="Month", ylabel="Adaptive Margin")
    for (s, sp) in enumerate(species_list)
        sub_df = filter(r -> r.species_name == sp, df_spec)
        CairoMakie.lines!(ax_A, sub_df.month, sub_df.A_t, label=sp, color=colors[s], linewidth=2, linestyle=linestyles[s])
    end
    CairoMakie.axislegend(ax_A)
    CairoMakie.save(joinpath(output_dir, "monthly_adaptive_margin.png"), fig_A)

    # Plot 6: Restoring Force
    fig_lam = CairoMakie.Figure(size=(800, 600))
    ax_lam = CairoMakie.Axis(fig_lam[1,1], title="Restoring Force (lambda_t)", xlabel="Month", ylabel="Restoring Force")
    for (s, sp) in enumerate(species_list)
        sub_df = filter(r -> r.species_name == sp, df_spec)
        CairoMakie.lines!(ax_lam, sub_df.month, sub_df.lambda_t, label=sp, color=colors[s], linewidth=2, linestyle=linestyles[s])
    end
    CairoMakie.axislegend(ax_lam)
    CairoMakie.save(joinpath(output_dir, "monthly_restoring_force.png"), fig_lam)

    # Plot 7: Amplification
    fig_F = CairoMakie.Figure(size=(800, 600))
    ax_F = CairoMakie.Axis(fig_F[1,1], title="Amplification (F_t)", xlabel="Month", ylabel="Amplification (Baseline=1.0)")
    for (s, sp) in enumerate(species_list)
        sub_df = filter(r -> r.species_name == sp, df_spec)
        CairoMakie.lines!(ax_F, sub_df.month, sub_df.F_t, label=sp, color=colors[s], linewidth=2, linestyle=linestyles[s])
    end
    CairoMakie.axislegend(ax_F)
    CairoMakie.save(joinpath(output_dir, "monthly_amplification.png"), fig_F)

    println("Tranche 6 complete. Wrote PNGs to ", output_dir)
end


main()

