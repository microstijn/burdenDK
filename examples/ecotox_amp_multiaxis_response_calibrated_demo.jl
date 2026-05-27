# examples/ecotox_amp_multiaxis_response_calibrated_demo.jl
#
# PURPOSE:
# This script demonstrates the response-calibrated multi-axis diagnostic.
# It uses existing runtime model components to show how different pressure
# geometries narrow adaptive margin and alter amplification.
#
# SCOPE & NON-GOALS:
# - Uses existing ECOTOX records only (no fabricated records).
# - Retains exact same-axis mixture aggregation.
# - Both `axis_toxic_unit_sum` and `additive_axis_burden` remain identical here.
# - No synergism or antagonism implemented.
# - No Independent Action implemented.
# - No physiological condition memory Z_t implemented.
# - No raster integration implemented.
# - DEB equations are NOT altered.
#
# CALIBRATION LOGIC:
# We select records that activate Maintenance, Growth, and Reproduction axes.
# A deterministic synthetic exposure scenario is generated to expose species to
# single-axis and multi-axis pressure pulses. Calibration targets bounded peak F_t
# values without forcing complete collapse.

using TwoTimescaleResilience
using CSV
using DataFrames
using JSON
using CairoMakie

function inspect_candidate_records()
    library = load_ecotox_library()

    candidates = []

    for record in library
        # Basic validity check
        try
            validate_ecotox_record(record)
        catch
            continue
        end

        # Check NOEC and EC50 validity for this context
        NOEC_median = record["NOEC_median"]
        EC50_median = record["EC50_median"]

        if ismissing(NOEC_median) || ismissing(EC50_median) ||
           NOEC_median < 0 || EC50_median <= NOEC_median
            continue
        end

        cas_norm = record["cas_norm"]
        if isempty(cas_norm)
            continue
        end

        try
            axis = ecotox_effect_to_deb_axis(record["effect_code"])
            push!(candidates, (
                axis = string(axis),
                effect_code = record["effect_code"],
                cas_norm = cas_norm,
                cas_hyphenated = record["cas_hyphenated"],
                chemical_name = get(record, "chemical_name", "Unknown"),
                taxon_class = record["taxon_class"],
                NOEC_median = NOEC_median,
                EC50_median = EC50_median,
                n_NOEC = record["n_NOEC"],
                n_EC50 = record["n_EC50"]
            ))
        catch
            continue
        end
    end

    out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
    mkpath(out_dir)

    df = DataFrame(candidates)

    # Sort by axis and then effect code for clean output
    sort!(df, [:axis, :effect_code, :cas_norm])

    out_path = joinpath(out_dir, "ecotox_axis_candidate_records.csv")
    CSV.write(out_path, df)

    println("Candidate records summary written to: ", out_path)

    # internal summary for stdout
    for ax in ["assimilation", "maintenance", "growth", "reproduction"]
        subset = filter(row -> row.axis == ax, df)
        println("Axis '$ax': $(nrow(subset)) valid candidate records")
    end

    return df
end

function select_species_and_records(candidates_df)
    amp_lib = load_amp_species_library()
    ecotox_lib = load_ecotox_library()

    # 1. Select Species
    species_keys = ["Abatus_cordatus", "Podarcis_muralis", "Thalia_democratica"]
    selected_species = []

    for key in species_keys
        try
            params = amp_species_deb_params(amp_lib, key)
            push!(selected_species, (key = key, params = params))
        catch e
            println("Warning: Could not load DEBAxisParams for $key: $e")
        end
    end

    println("\nSelected $(length(selected_species)) species with distinct capacity:")
    for sp in selected_species
        println("  $(sp.key) (A0 = $(sp.params.A0), lambda_max = $(sp.params.lambda_max))")
    end

    # 2. Select ECOTOX Records
    # Let's pick 3 records to target Maintenance, Growth, and Reproduction.
    # We will pick specific CAS and effect code combos that exist in candidates_df
    # We prefer records with decent n_NOEC / n_EC50 if possible, and different compounds to make it interesting.

    selected_record_specs = [
        # 1. Maintenance
        (cas = "7440-43-9", effect = "MOR"), # Cadmium, often used in examples

        # 2. Growth
        (cas = "7440-66-6", effect = "DVP"), # Zinc

        # 3. Reproduction
        (cas = "108-95-2", effect = "REP"),  # Phenol
    ]

    selected_records = []

    for spec in selected_record_specs
        # Find exact match in library using the filtering runtime function
        # We also want to make sure it matches the NOEC/EC50 from the candidate df to guarantee validity
        matches = ecotox_filter_records(ecotox_lib; cas=spec.cas, effect_code=spec.effect)

        for m in matches
            # check if it's our valid candidate
            if m["EC50_median"] > m["NOEC_median"]
                push!(selected_records, m)
                break
            end
        end
    end

    println("\nSelected $(length(selected_records)) ECOTOX records for multi-axis scenario:")
    for rec in selected_records
        ax = ecotox_effect_to_deb_axis(rec["effect_code"])
        println("  CAS: $(rec["cas_hyphenated"]), Effect: $(rec["effect_code"]) -> Axis: $ax")
    end

    return selected_species, selected_records
end

function define_calibrated_scenario(selected_records)
    # Define a 12-month scenario structure
    # month 1-2: no pressure
    # month 3-4: maintenance pulse
    # month 5-6: no pressure
    # month 7-8: multi-axis pulse (maintenance + growth + reproduction)
    # month 9-12: recovery (no added pressure, memory carryover)

    # Calibration k-values for C_t = NOEC + k * (EC50 - NOEC)
    # Start with k=1.0 and adjust to hit F_t ≈ 1.05-1.5 for at least one species
    # If all fail to hit it without collapse, we can bump k. Let's use k=5.0 for the multi-axis pulse
    # to guarantee some response but not immediate total collapse.

    # Using deterministic k values
    k_single = 2.0
    k_multi = 4.0

    # We will build a matrix or array of dicts for each month
    scenario = []

    for month in 1:12
        month_concs = Dict{String, Float64}()

        # Default all selected compounds to 0
        for rec in selected_records
            month_concs[rec["cas_norm"]] = 0.0
        end

        if month in 3:4
            # Single pulse on maintenance
            for rec in selected_records
                if rec["effect_code"] == "MOR"
                    NOEC = rec["NOEC_median"]
                    EC50 = rec["EC50_median"]
                    month_concs[rec["cas_norm"]] = NOEC + k_single * (EC50 - NOEC)
                end
            end
        elseif month in 7:8
            # Multi-axis pulse
            for rec in selected_records
                NOEC = rec["NOEC_median"]
                EC50 = rec["EC50_median"]
                month_concs[rec["cas_norm"]] = NOEC + k_multi * (EC50 - NOEC)
            end
        end

        push!(scenario, month_concs)
    end

    return scenario, (k_single=k_single, k_multi=k_multi)
end

function run_diagnostic_scenario(selected_species, selected_records, scenario_months)
    memory_lib = load_compound_memory_library()

    # We will track results for both mixture methods, though they should be identical
    mixture_methods = ["additive_axis_burden", "axis_toxic_unit_sum"]
    response_modes = ["raw_margin_subtraction", "ec50_anchored_fractional_impairment"]

    compound_results = []
    species_results = []

    for method in mixture_methods
        for r_mode in response_modes
            for sp in selected_species
            sp_key = sp.key
            sp_params = sp.params

            state = TwoTimescaleResilience.EcotoxExposureState()

            for (month_idx, month_concs) in enumerate(scenario_months)
                # 1. Update memory and compute individual burdens
                monthly_compound_burdens = []

                for rec in selected_records
                    cas_norm = rec["cas_norm"]
                    C = get(month_concs, cas_norm, 0.0)

                    rho = compound_retention(cas_norm; memory_library=memory_lib)
                    K_bio = compound_bioaccumulation_factor(cas_norm; memory_library=memory_lib)

                    # Update internal burden state
                    B = update_internal_burden!(state, cas_norm, C; retention=rho, bioaccumulation_factor=K_bio, memory_library=memory_lib)

                    # Calculate axis-specific burden
                    burden = ecotox_record_to_deb_burden(B, rec)

                    # Active stress unit x_t
                    NOEC = rec["NOEC_median"]
                    EC50 = rec["EC50_median"]
                    x_t = TwoTimescaleResilience.ecotox_active_stress(B, NOEC, EC50)

                    chem_name = get(rec, "chemical_name", "Unknown")
                    ax = TwoTimescaleResilience.ecotox_effect_to_deb_axis(rec["effect_code"])

                    push!(monthly_compound_burdens, (
                        chemical_name = chem_name,
                        burden_assimilation = burden.assimilation,
                        burden_maintenance = burden.maintenance,
                        burden_growth = burden.growth,
                        burden_reproduction = burden.reproduction
                    ))

                    push!(compound_results, (
                        scenario = "response_calibrated_multiaxis",
                        mixture_method = method,
                        response_mode = r_mode,
                        species_key = sp_key,
                        species_name = replace(sp_key, "_" => " "),
                        month = month_idx,
                        axis = string(ax),
                        effect_code = rec["effect_code"],
                        cas_norm = cas_norm,
                        cas_hyphenated = rec["cas_hyphenated"],
                        chemical_name = chem_name,
                        C_t = C,
                        B_t = B,
                        x_t = x_t,
                        burden_assimilation = burden.assimilation,
                        burden_maintenance = burden.maintenance,
                        burden_growth = burden.growth,
                        burden_reproduction = burden.reproduction
                    ))
                end

                # 2. Aggregate burdens using the runtime function
                agg = aggregate_deb_axis_burdens(monthly_compound_burdens; mixture_method=method)

                # 3. Compute contribution diagnostics
                diag = mixture_contribution_diagnostics(monthly_compound_burdens)

                # 4. Compute physiological response
                mapped_agg = (
                    assimilation = agg.total_burden_assimilation,
                    maintenance = agg.total_burden_maintenance,
                    growth = agg.total_burden_growth,
                    reproduction = agg.total_burden_reproduction
                )
                response = TwoTimescaleResilience.compute_adaptive_margin_response(mapped_agg, sp_params, response_mode=r_mode)

                push!(species_results, (
                    scenario = "response_calibrated_multiaxis",
                    mixture_method = method,
                    response_mode = r_mode,
                    species_key = sp_key,
                    species_name = replace(sp_key, "_" => " "),
                    month = month_idx,

                    total_burden_assimilation = agg.total_burden_assimilation,
                    total_burden_maintenance = agg.total_burden_maintenance,
                    total_burden_growth = agg.total_burden_growth,
                    total_burden_reproduction = agg.total_burden_reproduction,

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

                    axis_weight_method = response.axis_weight_method,
                    axis_weight_scope = response.axis_weight_scope,

                    X_assimilation = response.X_assimilation,
                    X_maintenance = response.X_maintenance,
                    X_growth = response.X_growth,
                    X_reproduction = response.X_reproduction,
                    
                    E_assimilation = response.E_assimilation,
                    E_maintenance = response.E_maintenance,
                    E_growth = response.E_growth,
                    E_reproduction = response.E_reproduction,
                    
                    w_assimilation = response.w_assimilation,
                    w_maintenance = response.w_maintenance,
                    w_growth = response.w_growth,
                    w_reproduction = response.w_reproduction,
                    
                    Q_t = response.Q_t,
                    A0 = response.A0,
                    A_t = response.A,
                    lambda0 = response.lambda0,
                    lambda_t = response.lambda,
                    F_t = response.amplification
                ))
            end
        end
    end
    end

    out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
    mkpath(out_dir)

    df_comp = DataFrame(compound_results)
    df_spec = DataFrame(species_results)

    CSV.write(joinpath(out_dir, "multiaxis_compound_summary.csv"), df_comp)
    CSV.write(joinpath(out_dir, "multiaxis_species_summary.csv"), df_spec)

    return df_comp, df_spec
end

function perform_calibration_checks(df_spec)
    out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))

    # We only need to check one baseline mixture method
    df = filter(row -> row.mixture_method == "additive_axis_burden", df_spec)
    
    # Keeping old calibration summary logic for backwards compatibility tests, but scoping to raw mode if possible
    df_raw = filter(row -> row.response_mode == "raw_margin_subtraction", df)

    species_keys = unique(df_raw.species_key)

    calibration_summary = []
    target_hit = false

    for sp in species_keys
        sp_df = filter(row -> row.species_key == sp, df_raw)

        max_F_t_idx = argmax(sp_df.F_t)
        max_F_t = sp_df.F_t[max_F_t_idx]
        month_max_F_t = sp_df.month[max_F_t_idx]

        min_A_t_idx = argmin(sp_df.A_t)
        min_A_t = sp_df.A_t[min_A_t_idx]
        month_min_A_t = sp_df.month[min_A_t_idx]

        min_lambda_t_idx = argmin(sp_df.lambda_t)
        min_lambda_t = sp_df.lambda_t[min_lambda_t_idx]
        month_min_lambda_t = sp_df.month[min_lambda_t_idx]

        active_axes = Set{String}()
        for row in eachrow(sp_df)
            if row.total_burden_assimilation > 0; push!(active_axes, "assimilation"); end
            if row.total_burden_maintenance > 0; push!(active_axes, "maintenance"); end
            if row.total_burden_growth > 0; push!(active_axes, "growth"); end
            if row.total_burden_reproduction > 0; push!(active_axes, "reproduction"); end
        end

        notes = []
        if max_F_t >= 1.05
            push!(notes, "Target F_t >= 1.05 achieved")
            target_hit = true
        end
        if max_F_t > 1.5
            push!(notes, "F_t > 1.5 (High response)")
        end
        if min_A_t <= 0
            push!(notes, "Collapse (A_t <= 0)")
        end

        push!(calibration_summary, (
            species_key = sp,
            species_name = replace(sp, "_" => " "),
            peak_F_t = max_F_t,
            month_peak_F_t = month_max_F_t,
            min_A_t = min_A_t,
            month_min_A_t = month_min_A_t,
            min_lambda_t = min_lambda_t,
            month_min_lambda_t = month_min_lambda_t,
            active_axes_count = length(active_axes),
            activated_axes = join(sort(collect(active_axes)), ", "),
            notes = isempty(notes) ? "Target F_t not reached" : join(notes, "; ")
        ))
    end

    cal_df = DataFrame(calibration_summary)
    CSV.write(joinpath(out_dir, "multiaxis_response_calibration_summary.csv"), cal_df)
    
    # -------------------------------------------------------------
    # Tranche 5 logic: New Comparison Summary
    # -------------------------------------------------------------
    comparison_summary = []
    
    response_modes = unique(df.response_mode)
    
    for sp in species_keys
        sp_df_all = filter(row -> row.species_key == sp, df)
        
        # Raw results to compare against
        raw_df = filter(row -> row.response_mode == "raw_margin_subtraction", sp_df_all)
        max_F_t_raw = maximum(raw_df.F_t)
        min_A_t_raw = minimum(raw_df.A_t)
        
        for r_mode in response_modes
            sp_mode_df = filter(row -> row.response_mode == r_mode, sp_df_all)
            
            max_F_t_idx = argmax(sp_mode_df.F_t)
            max_F_t = sp_mode_df.F_t[max_F_t_idx]
            month_max_F_t = sp_mode_df.month[max_F_t_idx]
            
            min_A_t_idx = argmin(sp_mode_df.A_t)
            min_A_t = sp_mode_df.A_t[min_A_t_idx]
            month_min_A_t = sp_mode_df.month[min_A_t_idx]
            
            max_Q_t_idx = argmax(sp_mode_df.Q_t)
            max_Q_t = sp_mode_df.Q_t[max_Q_t_idx]
            month_max_Q_t = sp_mode_df.month[max_Q_t_idx]
            
            active_axes = Set{String}()
            for row in eachrow(sp_mode_df)
                if row.total_burden_assimilation > 0; push!(active_axes, "assimilation"); end
                if row.total_burden_maintenance > 0; push!(active_axes, "maintenance"); end
                if row.total_burden_growth > 0; push!(active_axes, "growth"); end
                if row.total_burden_reproduction > 0; push!(active_axes, "reproduction"); end
            end
            
            # Compare to raw mode
            if r_mode == "ec50_anchored_fractional_impairment"
                delta_F = max_F_t - max_F_t_raw
                delta_A = min_A_t - min_A_t_raw
            else
                delta_F = 0.0
                delta_A = 0.0
            end
            
            push!(comparison_summary, (
                species_key = sp,
                response_mode = r_mode,
                max_F_t = max_F_t,
                month_max_F_t = month_max_F_t,
                min_A_t = min_A_t,
                month_min_A_t = month_min_A_t,
                max_Q_t = max_Q_t,
                month_max_Q_t = month_max_Q_t,
                activated_axes = join(sort(collect(active_axes)), ", "),
                axis_weight_method = sp_mode_df.axis_weight_method[1],
                axis_weight_scope = sp_mode_df.axis_weight_scope[1],
                delta_max_F_t_ec50_minus_raw = delta_F,
                delta_min_A_t_ec50_minus_raw = delta_A
            ))
        end
    end
    
    comp_sum_df = DataFrame(comparison_summary)
    CSV.write(joinpath(out_dir, "multiaxis_response_mode_comparison_summary.csv"), comp_sum_df)

    println("\nCalibration Checks (Raw Mode):")
    for row in eachrow(cal_df)
        println("  $(row.species_key): Peak F_t = $(round(row.peak_F_t, digits=3)) (Month $(row.month_peak_F_t)), Min A_t = $(round(row.min_A_t, digits=1)), Notes: $(row.notes)")
    end

    return cal_df, target_hit
end

function generate_plots(df_spec, df_comp)
    out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))

    # We plot the baseline mixture method since they are expected to overlap
    df = filter(row -> row.mixture_method == "additive_axis_burden", df_spec)

    species_keys = unique(df.species_key)
    colors = cgrad(:tab10, length(species_keys), categorical=true)

    # Styling configuration
    ec50_linestyle = nothing      # solid
    raw_linestyle = :dash         # dashed
    marker_ec50 = :rect           # square
    marker_raw = :circle          # circle

    # Plot A: Multi-axis Axis Burdens
    # Filter to one response mode because burdens are independent of response mode
    df_ec50 = filter(row -> row.response_mode == "ec50_anchored_fractional_impairment", df)
    # We pick the first species to show the axis breakdown, as stress/burden is identical across species
    sp_show = species_keys[1]
    sp_df_ec50 = filter(row -> row.species_key == sp_show, df_ec50)

    fig_B = Figure(size=(800, 400))
    ax_B = Axis(fig_B[1, 1], xlabel="Month", ylabel="Total Axis Burden (s_a)", title="DEB Axis Burdens (Species: $(replace(sp_show, "_" => " ")))")

    lines!(ax_B, sp_df_ec50.month, sp_df_ec50.total_burden_assimilation, label="Assimilation", linewidth=2)
    lines!(ax_B, sp_df_ec50.month, sp_df_ec50.total_burden_maintenance, label="Maintenance", linewidth=2)
    lines!(ax_B, sp_df_ec50.month, sp_df_ec50.total_burden_growth, label="Growth", linewidth=2)
    lines!(ax_B, sp_df_ec50.month, sp_df_ec50.total_burden_reproduction, label="Reproduction", linewidth=2)

    axislegend(ax_B, position=:lt)
    save(joinpath(out_dir, "multiaxis_axis_burdens.png"), fig_B)

    # Plot B: Multi-axis Axis Impairments
    # Only meaningful for EC50 mode
    fig_I = Figure(size=(800, 400))
    ax_I = Axis(fig_I[1, 1], xlabel="Month", ylabel="Fractional Impairment (E_a)", title="Axis Impairments (EC50 Mode, Species: $(replace(sp_show, "_" => " ")))")
    
    lines!(ax_I, sp_df_ec50.month, sp_df_ec50.E_assimilation, label="E_assimilation", linewidth=2)
    lines!(ax_I, sp_df_ec50.month, sp_df_ec50.E_maintenance, label="E_maintenance", linewidth=2)
    lines!(ax_I, sp_df_ec50.month, sp_df_ec50.E_growth, label="E_growth", linewidth=2)
    lines!(ax_I, sp_df_ec50.month, sp_df_ec50.E_reproduction, label="E_reproduction", linewidth=2)

    axislegend(ax_I, position=:lt)
    save(joinpath(out_dir, "multiaxis_axis_impairments.png"), fig_I)

    # Plot C: Multi-axis Weighted Impairment Q_t
    # Only meaningful for EC50 mode
    fig_Q = Figure(size=(800, 400))
    ax_Q = Axis(fig_Q[1, 1], xlabel="Month", ylabel="Weighted Impairment (Q_t)", title="Weighted Impairment (Q_t) by Species (EC50 Mode)")
    for (i, sp) in enumerate(species_keys)
        sp_df = filter(row -> row.species_key == sp, df_ec50)
        lines!(ax_Q, sp_df.month, sp_df.Q_t, label=replace(sp, "_" => " "), color=colors[i], linewidth=2)
    end
    axislegend(ax_Q, position=:lt)
    save(joinpath(out_dir, "multiaxis_weighted_impairment_Q.png"), fig_Q)

    # Plot D: Multi-axis Adaptive Margin
    fig_A = Figure(size=(800, 400))
    ax_A = Axis(fig_A[1, 1], xlabel="Month", ylabel="Adaptive Margin (A_t)", title="Adaptive Margin Response by Species")
    for (i, sp) in enumerate(species_keys)
        sp_df_ec50 = filter(row -> row.species_key == sp && row.response_mode == "ec50_anchored_fractional_impairment", df)
        sp_df_raw = filter(row -> row.species_key == sp && row.response_mode == "raw_margin_subtraction", df)
        
        name = replace(sp, "_" => " ")
        lines!(ax_A, sp_df_ec50.month, sp_df_ec50.A_t, label="$name — EC50", color=colors[i], linestyle=ec50_linestyle, linewidth=2)
        lines!(ax_A, sp_df_raw.month, sp_df_raw.A_t, label="$name — Raw", color=colors[i], linestyle=raw_linestyle, linewidth=2)
    end
    axislegend(ax_A, position=:lt)
    save(joinpath(out_dir, "multiaxis_adaptive_margin.png"), fig_A)

    # Plot E: Multi-axis Restoring Force
    fig_L = Figure(size=(800, 400))
    ax_L = Axis(fig_L[1, 1], xlabel="Month", ylabel="Restoring Force (lambda_t)", title="Restoring Force Response by Species")
    for (i, sp) in enumerate(species_keys)
        sp_df_ec50 = filter(row -> row.species_key == sp && row.response_mode == "ec50_anchored_fractional_impairment", df)
        sp_df_raw = filter(row -> row.species_key == sp && row.response_mode == "raw_margin_subtraction", df)
        
        name = replace(sp, "_" => " ")
        lines!(ax_L, sp_df_ec50.month, sp_df_ec50.lambda_t, label="$name — EC50", color=colors[i], linestyle=ec50_linestyle, linewidth=2)
        lines!(ax_L, sp_df_raw.month, sp_df_raw.lambda_t, label="$name — Raw", color=colors[i], linestyle=raw_linestyle, linewidth=2)
    end
    axislegend(ax_L, position=:lt)
    save(joinpath(out_dir, "multiaxis_restoring_force.png"), fig_L)

    # Plot F: Multi-axis Amplification
    fig_F = Figure(size=(800, 400))
    ax_F = Axis(fig_F[1, 1], xlabel="Month", ylabel="Amplification (F_t)", title="Amplification Response by Species")
    for (i, sp) in enumerate(species_keys)
        sp_df_ec50 = filter(row -> row.species_key == sp && row.response_mode == "ec50_anchored_fractional_impairment", df)
        sp_df_raw = filter(row -> row.species_key == sp && row.response_mode == "raw_margin_subtraction", df)
        
        name = replace(sp, "_" => " ")
        lines!(ax_F, sp_df_ec50.month, sp_df_ec50.F_t, label="$name — EC50", color=colors[i], linestyle=ec50_linestyle, linewidth=2)
        lines!(ax_F, sp_df_raw.month, sp_df_raw.F_t, label="$name — Raw", color=colors[i], linestyle=raw_linestyle, linewidth=2)
    end
    axislegend(ax_F, position=:lt)
    save(joinpath(out_dir, "multiaxis_amplification.png"), fig_F)

    # Existing Dominant Compounds Scatter
    fig_D = Figure(size=(800, 600))
    ax_D = Axis(fig_D[1, 1], xlabel="Month", ylabel="Dominant Compound Fraction", title="Dominant Compound by Axis")

    # We plot the fraction for each axis using the ec50 df for the first species, as dominant fractions are independent of response mode
    sp_df_ec50_first = filter(row -> row.species_key == species_keys[1], df_ec50)
    
    axes_list = [
        (:maintenance, :max_single_compound_fraction_maintenance),
        (:growth, :max_single_compound_fraction_growth),
        (:reproduction, :max_single_compound_fraction_reproduction)
    ]

    for (ax_name, frac_col) in axes_list
        scatter!(ax_D, sp_df_ec50_first.month, sp_df_ec50_first[!, frac_col], label=string(ax_name), markersize=10)
    end
    axislegend(ax_D, position=:lt)
    save(joinpath(out_dir, "multiaxis_dominant_compounds.png"), fig_D)

    println("\nPlots generated in: ", out_dir)
end

function main()
    println("Running Response-Calibrated Multi-Axis Diagnostic...")

    # Tranche 1
    candidates_df = inspect_candidate_records()

    # Tranche 2
    selected_species, selected_records = select_species_and_records(candidates_df)

    # Tranche 3
    scenario_months, calibration_k = define_calibrated_scenario(selected_records)

    println("\nScenario Defined:")
    println("  k_single: $(calibration_k.k_single)")
    println("  k_multi: $(calibration_k.k_multi)")
    for (m, concs) in enumerate(scenario_months)
        println("  Month $m: ", concs)
    end

    # Tranche 4
    df_comp, df_spec = run_diagnostic_scenario(selected_species, selected_records, scenario_months)
    println("\nGenerated summaries:")
    println("  multiaxis_compound_summary.csv ($(nrow(df_comp)) rows)")
    println("  multiaxis_species_summary.csv ($(nrow(df_spec)) rows)")

    # Tranche 5
    cal_df, target_hit = perform_calibration_checks(df_spec)

    # Tranche 6
    generate_plots(df_spec, df_comp)
end


main()

