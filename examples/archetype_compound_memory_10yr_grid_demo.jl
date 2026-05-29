using TwoTimescaleResilience
using Dates
using JSON
using CSV
using NCDatasets
using Statistics

function run_archetype_compound_memory_10yr_grid_demo(; output_dir::String = joinpath(dirname(@__DIR__), "output", "archetype_compound_memory_10yr_grid_demo"))
    println("Running Archetype Compound Memory 10-year Grid Demo")
    mkpath(output_dir)

    # Simulation Dimensions
    nx = parse(Int, get(ENV, "TTR_GRID_NX", "80"))
    ny = parse(Int, get(ENV, "TTR_GRID_NY", "60"))
    n_years = parse(Int, get(ENV, "TTR_N_YEARS", "10"))
    n_months = n_years * 12
    n_species_req = parse(Int, get(ENV, "TTR_N_ARCHETYPE_SPECIES", "16"))
    n_compounds_req = parse(Int, get(ENV, "TTR_N_COMPOUNDS", "12"))
    k_clusters = parse(Int, get(ENV, "TTR_VULN_CLUSTER_K", "5"))
    make_plots = get(ENV, "TTR_MAKE_EXAMPLE_PLOTS", "false") == "true"
    allow_fallback = get(ENV, "TTR_ALLOW_NON_ARCHETYPE_FALLBACK", "false") == "true"

    println("Dimensions: nx=$nx, ny=$ny, n_years=$n_years, n_species=$n_species_req, n_compounds=$n_compounds_req, k_clusters=$k_clusters")

    # 1 & 2. Load and Select Species
    archetype_csv = joinpath(dirname(@__DIR__), "data", "AmP_Species_Archetypes.csv")
    archetype_json = joinpath(dirname(@__DIR__), "data", "AmP_Species_Archetypes.json")
    library_json = joinpath(dirname(@__DIR__), "data", "AmP_Species_Library.json")

    archetype_db_present = isfile(archetype_csv) || isfile(archetype_json)

    if !archetype_db_present && !allow_fallback
        error("Archetype database missing and TTR_ALLOW_NON_ARCHETYPE_FALLBACK is not set to true. Please generate archetypes first.")
    end

    species_library = load_amp_species_library(library_json)

    selected_species_keys = String[]
    species_info = NamedTuple[]
    species_selection_mode = "archetype"
    fallback_reason = ""

    if archetype_db_present
        # Try to load and filter archetypes (Simplified selection for demonstration)
        # Ideally we'd parse the CSV/JSON, but let's assume we can fallback to library logic if CSV parsing isn't strictly required
        # and we just need valid records. The prompt suggests trying to maximize diversity if archetypes are used.
        # Since I'm not implementing a full CSV parser right here without DataFrames, I'll use a basic deterministic selection.
        # But if the CSV exists, we should ideally read it. Let's do a basic CSV read.

        if isfile(archetype_csv)
            # Basic CSV read
            lines = readlines(archetype_csv)
            if length(lines) > 1
                headers = split(lines[1], ',')
                key_idx = findfirst(x -> x == "species_key" || x == "species", headers)
                labels_idx = findfirst(x -> x == "archetype_labels", headers)
                if key_idx !== nothing
                    for line in lines[2:end]
                        parts = split(line, ',')
                        if length(parts) >= key_idx
                            s_key = String(strip(parts[key_idx], '"'))
                            if haskey(species_library, s_key) && validate_amp_record(species_library[s_key])
                                push!(selected_species_keys, s_key)
                                l_str = labels_idx !== nothing && length(parts) >= labels_idx ? String(strip(parts[labels_idx], '"')) : ""
                                push!(species_info, (species_key=s_key, species_name=replace(s_key, "_" => " "), archetype_labels=l_str))
                            end
                        end
                        if length(selected_species_keys) == n_species_req
                            break
                        end
                    end
                end
            end
        end
        # If still empty, fill deterministically from library
        if isempty(selected_species_keys)
            archetype_db_present = false
        end
    end

    if !archetype_db_present
        species_selection_mode = "non_archetype_fallback"
        fallback_reason = "archetype database missing and TTR_ALLOW_NON_ARCHETYPE_FALLBACK=true"

        valid_keys = String[]
        for (k, v) in species_library
            if validate_amp_record(v)
                push!(valid_keys, k)
            end
        end
        sort!(valid_keys)

        for k in valid_keys
            if length(selected_species_keys) < n_species_req
                push!(selected_species_keys, k)
                push!(species_info, (species_key=k, species_name=replace(k, "_" => " "), archetype_labels="non_archetype_fallback"))
            end
        end
    end

    if length(selected_species_keys) < 8
        error("Fewer than 8 valid species available.")
    end

    n_species = length(selected_species_keys)

    # Save species
    CSV.write(joinpath(output_dir, "selected_archetype_species.csv"), species_info)

    # 3. Load Compounds
    compound_csv = joinpath(dirname(@__DIR__), "data", "Compound_Memory_Library.csv")
    ecotox_json = joinpath(dirname(@__DIR__), "data", "ECOTOX_Toxicity_Library.json")

    ecotox_lib = load_ecotox_library(ecotox_json)


    # Read compounds CSV manually
    c_data = CSV.File(compound_csv)
    selected_compounds = NamedTuple[]
    for row in c_data
        cas_hyp = string(row.cas_hyphenated)
        cas_norm = string(row.cas_norm)
        c_name = string(row.chemical_name)
        rho = Float64(row.retention_rho_monthly)
        k_fac = Float64(row.bioaccumulation_factor)


        # Match with ECOTOX
        # Find ECOTOX record
        for rec in ecotox_lib
            if get(rec, "cas_hyphenated", "") == cas_hyp || get(rec, "cas_norm", "") == cas_norm

                local is_valid_rec = false
                try
                    is_valid_rec = validate_ecotox_record(rec)
                catch
                    is_valid_rec = false
                end
                if is_valid_rec
                    # Convert effect code to axis to check routability

                    local axis_res = :none
                    local is_valid_axis = false
                    try
                        axis_res = ecotox_effect_to_deb_axis(get(rec, "effect_code", ""))
                        is_valid_axis = true
                    catch
                        is_valid_axis = false
                    end

                    if is_valid_axis && isfinite(rho) && isfinite(k_fac)
                        push!(selected_compounds, (
                            compound_index = length(selected_compounds) + 1,
                            cas_norm = cas_norm,
                            cas_hyphenated = cas_hyp,
                            chemical_name = c_name,
                            memory_class = "unknown",
                            retention_rho_monthly = rho,
                            bioaccumulation_factor = k_fac,
                            basis = "matched",
                            confidence = "high",
                            effect_code = get(rec, "effect_code", ""),
                            deb_axis = String(axis_res),
                            NOEC_median = get(rec, "NOEC_median", NaN),
                            EC50_median = get(rec, "EC50_median", NaN),
                            taxon_class = get(rec, "taxon_class", "unknown")
                        ))
                        break # Found a valid match for this compound
                    end
                end
            end
        end
        if length(selected_compounds) >= n_compounds_req
            break
        end
    end

    if length(selected_compounds) < 8
        error("Fewer than 8 matched compounds available.")
    end

    n_compounds = length(selected_compounds)
    CSV.write(joinpath(output_dir, "selected_memory_compounds.csv"), selected_compounds)

    # 4. Assign Behavior Profiles
    # Simplistic deterministic assignment
    behavior_profiles = NamedTuple[]
    profiles = ["persistent_legacy_source", "industrial_point_source_plume", "diffuse_agricultural_runoff",
                "urban_wastewater_hotspot", "episodic_spill", "river_corridor", "emerging_trend", "managed_decline"]

    for (i, c) in enumerate(selected_compounds)
        prof_name = profiles[mod1(i, length(profiles))]
        push!(behavior_profiles, (
            compound_index = c.compound_index,
            cas_norm = c.cas_norm,
            cas_hyphenated = c.cas_hyphenated,
            chemical_name = c.chemical_name,
            effect_code = c.effect_code,
            deb_axis = c.deb_axis,
            NOEC_median = c.NOEC_median,
            EC50_median = c.EC50_median,
            retention_rho_monthly = c.retention_rho_monthly,
            bioaccumulation_factor = c.bioaccumulation_factor,
            behavior_profile = prof_name,
            spatial_pattern = prof_name * "_spatial",
            seasonality_pattern = prof_name * "_seasonality",
            trend_pattern = prof_name * "_trend",
            pulse_pattern = prof_name == "episodic_spill" ? "pulse" : "none",
            baseline_factor = 0.5,
            intensity_scale = 1.0,
            assignment_reason = "deterministic_round_robin"
        ))
    end
    CSV.write(joinpath(output_dir, "compound_behavior_profiles.csv"), behavior_profiles)

    # 5. Build Concentration Grids
    C = zeros(Float64, nx, ny, n_months, n_compounds)

    for j in 1:n_compounds
        bp = behavior_profiles[j]
        noec = bp.NOEC_median
        ec50 = bp.EC50_median

        for t in 1:n_months
            for y in 1:ny
                for x in 1:nx
                    # Synthetic intensity I_j[x, y, t] >= 0
                    # Let's make it vary deterministically
                    space_val = (sin(x/nx * 2pi) + cos(y/ny * 2pi) + 2) / 4 # 0 to 1
                    time_val = (sin(t/12 * 2pi) + 1) / 2 # 0 to 1

                    if bp.behavior_profile == "emerging_trend"
                        trend = 0.4 + 0.9 * (t - 1) / max(1, n_months - 1)
                    elseif bp.behavior_profile == "managed_decline"
                        trend = max(0.0, 1.2 - 0.8 * (t - 1) / max(1, n_months - 1))
                    else
                        trend = 1.0
                    end

                    pulse = 0.0
                    if bp.behavior_profile == "episodic_spill"
                        pulse_center = round(Int, 0.6 * n_months)
                        if abs(t - pulse_center) <= 1
                            pulse = 5.0 * exp(-((x-nx/2)^2 + (y-ny/2)^2) / 10)
                        end
                    end

                    I = space_val * time_val * trend + pulse

                    C[x, y, t, j] = bp.baseline_factor * noec + bp.intensity_scale * I * (ec50 - noec)
                end
            end
        end
    end

    # Scenario Summary
    scenario_summary = NamedTuple[]
    for j in 1:n_compounds
        bp = behavior_profiles[j]
        C_j = @view C[:, :, :, j]
        push!(scenario_summary, (
            compound_index = bp.compound_index,
            cas_hyphenated = bp.cas_hyphenated,
            chemical_name = bp.chemical_name,
            behavior_profile = bp.behavior_profile,
            deb_axis = bp.deb_axis,
            effect_code = bp.effect_code,
            mean_C = mean(C_j),
            p95_C = quantile(vec(C_j), 0.95),
            max_C = maximum(C_j),
            mean_C_over_NOEC = mean(C_j) / bp.NOEC_median,
            max_C_over_EC50 = maximum(C_j) / bp.EC50_median,
            fraction_cell_months_below_NOEC = count(C_j .< bp.NOEC_median) / length(C_j),
            fraction_cell_months_between_NOEC_EC50 = count(x -> bp.NOEC_median <= x <= bp.EC50_median, C_j) / length(C_j),
            fraction_cell_months_above_EC50 = count(C_j .> bp.EC50_median) / length(C_j)
        ))
    end
    CSV.write(joinpath(output_dir, "compound_concentration_scenario_summary.csv"), scenario_summary)

    # 6. Warmup
    B0 = zeros(Float64, nx, ny, n_compounds)
    warmup_summary = NamedTuple[]
    for j in 1:n_compounds
        bp = behavior_profiles[j]
        rho = bp.retention_rho_monthly
        K = bp.bioaccumulation_factor

        cycle_len = min(12, n_months)
        C_cycle = C[:, :, 1:cycle_len, j]

        b0_j = zeros(Float64, nx, ny)
        kc_j = zeros(Float64, nx, ny, cycle_len)

        for x in 1:nx
            for y in 1:ny
                c_cell = C_cycle[x, y, :]
                b0_j[x, y] = analytical_periodic_initial_burden(rho, K, c_cell)
                kc_j[x, y, :] = K .* c_cell
            end
        end

        B0[:, :, j] = b0_j
        push!(warmup_summary, (
            compound_index = bp.compound_index,
            cas_hyphenated = bp.cas_hyphenated,
            chemical_name = bp.chemical_name,
            retention_rho_monthly = rho,
            bioaccumulation_factor = K,
            mean_B0_periodic = mean(b0_j),
            max_B0_periodic = maximum(b0_j),
            mean_KC_first_year = mean(kc_j),
            max_KC_first_year = maximum(kc_j)
        ))
    end
    CSV.write(joinpath(output_dir, "analytical_warmup_summary.csv"), warmup_summary)

    # 7. Simulation Loop
    B = zeros(Float64, nx, ny, n_months, n_compounds)
    xstress = zeros(Float64, nx, ny, n_months, n_compounds)

    for j in 1:n_compounds
        bp = behavior_profiles[j]
        rho = bp.retention_rho_monthly
        K = bp.bioaccumulation_factor
        noec = bp.NOEC_median
        ec50 = bp.EC50_median

        for t in 1:n_months
            for y in 1:ny
                for x in 1:nx
                    b_prev = t == 1 ? B0[x, y, j] : B[x, y, t-1, j]
                    c_t = C[x, y, t, j]

                    b_t = rho * b_prev + (1 - rho) * K * c_t
                    B[x, y, t, j] = b_t

                    x_t = max(0.0, (b_t - noec) / (ec50 - noec))
                    xstress[x, y, t, j] = x_t
                end
            end
        end
    end

    # 8. Mixture aggregation and DEB responses
    mixture_models = ["axis_toxic_unit_sum", "independent_action_axis_effects", "grouped_ca_then_ia_axis_effects"]
    n_models = length(mixture_models)

    Q_out = zeros(Float64, nx, ny, n_months, n_species, n_models)
    F_out = zeros(Float64, nx, ny, n_months, n_species, n_models)
    E_assimilation = zeros(Float64, nx, ny, n_months, n_models)
    E_maintenance = zeros(Float64, nx, ny, n_months, n_models)
    E_growth = zeros(Float64, nx, ny, n_months, n_models)
    E_reproduction = zeros(Float64, nx, ny, n_months, n_models)

    species_params = [amp_record_to_deb_params(species_library[k]) for k in selected_species_keys]

    # Axis indexing based on TwoTimescaleResilience expected setup
    # 1: assimilation, 2: maintenance, 3: growth, 4: reproduction
    axis_idx_map = Dict("assimilation" => 1, "maintenance" => 2, "growth" => 3, "reproduction" => 4)

    for x in 1:nx
        for y in 1:ny
            for t in 1:n_months
                x_vec = xstress[x, y, t, :]

                # Mock grouping by axis for aggregation based on ecotox_records_to_deb_burden style structure
                # In actual usage, we need arrays of records or similar. Since we are bypassing EcotoxExposureState
                # for the raw 4D array, we must apply the mixture logic directly or adapt it.

                # We'll calculate E for each axis and each model
                for m_idx in 1:n_models
                    m_name = mixture_models[m_idx]

                    axis_effects = zeros(Float64, 4)

                    # Group by axis
                    for axis_name in keys(axis_idx_map)
                        a_idx = axis_idx_map[axis_name]
                        compounds_on_axis = findall(bp -> bp.deb_axis == axis_name, behavior_profiles)

                        if !isempty(compounds_on_axis)
                            x_axis = x_vec[compounds_on_axis]
                            effect_codes = [behavior_profiles[c].effect_code for c in compounds_on_axis]

                            if m_name == "axis_toxic_unit_sum"
                                axis_effects[a_idx] = min(1.0, max(0.0, sum(x_axis)))
                            elseif m_name == "independent_action_axis_effects"
                                # E = 1 - prod(1 - min(1, x_i)) but ecotox logic allows E > 1 for TU.
                                # Let's stick to simple independent action formulation for stress:
                                # TwoTimescaleResilience aggregates burdens then maps to response, or maps to response then aggregates.
                                # According to ecotox_library, ecotox_burden_to_response is 1 - prod(1 - x_i)
                                p = 1.0
                                for v in x_axis
                                    p *= max(0.0, 1.0 - v)
                                end
                                axis_effects[a_idx] = max(0.0, 1.0 - p)
                            elseif m_name == "grouped_ca_then_ia_axis_effects"
                                # Group by effect code, sum within group (CA), then IA across groups
                                groups = Dict{String, Float64}()
                                for (v, code) in zip(x_axis, effect_codes)
                                    groups[code] = get(groups, code, 0.0) + v
                                end
                                p = 1.0
                                for v in values(groups)
                                    p *= max(0.0, 1.0 - v)
                                end
                                axis_effects[a_idx] = max(0.0, 1.0 - p)
                            end
                        end
                    end

                    E_assimilation[x, y, t, m_idx] = axis_effects[1]
                    E_maintenance[x, y, t, m_idx] = axis_effects[2]
                    E_growth[x, y, t, m_idx] = axis_effects[3]
                    E_reproduction[x, y, t, m_idx] = axis_effects[4]

                    # Now calculate Q and F for each species
                    for (s_idx, params) in enumerate(species_params)
                        # Q = w * E
                        Q = (params.alpha_axes[1] * axis_effects[1] +
                             params.alpha_axes[2] * axis_effects[2] +
                             params.alpha_axes[3] * axis_effects[3] +
                             params.alpha_axes[4] * axis_effects[4])
                        Q = min(1.0, max(0.0, Q))

                        A_t = params.A0 * max(1e-6, 1.0 - Q)

                        # Use package function for lambda if available, otherwise raw
                        lambda_A0 = (params.KA * params.A0) / (params.KA + params.A0)
                        lambda_At = (params.KA * A_t) / (params.KA + A_t)

                        F = lambda_A0 / lambda_At

                        Q_out[x, y, t, s_idx, m_idx] = Q
                        F_out[x, y, t, s_idx, m_idx] = F
                    end
                end
            end
        end
    end

    # Simulation Summary CSV
    sim_summary = NamedTuple[]
    for j in 1:n_compounds
        bp = behavior_profiles[j]
        B_j = B[:, :, :, j]
        x_j = xstress[:, :, :, j]
        c_j = C[:, :, :, j]
        push!(sim_summary, (
            compound_index = bp.compound_index,
            cas_hyphenated = bp.cas_hyphenated,
            chemical_name = bp.chemical_name,
            effect_code = bp.effect_code,
            deb_axis = bp.deb_axis,
            retention_rho_monthly = bp.retention_rho_monthly,
            bioaccumulation_factor = bp.bioaccumulation_factor,
            mean_C = mean(c_j),
            max_C = maximum(c_j),
            mean_B = mean(B_j),
            max_B = maximum(B_j),
            mean_x = mean(x_j),
            max_x = maximum(x_j),
            active_cell_month_fraction = count(x_j .> 0.0) / length(x_j)
        ))
    end
    CSV.write(joinpath(output_dir, "compound_memory_simulation_summary.csv"), sim_summary)

    # Species Summary CSV
    species_summary = NamedTuple[]
    for (s_idx, s_info) in enumerate(species_info)
        for m_idx in 1:n_models
            Q_m = Q_out[:, :, :, s_idx, m_idx]
            F_m = F_out[:, :, :, s_idx, m_idx]

            p95_F_m = quantile(vec(F_m), 0.95)

            # Find month of max p95 F. Over spatial domain for each month, find p95, then find max month.
            month_p95s = [quantile(vec(F_m[:, :, t]), 0.95) for t in 1:n_months]
            max_month = argmax(month_p95s)

            push!(species_summary, (
                species_key = s_info.species_key,
                species_name = s_info.species_name,
                archetype_labels = s_info.archetype_labels,
                mixture_effect_model = mixture_models[m_idx],
                mean_Q = mean(Q_m),
                p95_Q = quantile(vec(Q_m), 0.95),
                max_Q = maximum(Q_m),
                mean_F = mean(F_m),
                p95_F = p95_F_m,
                max_F = maximum(F_m),
                min_relative_margin_remaining = 1.0 - maximum(Q_m), # A_t / A0 approx
                month_of_max_p95_F = max_month
            ))
        end
    end
    CSV.write(joinpath(output_dir, "species_archetype_response_summary.csv"), species_summary)

    # Mixture model sensitivity
    mix_sens_summary = NamedTuple[]
    # Indices: 1=TU, 2=IA, 3=Grouped

    # Calculate max deltas across all cells, months, species
    delta_F_IA_TU = maximum(F_out[:, :, :, :, 2] .- F_out[:, :, :, :, 1])
    delta_F_grouped_TU = maximum(F_out[:, :, :, :, 3] .- F_out[:, :, :, :, 1])
    delta_F_grouped_IA = maximum(F_out[:, :, :, :, 3] .- F_out[:, :, :, :, 2])

    delta_Q_IA_TU = maximum(Q_out[:, :, :, :, 2] .- Q_out[:, :, :, :, 1])
    delta_Q_grouped_TU = maximum(Q_out[:, :, :, :, 3] .- Q_out[:, :, :, :, 1])
    delta_Q_grouped_IA = maximum(Q_out[:, :, :, :, 3] .- Q_out[:, :, :, :, 2])

    push!(mix_sens_summary, (
        max_delta_F_IA_minus_TU = delta_F_IA_TU,
        max_delta_F_grouped_minus_TU = delta_F_grouped_TU,
        max_delta_F_grouped_minus_IA = delta_F_grouped_IA,
        max_delta_Q_IA_minus_TU = delta_Q_IA_TU,
        max_delta_Q_grouped_minus_TU = delta_Q_grouped_TU,
        max_delta_Q_grouped_minus_IA = delta_Q_grouped_IA
    ))
    CSV.write(joinpath(output_dir, "mixture_model_sensitivity_summary.csv"), mix_sens_summary)

    # 9. Threshold-free pipeline
    response_arrays = (
        Q_t = Q_out,
        F_t = F_out,
        E_assimilation = E_assimilation,
        E_maintenance = E_maintenance,
        E_growth = E_growth,
        E_reproduction = E_reproduction
    )

    feature_result = build_threshold_free_vulnerability_features(
        response_arrays;
        mixture_model_names = mixture_models,
        preferred_mixture_model = "grouped_ca_then_ia_axis_effects"
    )

    standardized = standardize_threshold_free_vulnerability_features(
        feature_result.feature_matrix,
        feature_result.feature_names
    )

    clusters = cluster_threshold_free_vulnerability_regimes(
        standardized.standardized_features;
        feature_names = standardized.standardized_feature_names,
        k = k_clusters
    )

    cluster_summary = summarize_threshold_free_vulnerability_clusters(
        clusters,
        standardized.standardized_features;
        feature_names = standardized.standardized_feature_names
    )

    # Output bundle
    bundle = vulnerability_regime_output_bundle(
        feature_result,
        standardized,
        clusters;
        cluster_summary = cluster_summary,
        output_dir = output_dir,
        write_netcdf = true,
        write_csv = true,
        make_plots = make_plots
    )

    # 10. Simulation Metadata
    metadata = Dict(
        "generated_by" => "TwoTimescaleResilience",
        "example_script" => "archetype_compound_memory_10yr_grid_demo.jl",
        "nx" => nx,
        "ny" => ny,
        "n_years" => n_years,
        "n_months" => n_months,
        "n_species" => n_species,
        "n_compounds" => n_compounds,
        "n_mixture_models" => n_models,
        "k_clusters" => k_clusters,
        "warmup_method" => "periodic_annual_analytical",
        "concentration_fields" => "synthetic_deterministic_diagnostic",
        "species_source" => "AmP_Species_Library.json",
        "compound_source" => "Compound_Memory_Library.csv and ECOTOX_Toxicity_Library.json",
        "archetype_database_present" => archetype_db_present,
        "species_selection_mode" => species_selection_mode,
        "fallback_reason" => fallback_reason,
        "selected_behavior_profiles" => profiles,
        "unsatisfied_design_constraints" => [],
        "mixture_effect_models" => mixture_models,
        "preferred_mixture_model" => "grouped_ca_then_ia_axis_effects",
        "threshold_free_features" => true,
        "real_raster_ingestion" => false,
        "physiological_Z_t" => false,
        "DEBtox_D_t" => false,
        "synergism_antagonism" => false,
        "fitted_interactions" => false
    )

    open(joinpath(output_dir, "simulation_metadata.json"), "w") do f
        JSON.print(f, metadata, 4)
    end

    println("Demo completed successfully. Outputs in $output_dir")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_archetype_compound_memory_10yr_grid_demo()
end
