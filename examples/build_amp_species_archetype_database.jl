# ===========================================================================
# Build AmP Species Archetype Database
#
# This script classifies AmP species into reusable response-capacity
# archetypes (e.g., amplification potential, vulnerability thresholds,
# dominant sensitive axes) to support downstream simulations and examples.
# ===========================================================================

using TwoTimescaleResilience
using JSON
using CSV
using DataFrames
using Dates
using Statistics

# ===========================================================================
# Configuration
# ===========================================================================
DATA_DIR = joinpath(@__DIR__, "..", "data")
AMP_LIBRARY_PATH = joinpath(DATA_DIR, "AmP_Species_Library.json")
DIAGNOSTIC_CSV_PATH = joinpath(@__DIR__, "..", "output", "amp_species_response_capacity_diagnostics", "amp_species_response_capacity.csv")
OUTPUT_DIR_SUMMARY = joinpath(@__DIR__, "..", "output", "amp_species_archetype_database")
OUTPUT_CSV_PATH = joinpath(DATA_DIR, "AmP_Species_Archetypes.csv")
OUTPUT_JSON_PATH = joinpath(DATA_DIR, "AmP_Species_Archetypes.json")
OUTPUT_SUMMARY_PATH = joinpath(OUTPUT_DIR_SUMMARY, "amp_species_archetype_database_build_summary.csv")

mkpath(DATA_DIR)
mkpath(OUTPUT_DIR_SUMMARY)

# ===========================================================================
# Helper Functions for Metric Calculation (Fallback)
# ===========================================================================
function safe_divide(num, den)
    if den == 0.0 || isnan(den)
        return Inf
    end
    return num / den
end

function get_x_collapse(A0, alpha)
    if alpha <= 0.0 || isnan(alpha)
        return Inf
    end
    return A0 / alpha
end

function get_Fmax(lambda0, lambda_min)
    if lambda_min <= 0.0
        return Inf
    end
    return lambda0 / lambda_min
end

# ===========================================================================
# Main Routine
# ===========================================================================
function main()
    println("Loading AmP Species Library...")
    if !isfile(AMP_LIBRARY_PATH)
        error("AmP library not found at $AMP_LIBRARY_PATH")
    end
    amp_lib = load_amp_species_library(AMP_LIBRARY_PATH)

    # 1. Load or compute metrics
    metrics_df = nothing
    source_diagnostic = ""

    if isfile(DIAGNOSTIC_CSV_PATH)
        println("Loading existing diagnostics from CSV...")
        metrics_df = CSV.read(DIAGNOSTIC_CSV_PATH, DataFrame)
        source_diagnostic = DIAGNOSTIC_CSV_PATH
    else
        println("Diagnostics CSV not found. Computing metrics dynamically...")
        source_diagnostic = "Computed dynamically"
        rows = []
        for (key, rec) in amp_lib
            try
                params = amp_record_to_deb_params(rec)
                A0 = params.A0
                lambda_min = params.lambda_min
                lambda_max = params.lambda_max
                lambda0 = TwoTimescaleResilience.restoring_force_from_margin(A0, params)

                Fmax = get_Fmax(lambda0, lambda_min)

                x_c_a = get_x_collapse(A0, params.alpha_axes[1])
                x_c_m = get_x_collapse(A0, params.alpha_axes[2])
                x_c_g = get_x_collapse(A0, params.alpha_axes[3])
                x_c_r = get_x_collapse(A0, params.alpha_axes[4])

                c_vals = (assimilation=x_c_a, maintenance=x_c_m, growth=x_c_g, reproduction=x_c_r)
                min_c = minimum(c_vals)
                most_sens = if min_c == x_c_a
                    :assimilation
                elseif min_c == x_c_m
                    :maintenance
                elseif min_c == x_c_g
                    :growth
                else
                    :reproduction
                end

                push!(rows, (
                    species_key = key,
                    species_name = replace(key, "_" => " "),
                    A0 = A0,
                    lambda_min = lambda_min,
                    lambda_max = lambda_max,
                    lambda0 = lambda0,
                    Fmax = Fmax,
                    lambda0_relative_position = (lambda0 - lambda_min) / (lambda_max - lambda_min),
                    alpha_assimilation = params.alpha_axes[1],
                    alpha_maintenance = params.alpha_axes[2],
                    alpha_growth = params.alpha_axes[3],
                    alpha_reproduction = params.alpha_axes[4],
                    x_collapse_assimilation = x_c_a,
                    x_collapse_maintenance = x_c_m,
                    x_collapse_growth = x_c_g,
                    x_collapse_reproduction = x_c_r,
                    min_x_collapse = min_c,
                    most_sensitive_axis = String(most_sens)
                ))
            catch
                continue # Skip invalid
            end
        end
        metrics_df = DataFrame(rows)
    end

    # Filter to valid species
    valid_species_keys = Set(keys(amp_lib))
    metrics_df = filter(row -> row.species_key in valid_species_keys, metrics_df)

    # Sort deterministically
    sort!(metrics_df, :species_key)

    n_species = nrow(metrics_df)
    println("Processing $n_species species...")

    # 2. Compute Quantiles
    finite_fmax = filter(x -> isfinite(x) && x > 0, metrics_df.Fmax)
    fmax_p20 = isempty(finite_fmax) ? NaN : quantile(finite_fmax, 0.20)
    fmax_p40 = isempty(finite_fmax) ? NaN : quantile(finite_fmax, 0.40)
    fmax_p60 = isempty(finite_fmax) ? NaN : quantile(finite_fmax, 0.60)
    fmax_p80 = isempty(finite_fmax) ? NaN : quantile(finite_fmax, 0.80)
    fmax_p95 = isempty(finite_fmax) ? NaN : quantile(finite_fmax, 0.95)

    finite_xc = filter(x -> isfinite(x) && x > 0, metrics_df.min_x_collapse)
    xc_p20 = isempty(finite_xc) ? NaN : quantile(finite_xc, 0.20)
    xc_p40 = isempty(finite_xc) ? NaN : quantile(finite_xc, 0.40)
    xc_p60 = isempty(finite_xc) ? NaN : quantile(finite_xc, 0.60)
    xc_p80 = isempty(finite_xc) ? NaN : quantile(finite_xc, 0.80)
    xc_p95 = isempty(finite_xc) ? NaN : quantile(finite_xc, 0.95)

    # Calculate percentiles relative to distribution for columns
    function get_quantile(val, arr)
        if !isfinite(val) || val <= 0 || isempty(arr)
            return NaN
        end
        return count(x -> x <= val, arr) / length(arr)
    end

    # 3. Assign Archetypes
    archetypes = DataFrame()
    archetypes.species_key = metrics_df.species_key
    archetypes.species_name = metrics_df.species_name
    archetypes.A0 = metrics_df.A0
    archetypes.lambda_min = metrics_df.lambda_min
    archetypes.lambda_max = metrics_df.lambda_max
    archetypes.lambda0 = metrics_df.lambda0
    archetypes.Fmax = metrics_df.Fmax
    archetypes.lambda0_relative_position = metrics_df.lambda0_relative_position
    archetypes.alpha_assimilation = metrics_df.alpha_assimilation
    archetypes.alpha_maintenance = metrics_df.alpha_maintenance
    archetypes.alpha_growth = metrics_df.alpha_growth
    archetypes.alpha_reproduction = metrics_df.alpha_reproduction
    archetypes.x_collapse_assimilation = metrics_df.x_collapse_assimilation
    archetypes.x_collapse_maintenance = metrics_df.x_collapse_maintenance
    archetypes.x_collapse_growth = metrics_df.x_collapse_growth
    archetypes.x_collapse_reproduction = metrics_df.x_collapse_reproduction
    archetypes.min_x_collapse = metrics_df.min_x_collapse
    archetypes.most_sensitive_axis = metrics_df.most_sensitive_axis

    # Derived Metrics
    alpha_sum = archetypes.alpha_assimilation .+ archetypes.alpha_maintenance .+ archetypes.alpha_growth .+ archetypes.alpha_reproduction

    archetypes.alpha_fraction_assimilation = ifelse.(alpha_sum .> 0, archetypes.alpha_assimilation ./ alpha_sum, NaN)
    archetypes.alpha_fraction_maintenance = ifelse.(alpha_sum .> 0, archetypes.alpha_maintenance ./ alpha_sum, NaN)
    archetypes.alpha_fraction_growth = ifelse.(alpha_sum .> 0, archetypes.alpha_growth ./ alpha_sum, NaN)
    archetypes.alpha_fraction_reproduction = ifelse.(alpha_sum .> 0, archetypes.alpha_reproduction ./ alpha_sum, NaN)

    archetypes.Fmax_quantile = [get_quantile(v, finite_fmax) for v in archetypes.Fmax]
    archetypes.min_x_collapse_quantile = [get_quantile(v, finite_xc) for v in archetypes.min_x_collapse]

    # Initialize Indicators
    indicators = [
        :is_low_amplification_potential, :is_intermediate_amplification_potential,
        :is_high_amplification_potential, :is_extreme_amplification_potential,
        :is_high_sensitivity_low_collapse_threshold, :is_intermediate_collapse_threshold,
        :is_robust_high_collapse_threshold, :is_extreme_robust_high_collapse_threshold,
        :is_assimilation_sensitive, :is_maintenance_sensitive, :is_growth_sensitive,
        :is_reproduction_sensitive, :is_balanced_axis_sensitivity,
        :is_vulnerable_high_amplification_low_threshold, :is_robust_low_amplification
    ]
    for ind in indicators
        archetypes[!, ind] .= false
    end

    labels_list = Vector{String}[]

    for i in 1:n_species
        lbls = String[]
        row = archetypes[i, :]

        # Fmax
        has_high_amp = false
        if isfinite(row.Fmax) && row.Fmax > 0
            if row.Fmax <= fmax_p20
                push!(lbls, "low_amplification_potential")
                archetypes[i, :is_low_amplification_potential] = true
            end
            if fmax_p40 <= row.Fmax <= fmax_p60
                push!(lbls, "intermediate_amplification_potential")
                archetypes[i, :is_intermediate_amplification_potential] = true
            end
            if row.Fmax >= fmax_p80
                push!(lbls, "high_amplification_potential")
                archetypes[i, :is_high_amplification_potential] = true
                has_high_amp = true
            end
            if row.Fmax >= fmax_p95
                push!(lbls, "extreme_amplification_potential")
                archetypes[i, :is_extreme_amplification_potential] = true
                has_high_amp = true
            end
        end

        # min_x_collapse
        has_low_thresh = false
        has_robust_thresh = false
        if isfinite(row.min_x_collapse) && row.min_x_collapse > 0
            if row.min_x_collapse <= xc_p20
                push!(lbls, "high_sensitivity_low_collapse_threshold")
                archetypes[i, :is_high_sensitivity_low_collapse_threshold] = true
                has_low_thresh = true
            end
            if xc_p40 <= row.min_x_collapse <= xc_p60
                push!(lbls, "intermediate_collapse_threshold")
                archetypes[i, :is_intermediate_collapse_threshold] = true
            end
            if row.min_x_collapse >= xc_p80
                push!(lbls, "robust_high_collapse_threshold")
                archetypes[i, :is_robust_high_collapse_threshold] = true
                has_robust_thresh = true
            end
            if row.min_x_collapse >= xc_p95
                push!(lbls, "extreme_robust_high_collapse_threshold")
                archetypes[i, :is_extreme_robust_high_collapse_threshold] = true
                has_robust_thresh = true
            end
        elseif isinf(row.min_x_collapse) && row.min_x_collapse > 0
            # Inf means extremely robust
            push!(lbls, "robust_high_collapse_threshold")
            archetypes[i, :is_robust_high_collapse_threshold] = true
            push!(lbls, "extreme_robust_high_collapse_threshold")
            archetypes[i, :is_extreme_robust_high_collapse_threshold] = true
            has_robust_thresh = true
        end

        # Sensitivity axis
        ms = String(row.most_sensitive_axis)
        if ms == "assimilation"
            push!(lbls, "assimilation_sensitive")
            archetypes[i, :is_assimilation_sensitive] = true
        elseif ms == "maintenance"
            push!(lbls, "maintenance_sensitive")
            archetypes[i, :is_maintenance_sensitive] = true
        elseif ms == "growth"
            push!(lbls, "growth_sensitive")
            archetypes[i, :is_growth_sensitive] = true
        elseif ms == "reproduction"
            push!(lbls, "reproduction_sensitive")
            archetypes[i, :is_reproduction_sensitive] = true
        end

        # Fractions
        if isfinite(row.alpha_fraction_assimilation)
            max_frac = max(row.alpha_fraction_assimilation, row.alpha_fraction_maintenance, row.alpha_fraction_growth, row.alpha_fraction_reproduction)
            if row.alpha_fraction_assimilation >= 0.60 push!(lbls, "strongly_assimilation_weighted") end
            if row.alpha_fraction_maintenance >= 0.60 push!(lbls, "strongly_maintenance_weighted") end
            if row.alpha_fraction_growth >= 0.60 push!(lbls, "strongly_growth_weighted") end
            if row.alpha_fraction_reproduction >= 0.60 push!(lbls, "strongly_reproduction_weighted") end
            if max_frac <= 0.40
                push!(lbls, "balanced_axis_sensitivity")
                archetypes[i, :is_balanced_axis_sensitivity] = true
            end
        end

        # Combined
        if has_high_amp && ms == "maintenance" push!(lbls, "high_Fmax_maintenance_sensitive") end
        if has_high_amp && ms == "growth" push!(lbls, "high_Fmax_growth_sensitive") end
        if has_high_amp && ms == "reproduction" push!(lbls, "high_Fmax_reproduction_sensitive") end

        if has_low_thresh && ms == "maintenance" push!(lbls, "low_threshold_maintenance_sensitive") end
        if has_low_thresh && ms == "growth" push!(lbls, "low_threshold_growth_sensitive") end
        if has_low_thresh && ms == "reproduction" push!(lbls, "low_threshold_reproduction_sensitive") end

        if has_robust_thresh && archetypes[i, :is_low_amplification_potential]
            push!(lbls, "robust_low_amplification")
            archetypes[i, :is_robust_low_amplification] = true
        end
        if has_high_amp && has_low_thresh
            push!(lbls, "vulnerable_high_amplification_low_threshold")
            archetypes[i, :is_vulnerable_high_amplification_low_threshold] = true
        end

        push!(labels_list, unique(lbls))
    end

    archetypes.archetype_labels = [join(lbls, ";") for lbls in labels_list]

    # Provenance
    gen_time = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-dd\\THH:MM:SS\\Z")
    archetypes.source_amp_library .= "data/AmP_Species_Library.json"
    archetypes.source_response_capacity_diagnostic .= source_diagnostic
    archetypes.archetype_algorithm_version .= "1.0"
    archetypes.generated_at_utc .= gen_time

    # Order columns
    cols = ["species_key", "species_name", "archetype_labels",
            "A0", "lambda_min", "lambda_max", "lambda0", "Fmax", "lambda0_relative_position",
            "alpha_assimilation", "alpha_maintenance", "alpha_growth", "alpha_reproduction",
            "alpha_fraction_assimilation", "alpha_fraction_maintenance", "alpha_fraction_growth", "alpha_fraction_reproduction",
            "x_collapse_assimilation", "x_collapse_maintenance", "x_collapse_growth", "x_collapse_reproduction",
            "min_x_collapse", "most_sensitive_axis", "Fmax_quantile", "min_x_collapse_quantile"]
    for ind in indicators
        push!(cols, String(ind))
    end
    push!(cols, "source_amp_library", "source_response_capacity_diagnostic", "archetype_algorithm_version", "generated_at_utc")

    archetypes = archetypes[:, cols]

    # 4. Write CSV
    CSV.write(OUTPUT_CSV_PATH, archetypes)
    println("Wrote CSV to $OUTPUT_CSV_PATH")

    # 5. Write JSON
    json_species = []
    for i in 1:nrow(archetypes)
        row = archetypes[i, :]
        push!(json_species, Dict(
            "species_key" => row.species_key,
            "species_name" => row.species_name,
            "archetype_labels" => labels_list[i],
            "metrics" => Dict(
                "A0" => row.A0,
                "lambda_min" => row.lambda_min,
                "lambda_max" => row.lambda_max,
                "lambda0" => row.lambda0,
                "Fmax" => row.Fmax,
                "min_x_collapse" => row.min_x_collapse,
                "most_sensitive_axis" => row.most_sensitive_axis
            )
        ))
    end

    # Collect all unique archetypes
    all_archetypes = unique(vcat(labels_list...))

    json_data = Dict(
        "metadata" => Dict(
            "generated_by" => "examples/build_amp_species_archetype_database.jl",
            "archetype_algorithm_version" => "1.0",
            "source_amp_library" => "data/AmP_Species_Library.json",
            "source_response_capacity_diagnostic" => source_diagnostic,
            "generated_at_utc" => gen_time,
            "n_species" => n_species,
            "n_archetypes" => length(all_archetypes)
        ),
        "archetype_definitions" => Dict(
            "low_amplification_potential" => "Fmax <= 20th percentile",
            "high_amplification_potential" => "Fmax >= 80th percentile",
            "high_sensitivity_low_collapse_threshold" => "min_x_collapse <= 20th percentile",
            "robust_high_collapse_threshold" => "min_x_collapse >= 80th percentile",
            "assimilation_sensitive" => "assimilation is most sensitive axis",
            "maintenance_sensitive" => "maintenance is most sensitive axis",
            "growth_sensitive" => "growth is most sensitive axis",
            "reproduction_sensitive" => "reproduction is most sensitive axis",
            # Additional keys can go here if needed
        ),
        "species" => json_species
    )

    open(OUTPUT_JSON_PATH, "w") do f
        JSON.print(f, json_data, 2)
    end
    println("Wrote JSON to $OUTPUT_JSON_PATH")

    # 6. Write Build Summary
    summary_rows = []
    for arch in all_archetypes
        arch_indices = findall(labels -> arch in labels, labels_list)
        arch_df = archetypes[arch_indices, :]

        # Example species
        example_species = arch_df[1, :species_key]
        example_name = arch_df[1, :species_name]

        # Fmax stats
        f_vals = filter(x -> isfinite(x) && x > 0, arch_df.Fmax)
        f_min = isempty(f_vals) ? NaN : minimum(f_vals)
        f_med = isempty(f_vals) ? NaN : median(f_vals)
        f_max = isempty(f_vals) ? NaN : maximum(f_vals)

        # xc stats
        xc_vals = filter(x -> isfinite(x) && x > 0, arch_df.min_x_collapse)
        xc_min = isempty(xc_vals) ? NaN : minimum(xc_vals)
        xc_med = isempty(xc_vals) ? NaN : median(xc_vals)
        xc_max = isempty(xc_vals) ? NaN : maximum(xc_vals)

        push!(summary_rows, (
            archetype_label = arch,
            n_species = length(arch_indices),
            example_species_key = example_species,
            example_species_name = example_name,
            min_Fmax = f_min,
            median_Fmax = f_med,
            max_Fmax = f_max,
            min_min_x_collapse = xc_min,
            median_min_x_collapse = xc_med,
            max_min_x_collapse = xc_max
        ))
    end
    summary_df = DataFrame(summary_rows)
    sort!(summary_df, :n_species, rev=true)
    CSV.write(OUTPUT_SUMMARY_PATH, summary_df)
    println("Wrote Build Summary to $OUTPUT_SUMMARY_PATH")

    # Console Summary
    n_assigned = count(labels -> !isempty(labels), labels_list)
    println("\n========================================")
    println("BUILD SUMMARY")
    println("========================================")
    println("Species processed: $n_species")
    println("Species assigned at least one archetype: $n_assigned")
    println("Number of species with missing/non-finite Fmax: ", count(x -> !isfinite(x) || x <= 0, archetypes.Fmax))
    println("\nCounts by archetype:")
    for row in eachrow(summary_df)
        println("  $(row.archetype_label): $(row.n_species)")
    end
    println("========================================")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
