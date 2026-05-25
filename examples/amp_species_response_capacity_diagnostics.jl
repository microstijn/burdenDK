# ===========================================================================
# AmP Species Response Capacity Diagnostics
#
# This script computes species-level response capacity directly from
# AmP-derived DEBAxisParams. It analytically answers whether large
# amplification is possible under the current parameterization.
# ===========================================================================

using Pkg
try
    Pkg.activate(joinpath(@__DIR__, ".."))
catch err
    @warn "Could not activate project" exception=(err, catch_backtrace())
end

using TwoTimescaleResilience
using Printf
using Statistics
using CairoMakie

# Optional DataFrames/CSV
csv_available = false
try
    eval(:(using DataFrames))
    eval(:(using CSV))
    global csv_available = true
catch
    global csv_available = false
end

# Ensure output directory exists
const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "amp_species_response_capacity_diagnostics")
mkpath(OUTPUT_DIR)

# ===========================================================================
# Helper Functions
# ===========================================================================

human_name(species_key::String) = replace(species_key, "_" => " ")

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
    return safe_divide(A0, alpha)
end

function get_Fmax(lambda0, lambda_min)
    if lambda_min <= 0.0 || isnan(lambda_min)
        @warn "Invalid lambda_min ($lambda_min). Setting Fmax to NaN."
        return NaN
    end
    return lambda0 / lambda_min
end

# ===========================================================================
# Data Loading and Processing Loop
# ===========================================================================

amp_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
amp = load_amp_species_library(amp_path)
species_keys = sort(collect(keys(amp)))

println("\n============================================================")
println("AmP Species Response Capacity Diagnostics")
println("============================================================")
println("Total species analyzed: ", length(species_keys))
println("\nAnalytical explanation:")
println("  Fmax = lambda(A0) / lambda_min")
println("  x_collapse_axis = A0 / alpha_axis")

results = []

for species_key in species_keys
    params = amp_species_deb_params(amp, species_key)

    A0 = params.A0
    alpha_axes = params.alpha_axes
    lambda_min = params.lambda_min
    lambda_max = params.lambda_max
    KA = params.KA

    lambda0 = restoring_force_from_margin(A0, params)
    Fmax = get_Fmax(lambda0, lambda_min)

    lambda_dynamic_range = lambda_max - lambda_min
    lambda0_relative_position = lambda_dynamic_range > 0.0 ? (lambda0 - lambda_min) / lambda_dynamic_range : NaN

    A0_over_KA = KA > 0.0 ? A0 / KA : NaN

    x_collapse_assimilation = get_x_collapse(A0, alpha_axes[1])
    x_collapse_maintenance = get_x_collapse(A0, alpha_axes[2])
    x_collapse_growth = get_x_collapse(A0, alpha_axes[3])
    x_collapse_reproduction = get_x_collapse(A0, alpha_axes[4])

    collapses = [x_collapse_assimilation, x_collapse_maintenance, x_collapse_growth, x_collapse_reproduction]
    axis_names = (:assimilation, :maintenance, :growth, :reproduction)

    finite_collapses = filter(x -> !isinf(x) && !isnan(x), collapses)
    min_x_collapse = isempty(finite_collapses) ? Inf : minimum(finite_collapses)

    most_sensitive_axis = :none
    if !isinf(min_x_collapse)
        idx = findfirst(x -> x == min_x_collapse, collapses)
        if !isnothing(idx)
            most_sensitive_axis = axis_names[idx]
        end
    end

    push!(results, (
        species_key = species_key,
        species_name = human_name(species_key),
        A0 = A0,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        KA = KA,
        lambda0 = lambda0,
        Fmax = Fmax,
        lambda0_relative_position = lambda0_relative_position,
        A0_over_KA = A0_over_KA,
        alpha_assimilation = alpha_axes[1],
        alpha_maintenance = alpha_axes[2],
        alpha_growth = alpha_axes[3],
        alpha_reproduction = alpha_axes[4],
        x_collapse_assimilation = x_collapse_assimilation,
        x_collapse_maintenance = x_collapse_maintenance,
        x_collapse_growth = x_collapse_growth,
        x_collapse_reproduction = x_collapse_reproduction,
        min_x_collapse = min_x_collapse,
        most_sensitive_axis = most_sensitive_axis
    ))
end

# ===========================================================================
# Printed Summaries and Rankings
# ===========================================================================

println("\n------------------------------------------------------------")
println("Preferred Species Diagnostics")
println("------------------------------------------------------------")

preferred = ["Abatus_cordatus", "Podarcis_muralis", "Thalia_democratica"]

for sp_key in preferred
    idx = findfirst(r -> r.species_key == sp_key, results)
    if !isnothing(idx)
        r = results[idx]
        println("\nSpecies: ", r.species_name)
        @printf("  A0:                        %.4f\n", r.A0)
        @printf("  alpha axes (a, m, g, r):   (%.4f, %.4f, %.4f, %.4f)\n", r.alpha_assimilation, r.alpha_maintenance, r.alpha_growth, r.alpha_reproduction)
        @printf("  lambda_min:                %.4f\n", r.lambda_min)
        @printf("  lambda_max:                %.4f\n", r.lambda_max)
        @printf("  KA:                        %.4f\n", r.KA)
        @printf("  lambda0:                   %.4f\n", r.lambda0)
        @printf("  Fmax:                      %.4f\n", r.Fmax)
        @printf("  A0 / KA:                   %.4f\n", r.A0_over_KA)
        @printf("  lambda0 relative position: %.4f\n", r.lambda0_relative_position)
        @printf("  x_collapse (a, m, g, r):   (%.4f, %.4f, %.4f, %.4f)\n", r.x_collapse_assimilation, r.x_collapse_maintenance, r.x_collapse_growth, r.x_collapse_reproduction)
        println("  most sensitive axis:       ", r.most_sensitive_axis)
    end
end

println("\n------------------------------------------------------------")
println("Ranking 1: largest Fmax (Top 15)")
println("------------------------------------------------------------")
sorted_fmax_desc = sort(results, by = x -> isnan(x.Fmax) ? -Inf : x.Fmax, rev = true)
@printf("%-4s %-30s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n", "Rank", "Species", "Fmax", "lambda0", "lambda_min", "lambda_max", "A0", "KA", "A0/KA")
for i in 1:min(15, length(sorted_fmax_desc))
    r = sorted_fmax_desc[i]
    @printf("%-4d %-30s %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n", i, r.species_name, r.Fmax, r.lambda0, r.lambda_min, r.lambda_max, r.A0, r.KA, r.A0_over_KA)
end

println("\n------------------------------------------------------------")
println("Ranking 2: smallest Fmax (Bottom 15)")
println("------------------------------------------------------------")
sorted_fmax_asc = sort(results, by = x -> isnan(x.Fmax) ? Inf : x.Fmax)
@printf("%-4s %-30s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n", "Rank", "Species", "Fmax", "lambda0", "lambda_min", "lambda_max", "A0", "KA", "A0/KA")
for i in 1:min(15, length(sorted_fmax_asc))
    r = sorted_fmax_asc[i]
    @printf("%-4d %-30s %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n", i, r.species_name, r.Fmax, r.lambda0, r.lambda_min, r.lambda_max, r.A0, r.KA, r.A0_over_KA)
end

println("\n------------------------------------------------------------")
println("Ranking 3: lowest collapse burden (Top 15)")
println("------------------------------------------------------------")
sorted_collapse_asc = sort(results, by = x -> x.min_x_collapse)
@printf("%-4s %-30s %-15s %-20s %-10s %-10s %-10s\n", "Rank", "Species", "min_x_collapse", "most_sensitive_axis", "A0", "alpha", "Fmax")
for i in 1:min(15, length(sorted_collapse_asc))
    r = sorted_collapse_asc[i]
    alpha_val = if r.most_sensitive_axis == :assimilation
        r.alpha_assimilation
    elseif r.most_sensitive_axis == :maintenance
        r.alpha_maintenance
    elseif r.most_sensitive_axis == :growth
        r.alpha_growth
    elseif r.most_sensitive_axis == :reproduction
        r.alpha_reproduction
    else
        0.0
    end
    @printf("%-4d %-30s %-15.4f %-20s %-10.4f %-10.4f %-10.4f\n", i, r.species_name, r.min_x_collapse, String(r.most_sensitive_axis), r.A0, alpha_val, r.Fmax)
end

println("\n------------------------------------------------------------")
println("Ranking 4: highest collapse burden (Top 15)")
println("------------------------------------------------------------")
finite_collapse = filter(x -> !isinf(x.min_x_collapse), results)
sorted_collapse_desc = sort(finite_collapse, by = x -> x.min_x_collapse, rev = true)
@printf("%-4s %-30s %-15s %-20s %-10s %-10s %-10s\n", "Rank", "Species", "min_x_collapse", "most_sensitive_axis", "A0", "alpha", "Fmax")
for i in 1:min(15, length(sorted_collapse_desc))
    r = sorted_collapse_desc[i]
    alpha_val = if r.most_sensitive_axis == :assimilation
        r.alpha_assimilation
    elseif r.most_sensitive_axis == :maintenance
        r.alpha_maintenance
    elseif r.most_sensitive_axis == :growth
        r.alpha_growth
    elseif r.most_sensitive_axis == :reproduction
        r.alpha_reproduction
    else
        0.0
    end
    @printf("%-4d %-30s %-15.4f %-20s %-10.4f %-10.4f %-10.4f\n", i, r.species_name, r.min_x_collapse, String(r.most_sensitive_axis), r.A0, alpha_val, r.Fmax)
end

println("\n============================================================")
println("ABATUS ANALYTICAL EXPLANATION")
println("============================================================")
idx_abatus = findfirst(r -> r.species_key == "Abatus_cordatus", results)
if !isnothing(idx_abatus)
    params = amp_species_deb_params(amp, "Abatus_cordatus")
    x_example = 47.135
    depletion_M = params.alpha_axes[2] * x_example
    fraction_A0_depleted = depletion_M / params.A0
    x_collapse_M = get_x_collapse(params.A0, params.alpha_axes[2])
    fraction_to_collapse = x_example / x_collapse_M
    r_ab = results[idx_abatus]

    println("For Abatus cordatus:")
    @printf("- maintenance alpha = %.4f\n", params.alpha_axes[2])
    @printf("- A0 = %.4f\n", params.A0)
    @printf("- peak x example = %.3f\n", x_example)
    @printf("- A depletion from this maintenance burden = %.4f\n", depletion_M)
    @printf("- fraction of A0 depleted = %.4f\n", fraction_A0_depleted)
    @printf("- maintenance x needed to collapse A to zero = %.4f\n", x_collapse_M)
    @printf("- fraction of collapse burden reached = %.4f\n", fraction_to_collapse)
    @printf("- Fmax = %.4f\n\n", r_ab.Fmax)

    println("This explains why the monthly memory demo can show large internal burden and large active stress while F remains near 1. The active stress is still far below the stress required to collapse adaptive margin, and the maximum amplification for this parameterization is limited by lambda0 / lambda_min.")
end


# ===========================================================================
# CSV Output (Conditional)
# ===========================================================================

println("\n============================================================")
println("OPTIONAL CSV OUTPUT")
println("============================================================")
if csv_available
    try
        df = DataFrame(results)
        csv_path = joinpath(OUTPUT_DIR, "amp_species_response_capacity.csv")
        CSV.write(csv_path, df)
        println("Successfully wrote results to:")
        println("  $csv_path")
    catch e
        println("Failed to write CSV despite DataFrames/CSV being imported:")
        println(e)
    end
else
    println("CSV/DataFrames not available; skipping CSV output.")
end

# ===========================================================================
# Visualizations
# ===========================================================================

println("\n============================================================")
println("GENERATING PLOTS")
println("============================================================")

fmax_vals = Float64[]
for r in results
    if !isnan(r.Fmax) && !isinf(r.Fmax)
        push!(fmax_vals, r.Fmax)
    end
end

if !isempty(fmax_vals)
    # PLOT 1: Fmax distribution
    fig1 = Figure(size=(800, 600))
    ax1 = Axis(fig1[1,1], title="Distribution of maximum possible amplification Fmax",
               xlabel="Fmax", ylabel="Species count")
    hist!(ax1, fmax_vals, bins=50)
    save(joinpath(OUTPUT_DIR, "Fmax_distribution.png"), fig1)

    # PLOT 1b: Fmax distribution (clipped to 99th percentile)
    p99 = quantile(fmax_vals, 0.99)
    fmax_clipped = filter(x -> x <= p99, fmax_vals)
    fig1b = Figure(size=(800, 600))
    ax1b = Axis(fig1b[1,1], title="Distribution of Fmax (clipped to 99th percentile)",
                xlabel="Fmax", ylabel="Species count")
    hist!(ax1b, fmax_clipped, bins=50)
    save(joinpath(OUTPUT_DIR, "Fmax_distribution_clipped.png"), fig1b)
end

finite_collapse_vals = Float64[]
for r in results
    if !isinf(r.min_x_collapse) && !isnan(r.min_x_collapse) && r.min_x_collapse > 0.0
        push!(finite_collapse_vals, log10(r.min_x_collapse))
    end
end

if !isempty(finite_collapse_vals)
    # PLOT 2: min collapse burden distribution
    fig2 = Figure(size=(800, 600))
    ax2 = Axis(fig2[1,1], title="Distribution of minimum DEB-axis burden required for A collapse",
               xlabel="log10(min_x_collapse)", ylabel="Species count")
    hist!(ax2, finite_collapse_vals, bins=50)
    save(joinpath(OUTPUT_DIR, "min_x_collapse_distribution.png"), fig2)
end

# PLOT 3: Fmax vs min_x_collapse scatter
# PLOT 3: Fmax vs min_x_collapse scatter
plot3_data = filter(r -> !isinf(r.min_x_collapse) && !isnan(r.min_x_collapse) && r.min_x_collapse > 0.0 && !isnan(r.Fmax), results)
if !isempty(plot3_data)
    fig3 = Figure(size=(800, 600))
    ax3 = Axis(fig3[1,1], title="Fmax vs min_x_collapse",
               xlabel="log10(min_x_collapse)", ylabel="Fmax")

    xs = [log10(r.min_x_collapse) for r in plot3_data]
    ys = [r.Fmax for r in plot3_data]

    scatter!(ax3, xs, ys, color=(:gray, 0.5), markersize=8)

    # Highlight preferred species
    preferred_colors = [:red, :blue, :green]
    global pref_idx = 0
    for sp_key in preferred
        r_idx = findfirst(r -> r.species_key == sp_key, plot3_data)
        if !isnothing(r_idx)
            global pref_idx += 1
            r = plot3_data[r_idx]
            scatter!(ax3, [log10(r.min_x_collapse)], [r.Fmax],
                     color=preferred_colors[mod1(pref_idx, length(preferred_colors))], markersize=12, label=r.species_name)
        end
    end

    if pref_idx > 0
        axislegend(ax3, position=:rt)
    end

    save(joinpath(OUTPUT_DIR, "Fmax_vs_collapse_burden.png"), fig3)
end


# PLOT 4: lambda(A) curves for selected species
selected_species = String[]
for sp in preferred
    if any(r -> r.species_key == sp, results)
        push!(selected_species, sp)
    end
end
if !isempty(sorted_fmax_desc)
    push!(selected_species, sorted_fmax_desc[1].species_key)
end
if !isempty(sorted_collapse_asc)
    push!(selected_species, sorted_collapse_asc[1].species_key)
end
selected_species = unique(selected_species)

if !isempty(selected_species)
    fig4 = Figure(size=(800, 600))
    ax4 = Axis(fig4[1,1], title="Restoring-force curves for selected species",
               xlabel="A / A0", ylabel="lambda(A)")

    for (i, sp_key) in enumerate(selected_species)
        local params = amp_species_deb_params(amp, sp_key)
        A0 = params.A0
        if A0 > 0.0
            A_vals = range(0.0, A0, length=300)
            A_norm = A_vals ./ A0
            lambda_vals = [restoring_force_from_margin(A, params) for A in A_vals]
            lines!(ax4, A_norm, lambda_vals, label=human_name(sp_key), linewidth=3)
        end
    end
    axislegend(ax4, position=:lt)
    save(joinpath(OUTPUT_DIR, "lambda_curves_selected_species.png"), fig4)
end

# ===========================================================================
# Sanity Checks
# ===========================================================================

println("\n============================================================")
println("SANITY CHECKS")
println("============================================================")

tol = 1e-8
checks_passed = true

for r in results
    # 1. Fmax >= 1
    if !isnan(r.Fmax) && r.Fmax < 1.0 - tol
        @warn "Species $(r.species_key) has Fmax < 1: $(r.Fmax)"
        global checks_passed = false
    end

    # 2. finite x_collapse > 0
    if !isinf(r.min_x_collapse) && !isnan(r.min_x_collapse) && r.min_x_collapse <= 0.0
        @warn "Species $(r.species_key) has min_x_collapse <= 0: $(r.min_x_collapse)"
        global checks_passed = false
    end

    # 3. lambda0 bounds
    if r.lambda0 < r.lambda_min - tol || r.lambda0 > r.lambda_max + tol
        @warn "Species $(r.species_key) has lambda0 outside [lambda_min, lambda_max]: $(r.lambda0)"
        global checks_passed = false
    end
end

for sp in preferred
    if !any(r -> r.species_key == sp, results)
        @warn "Preferred species not found in results: $sp"
        # It's okay if not all are present in the library, but good to know
    end
end

if checks_passed
    println("Analytical response-capacity diagnostics completed.")
end

# ===========================================================================
# Interpretation Block
# ===========================================================================

println("\n============================================================")
println("INTERPRETATION")
println("============================================================")
println("1. Fmax is an analytical upper bound on amplification for a species under the current lambda(A) function.")
println("2. If Fmax is close to 1, no exposure scenario can produce large amplification unless DEBAxisParams are changed.")
println("3. x_collapse_axis tells how much single-axis DEB burden is required to drive A to zero.")
println("4. Large A0 or small alpha_axes imply large collapse thresholds.")
println("5. A0 >> KA often means lambda(A0) is near lambda_max and small A changes may not strongly affect lambda.")
println("6. This diagnostic is not calibration. It identifies which species/parameterizations are structurally capable of strong amplification.")
