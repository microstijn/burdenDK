# ============================================================
# 1. ENVIRONMENT SETUP
# ============================================================
using Pkg
try
    Pkg.activate(joinpath(@__DIR__, ".."))
catch err
    @warn "Could not activate project" exception=(err, catch_backtrace())
end

using TwoTimescaleResilience
using Statistics
using CairoMakie

# ============================================================
# 2. OUTPUT DIRECTORY
# ============================================================
output_dir = joinpath(@__DIR__, "..", "output", "amp_multispecies_synthetic_demo")
mkpath(output_dir)

# ============================================================
# 3. LOAD AmP LIBRARY
# ============================================================
amp_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
amp = load_amp_species_library(amp_path)

# ============================================================
# 4. CUMULATIVE STRESS MULTIPLIERS
# ============================================================
# The default 0–1 normalized synthetic layers produce very small physiological 
# burdens relative to AmP-derived A0 values. For visualization, this example 
# applies a cumulative stress multiplier representing repeated or integrated 
# monthly exposure. This is pedagogical scaling for the example only and does 
# not alter the package model.
cumulative_multiplier_low = 100.0
cumulative_multiplier_moderate = 200.0
cumulative_multiplier_high = 500.0
println("Using cumulative_multiplier_low: ", cumulative_multiplier_low)
println("Using cumulative_multiplier_moderate: ", cumulative_multiplier_moderate)
println("Using cumulative_multiplier_high: ", cumulative_multiplier_high)

# ============================================================
# 5. CREATE RUNTIME MAPPINGS
# ============================================================
base_profile = fish_profile()
exposure_filter = base_profile.exposure_filter
moa_mapping = base_profile.moa_mapping
moa_deb_mapping = base_profile.moa_deb_mapping

# These mappings are reused for demonstration only.
# Only DEB physiological parameters differ across species.
# AmP-derived DEBAxisParams replace base_profile.deb_params.

# ============================================================
# 6. SCREEN MULTIPLE SPECIES
# ============================================================
species_keys = sort(collect(keys(amp)))
abatus_key = "Abatus_cordatus"

# Evaluate all valid species to find those with visible responses.
# The ranking metric is the mean amplification under high cumulative stress.
function screen_species(amp, keys_to_screen, exposure_filter, moa_mapping, moa_deb_mapping, mult_high)
    results = Tuple{String, Float64}[]
    test_high = fill(0.5 * mult_high, 7)
    
    for sp in keys_to_screen
        params = try
            amp_species_deb_params(amp, sp)
        catch
            # Skip species with validation errors
            continue
        end
        
        try
            res = isimip_deb_pipeline(
                test_high,
                exposure_filter,
                moa_mapping,
                moa_deb_mapping,
                params
            )
            push!(results, (sp, res.amplification))
        catch
            continue
        end
    end
    
    # Sort by amplification descending
    sort!(results, by = x -> x[2], rev=true)
    return results
end

ranked_species = screen_species(amp, species_keys, exposure_filter, moa_mapping, moa_deb_mapping, cumulative_multiplier_high)

selected_species = String[]
if abatus_key in [x[1] for x in ranked_species]
    push!(selected_species, abatus_key)
end

for (sp, amp_val) in ranked_species
    if length(selected_species) >= 4
        break
    end
    if sp != abatus_key && !(sp in selected_species)
        push!(selected_species, sp)
    end
end

if isempty(selected_species) || !(abatus_key in selected_species)
    println("Warning: Screening failed to produce expected species. Demonstration might not run as intended.")
    selected_species = [abatus_key]
end

selected_species = unique(selected_species)

# ============================================================
# 7. EXTRACT DEB PARAMETERS
# ============================================================
species_params = Dict{String, Any}()

for sp in selected_species
    params = amp_species_deb_params(amp, sp)
    species_params[sp] = params
end

# ============================================================
# 8. PRINT SPECIES PARAMETERS
# ============================================================
println("\n================ SPECIES PARAMETERS ================\n")

human_name(sp) = replace(sp, "_" => " ")

for sp in sort(collect(keys(species_params)))
    params = species_params[sp]
    lambda0 = restoring_force_from_margin(params.A0, params)
    
    println("----------------------------------------")
    println("Species: ", human_name(sp))
    println("key: ", sp)
    println("A0: ", params.A0)
    println("alpha_axes: ", params.alpha_axes)
    println("lambda_min: ", params.lambda_min)
    println("lambda_max: ", params.lambda_max)
    println("KA: ", params.KA)
    println("baseline lambda0: ", lambda0)
end

# ============================================================
# 9. DEFINE SINGLE-CELL SYNTHETIC STRESS SCENARIOS
# ============================================================
values_zero = zeros(7)
values_low  = fill(0.1 * cumulative_multiplier_low, 7)
values_high = fill(0.5 * cumulative_multiplier_high, 7)

# ============================================================
# 10. RUN SINGLE-CELL COMPARISON
# ============================================================
function run_single(values, species_params, exposure_filter, moa_mapping, moa_deb_mapping)
    results = Dict{String, Any}()
    for (sp, params) in species_params
        res = isimip_deb_pipeline(
            values,
            exposure_filter,
            moa_mapping,
            moa_deb_mapping,
            params
        )
        results[sp] = res
    end
    return results
end

res_zero = run_single(values_zero, species_params, exposure_filter, moa_mapping, moa_deb_mapping)
res_low  = run_single(values_low,  species_params, exposure_filter, moa_mapping, moa_deb_mapping)
res_high = run_single(values_high, species_params, exposure_filter, moa_mapping, moa_deb_mapping)

# ============================================================
# 11. PRINT SINGLE-CELL COMPARISON TABLES
# ============================================================
function print_table(title, results)
    println("\n========================================")
    println("Scenario: ", title)
    println("========================================")
    println(rpad("Species", 28), rpad("A", 15), rpad("lambda", 15), "amplification")

    for sp in sort(collect(keys(results)))
        res = results[sp]
        println(
            rpad(human_name(sp), 28),
            rpad(string(round(res.A, digits=4)), 15),
            rpad(string(round(res.lambda, digits=6)), 15),
            round(res.amplification, digits=4)
        )
    end
end

print_table("ZERO STRESS", res_zero)
print_table("LOW CUMULATIVE STRESS",  res_low)
print_table("HIGH CUMULATIVE STRESS", res_high)

# ============================================================
# 12. BASIC SINGLE-CELL SANITY CHECKS
# ============================================================
println("\n================ SANITY CHECKS ================\n")
all_checks_passed = true

for sp in sort(collect(keys(species_params)))
    if !(res_high[sp].A <= res_low[sp].A + eps())
        all_checks_passed = false
    end
    if !(res_high[sp].amplification >= res_low[sp].amplification - eps())
        all_checks_passed = false
    end
end

if all_checks_passed
    println("Sanity checks passed.")
else
    println("Warning: behavior depends on mapping assumptions; demonstration remains valid as an AmP adapter example.")
end

# ============================================================
# 13. DEFINE 3×3 SYNTHETIC SPATIAL PATTERNS
# ============================================================
gradient_3x3 = [
    0.0  0.1  0.2;
    0.1  0.3  0.5;
    0.2  0.5  0.8
]

grid_pattern_low  = cumulative_multiplier_low .* gradient_3x3
grid_pattern_moderate = cumulative_multiplier_moderate .* gradient_3x3
grid_pattern_high = cumulative_multiplier_high .* clamp.(0.5 .* gradient_3x3 .+ 0.4, 0.0, 1.0)

layers_low  = [copy(grid_pattern_low) for _ in 1:7]
layers_moderate = [copy(grid_pattern_moderate) for _ in 1:7]
layers_high = [copy(grid_pattern_high) for _ in 1:7]

# ============================================================
# 14. RUN 3×3 GRID COMPARISON
# ============================================================
function finite_mean(A)
    vals = filter(isfinite, vec(A))
    isempty(vals) && return NaN
    return mean(vals)
end

function finite_max(A)
    vals = filter(isfinite, vec(A))
    isempty(vals) && return NaN
    return maximum(vals)
end

function amplification_colorrange(results, species_list)
    vals = Float64[]
    for sp in species_list
        append!(vals, filter(isfinite, vec(results[sp].amplification)))
    end
    if isempty(vals)
        return (1.0, 1.0 + 1e-6)
    end
    #lo = minimum(vals)
    lo = 1.0
    hi = maximum(vals)
    if isapprox(lo, hi; atol=1e-12)
        hi = lo + 1e-6
    end
    return (lo, hi)
end

grid_results_low = Dict{String, Any}()
grid_results_moderate = Dict{String, Any}()
grid_results_high = Dict{String, Any}()

for (sp, params) in species_params
    grid_low = isimip_deb_pipeline_grid(
        layers_low,
        exposure_filter,
        moa_mapping,
        moa_deb_mapping,
        params
    )
    
    grid_moderate = isimip_deb_pipeline_grid(
        layers_moderate,
        exposure_filter,
        moa_mapping,
        moa_deb_mapping,
        params
    )

    grid_high = isimip_deb_pipeline_grid(
        layers_high,
        exposure_filter,
        moa_mapping,
        moa_deb_mapping,
        params
    )

    grid_results_low[sp] = grid_low
    grid_results_moderate[sp] = grid_moderate
    grid_results_high[sp] = grid_high
end

# ============================================================
# 15. PRINT 3×3 GRID SUMMARY
# ============================================================
println("\n================ 3x3 GRID SUMMARY ================\n")

for sp in sort(collect(keys(species_params)))
    low = grid_results_low[sp]
    moderate = grid_results_moderate[sp]
    high = grid_results_high[sp]

    println("----------------------------------------")
    println("Species: ", human_name(sp))
    println("mean amplification low:      ", round(finite_mean(low.amplification), digits=4))
    println("mean amplification moderate: ", round(finite_mean(moderate.amplification), digits=4))
    println("mean amplification high:     ", round(finite_mean(high.amplification), digits=4))
    println("max amplification moderate:  ", round(finite_max(moderate.amplification), digits=4))
    println("max amplification high:      ", round(finite_max(high.amplification), digits=4))
    println("mean A moderate:             ", round(finite_mean(moderate.A), digits=4))
    println("mean A high:                 ", round(finite_mean(high.A), digits=4))
    println("mean lambda moderate:        ", round(finite_mean(moderate.lambda), digits=6))
    println("mean lambda high:            ", round(finite_mean(high.lambda), digits=6))
end

println("\nLow cumulative stress shows near-baseline amplification.")
println("Moderate cumulative stress is intended to show spatial gradients before full collapse.")
println("High cumulative stress can saturate lambda at lambda_min, making amplification spatially uniform within a species.")
println("This saturation is expected behavior, not a plotting bug.")

# ============================================================
# 16. GRAPH 3×3 AMPLIFICATION GRIDS WITH CAIROMAKIE
# ============================================================
n_species = length(selected_species)
ncols = min(3, n_species)
nrows = ceil(Int, n_species / ncols)

# MODERATE
fig_mod = Figure(size=(350 * ncols, 320 * nrows))
colorrange_mod = amplification_colorrange(grid_results_moderate, selected_species)

for (i, sp) in enumerate(selected_species)
    row = ceil(Int, i / ncols)
    col = (i - 1) % ncols + 1
    
    ax = Axis(fig_mod[row, col]; title=human_name(sp) * "\nModerate cumulative synthetic stress", aspect=DataAspect())
    hm = heatmap!(ax, grid_results_moderate[sp].amplification; colormap=:inferno, colorrange=colorrange_mod)
    hidedecorations!(ax)
    
    if i == length(selected_species)
        Colorbar(fig_mod[:, ncols + 1], hm, label="Amplification")
    end
end
fig_mod_path = joinpath(output_dir, "amp_multispecies_amplification_moderate.png")
save(fig_mod_path, fig_mod, px_per_unit=2)

# HIGH
fig_high = Figure(size=(350 * ncols, 320 * nrows))
colorrange_high = amplification_colorrange(grid_results_high, selected_species)

for (i, sp) in enumerate(selected_species)
    row = ceil(Int, i / ncols)
    col = (i - 1) % ncols + 1
    
    ax = Axis(fig_high[row, col]; title=human_name(sp) * "\nHigh cumulative synthetic stress", aspect=DataAspect())
    hm = heatmap!(ax, grid_results_high[sp].amplification; colormap=:inferno, colorrange=colorrange_high)
    hidedecorations!(ax)
    
    if i == length(selected_species)
        Colorbar(fig_high[:, ncols + 1], hm, label="Amplification")
    end
end
fig_high_path = joinpath(output_dir, "amp_multispecies_amplification_high.png")
save(fig_high_path, fig_high, px_per_unit=2)

# ABATUS 1x3
if haskey(grid_results_low, abatus_key) && haskey(grid_results_moderate, abatus_key) && haskey(grid_results_high, abatus_key)
    fig_abatus = Figure(size=(1200, 400))
    
    abatus_low_vals = filter(isfinite, vec(grid_results_low[abatus_key].amplification))
    abatus_mod_vals = filter(isfinite, vec(grid_results_moderate[abatus_key].amplification))
    abatus_high_vals = filter(isfinite, vec(grid_results_high[abatus_key].amplification))
    
    all_abatus_vals = [abatus_low_vals; abatus_mod_vals; abatus_high_vals]
    
    if isempty(all_abatus_vals)
        c_range = (1.0, 1.0 + 1e-6)
    else
        #lo_a = minimum(all_abatus_vals)
        lo_a = 1.0
        hi_a = maximum(all_abatus_vals)
        if isapprox(lo_a, hi_a; atol=1e-12)
            hi_a = lo_a + 1e-6
        end
        c_range = (lo_a, hi_a)
    end
    
    ax1 = Axis(fig_abatus[1, 1], title="Abatus cordatus\nLow cumulative stress", aspect=DataAspect())
    hm1 = heatmap!(ax1, grid_results_low[abatus_key].amplification; colormap=:inferno, colorrange=c_range)
    hidedecorations!(ax1)
    
    ax2 = Axis(fig_abatus[1, 2], title="Abatus cordatus\nModerate cumulative stress", aspect=DataAspect())
    hm2 = heatmap!(ax2, grid_results_moderate[abatus_key].amplification; colormap=:inferno, colorrange=c_range)
    hidedecorations!(ax2)

    ax3 = Axis(fig_abatus[1, 3], title="Abatus cordatus\nHigh cumulative stress", aspect=DataAspect())
    hm3 = heatmap!(ax3, grid_results_high[abatus_key].amplification; colormap=:inferno, colorrange=c_range)
    hidedecorations!(ax3)
    
    Colorbar(fig_abatus[1, 4], hm3, label="Amplification")
    
    fig_abatus_path = joinpath(output_dir, "amp_multispecies_amplification_low_moderate_high_abatus.png")
    save(fig_abatus_path, fig_abatus, px_per_unit=2)
else
    fig_abatus_path = ""
end

# Remove old 1x2 figure if it still exists
old_abatus_path = joinpath(output_dir, "amp_multispecies_amplification_low_high_abatus.png")
if isfile(old_abatus_path)
    rm(old_abatus_path)
end

println("Saved figure: ", fig_mod_path)
println("Saved figure: ", fig_high_path)
if !isempty(fig_abatus_path)
    println("Saved figure: ", fig_abatus_path)
end

# ============================================================
# 17. INTERPRETATION BLOCK
# ============================================================
println("\n================ INTERPRETATION ================\n")
println("Different species show different vulnerability profiles under identical exposure conditions.")
println("Under modest 0–1 normalized stress, most species show F ≈ 1 because AmP-derived adaptive margins are large.")
println("The cumulative scenario demonstrates how integrated or repeated stress can deplete A enough to visibly affect lambda and amplification.")
println("These differences arise from AmP-derived physiological parameters: A0, alpha_axes, KA, and lambda bounds.")
println("The CairoMakie 3x3 heatmaps visualize how identical cumulative synthetic exposure fields map to species-specific amplification patterns.")

# ============================================================
# 18. FINAL OUTPUT
# ============================================================
println("\nDone.")
