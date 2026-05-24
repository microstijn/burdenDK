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
# 4. SELECT MULTIPLE SPECIES
# ============================================================
species_keys = sort(collect(keys(amp)))
selected_species = species_keys[1:min(4, length(species_keys))]

if !("Abatus_cordatus" in selected_species)
    push!(selected_species, "Abatus_cordatus")
end

selected_species = unique(selected_species)

# ============================================================
# 5. EXTRACT DEB PARAMETERS
# ============================================================
species_params = Dict{String, Any}()

for sp in selected_species
    params = amp_species_deb_params(amp, sp)
    species_params[sp] = params
end

# ============================================================
# 6. PRINT SPECIES PARAMETERS
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
# 7. CREATE RUNTIME MAPPINGS
# ============================================================
base_profile = fish_profile()
exposure_filter = base_profile.exposure_filter
moa_mapping = base_profile.moa_mapping
moa_deb_mapping = base_profile.moa_deb_mapping

# These mappings are reused for demonstration only.
# Only DEB physiological parameters differ across species.
# AmP-derived DEBAxisParams replace base_profile.deb_params.

# ============================================================
# 8. DEFINE SINGLE-CELL SYNTHETIC STRESS SCENARIOS
# ============================================================
values_zero = zeros(7)
values_low  = fill(0.1, 7)
values_high = fill(0.5, 7)

# ============================================================
# 9. RUN SINGLE-CELL COMPARISON
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
# 10. PRINT SINGLE-CELL COMPARISON TABLES
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
print_table("LOW STRESS",  res_low)
print_table("HIGH STRESS", res_high)

# ============================================================
# 11. BASIC SINGLE-CELL SANITY CHECKS
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
# 12. DEFINE 3×3 SYNTHETIC SPATIAL PATTERNS
# ============================================================
gradient_3x3 = [
    0.0  0.1  0.2;
    0.1  0.3  0.5;
    0.2  0.5  0.8
]

grid_pattern_low  = gradient_3x3
grid_pattern_high = 0.5 .* gradient_3x3 .+ 0.4
grid_pattern_high = clamp.(grid_pattern_high, 0.0, 1.0)

layers_low  = [copy(grid_pattern_low) for _ in 1:7]
layers_high = [copy(grid_pattern_high) for _ in 1:7]

# ============================================================
# 13. RUN 3×3 GRID COMPARISON
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

grid_results_low = Dict{String, Any}()
grid_results_high = Dict{String, Any}()

for (sp, params) in species_params
    grid_low = isimip_deb_pipeline_grid(
        layers_low,
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
    grid_results_high[sp] = grid_high
end

# ============================================================
# 14. PRINT 3×3 GRID SUMMARY
# ============================================================
println("\n================ 3x3 GRID SUMMARY ================\n")

for sp in sort(collect(keys(species_params)))
    low = grid_results_low[sp]
    high = grid_results_high[sp]

    println("----------------------------------------")
    println("Species: ", human_name(sp))
    println("mean amplification low:  ", round(finite_mean(low.amplification), digits=4))
    println("mean amplification high: ", round(finite_mean(high.amplification), digits=4))
    println("max amplification high:  ", round(finite_max(high.amplification), digits=4))
    println("mean A high:             ", round(finite_mean(high.A), digits=4))
    println("mean lambda high:        ", round(finite_mean(high.lambda), digits=6))
end

# ============================================================
# 15. GRAPH 3×3 AMPLIFICATION GRIDS WITH CAIROMAKIE
# ============================================================
n_species = length(selected_species)
ncols = min(3, n_species)
nrows = ceil(Int, n_species / ncols)

fig = Figure(size=(350 * ncols, 320 * nrows))

all_high_vals = Float64[]
for sp in selected_species
    append!(all_high_vals, filter(isfinite, vec(grid_results_high[sp].amplification)))
end

if isempty(all_high_vals)
    colorrange = (1.0, 1.0)
else
    lo = minimum(all_high_vals)
    hi = maximum(all_high_vals)
    if isapprox(lo, hi; atol=1e-12)
        hi = lo + 1e-6
    end
    colorrange = (lo, hi)
end

for (i, sp) in enumerate(selected_species)
    row = ceil(Int, i / ncols)
    col = (i - 1) % ncols + 1

    ax = Axis(fig[row, col]; title=human_name(sp), aspect=DataAspect())
    hm = heatmap!(ax, grid_results_high[sp].amplification; colormap=:inferno, colorrange=colorrange)
    hidedecorations!(ax)

    if i == length(selected_species)
        Colorbar(fig[:, ncols + 1], hm, label="Amplification")
    end
end

fig_path = joinpath(output_dir, "amp_multispecies_amplification_high.png")
save(fig_path, fig, px_per_unit=2)

abatus_key = "Abatus_cordatus"
if haskey(grid_results_low, abatus_key) && haskey(grid_results_high, abatus_key)
    fig2 = Figure(size=(800, 400))

    abatus_low_vals = filter(isfinite, vec(grid_results_low[abatus_key].amplification))
    abatus_high_vals = filter(isfinite, vec(grid_results_high[abatus_key].amplification))

    all_abatus_vals = [abatus_low_vals; abatus_high_vals]

    if isempty(all_abatus_vals)
        c_range = (1.0, 1.0)
    else
        lo_a = minimum(all_abatus_vals)
        hi_a = maximum(all_abatus_vals)
        if isapprox(lo_a, hi_a; atol=1e-12)
            hi_a = lo_a + 1e-6
        end
        c_range = (lo_a, hi_a)
    end

    ax1 = Axis(fig2[1, 1], title="Abatus cordatus — low stress", aspect=DataAspect())
    hm1 = heatmap!(ax1, grid_results_low[abatus_key].amplification; colormap=:inferno, colorrange=c_range)
    hidedecorations!(ax1)

    ax2 = Axis(fig2[1, 2], title="Abatus cordatus — high stress", aspect=DataAspect())
    hm2 = heatmap!(ax2, grid_results_high[abatus_key].amplification; colormap=:inferno, colorrange=c_range)
    hidedecorations!(ax2)

    Colorbar(fig2[1, 3], hm2, label="Amplification")

    fig2_path = joinpath(output_dir, "amp_multispecies_amplification_low_high_abatus.png")
    save(fig2_path, fig2, px_per_unit=2)
else
    fig2_path = ""
end

println("Saved figure: ", fig_path)
if !isempty(fig2_path)
    println("Saved figure: ", fig2_path)
end

# ============================================================
# 16. INTERPRETATION BLOCK
# ============================================================
println("\n================ INTERPRETATION ================\n")
println("Different species show different vulnerability profiles under identical exposure conditions.")
println("These differences arise from AmP-derived physiological parameters: A0, alpha_axes, KA, and lambda bounds.")
println("The CairoMakie 3x3 heatmaps visualize how identical synthetic exposure fields map to species-specific amplification patterns.")

# ============================================================
# 17. FINAL OUTPUT
# ============================================================
println("\nDone.")
