# Deterministic Synthetic Grid Mixture Demo
#
# This is a deterministic synthetic NetCDF grid example.
# Concentrations are synthetic spatial patterns scaled using actual ECOTOX NOEC/EC50 records.
# Selected compounds are actual records from the local ECOTOX library.
# Selected species are actual records from the local AmP library.
# The example writes NetCDF input grids and reloads them before model evaluation.
#
# This is NOT real raster integration.
# This does NOT add ECOTOX records.
# This does NOT implement synergism or antagonism.
# This does NOT implement fitted interaction coefficients.
# This does NOT implement DEBtox scaled damage.
# This does NOT implement physiological Z_t.
# This is a stepping stone before real raster integration.

using TwoTimescaleResilience
using NCDatasets
using CairoMakie
using CSV
using DataFrames
using Statistics

println("Starting Deterministic Synthetic Grid Mixture Demo...")

# Get environment variables or defaults
nx = parse(Int, get(ENV, "TTR_GRID_NX", "160"))
ny = parse(Int, get(ENV, "TTR_GRID_NY", "120"))
n_months = 12
req_n_species = parse(Int, get(ENV, "TTR_GRID_N_SPECIES", "12"))
req_n_compounds = parse(Int, get(ENV, "TTR_GRID_N_COMPOUNDS", "12"))

out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_synthetic_grid_mixture_demo"))
mkpath(out_dir)

# -----------------------------------------------------------------------------
# 1. Species Selection
# -----------------------------------------------------------------------------
println("Selecting Species...")
amp_library = load_amp_species_library()

selected_species = []
amp_pairs = sort(collect(amp_library), by=x->x[1])
for (i, (key, rec)) in enumerate(amp_pairs)
    params = try
        amp_record_to_deb_params(rec)
    catch
        nothing
    end
    if params !== nothing
        w = axis_weights_for_species(params)
        lambda0 = restoring_force_from_margin(params.A0, params)
        lambda_min = params.KA
        lambda_max = params.KA + params.A0
        Fmax = lambda0 / lambda_min

        push!(selected_species, (
            species_key = get(rec, "species_name", "unknown_key"),
            species_name = replace(get(rec, "species_name", "unknown_key"), "_" => " "),
            A0 = params.A0,
            lambda_min = lambda_min,
            lambda_max = lambda_max,
            KA = params.KA,
            axis_weight_method = w.axis_weight_method,
            axis_weight_scope = w.axis_weight_scope,
            lambda0 = lambda0,
            Fmax = Fmax,
            params = params
        ))
    end
    if length(selected_species) >= req_n_species
        break
    end
end

if length(selected_species) < 10
    error("Failed to select at least 10 valid species.")
end

df_species = DataFrame([NamedTuple{Tuple(keys(x)[1:10])}(values(x)[1:10]) for x in selected_species])
CSV.write(joinpath(out_dir, "selected_species.csv"), df_species)
println("Selected $(length(selected_species)) species.")


# -----------------------------------------------------------------------------
# 2. Compound Selection
# -----------------------------------------------------------------------------
println("Selecting Compounds...")
ecotox_library = load_ecotox_library()
memory_library = load_compound_memory_library()

# Function to get norm cas
function get_cas_norm(c)
    if ismissing(c) || isnothing(c)
        return ""
    end
    return replace(string(c), r"[^\d]" => "")
end

valid_compounds = []
for rec in ecotox_library
    noec = get(rec, "NOEC_median", missing)
    ec50 = get(rec, "EC50_median", missing)

    if ismissing(noec) || ismissing(ec50) || isnothing(noec) || isnothing(ec50) || !isfinite(noec) || !isfinite(ec50) || ec50 <= noec
        continue
    end

    try
        axis = ecotox_effect_to_deb_axis(get(rec, "effect_code", ""))
        rho = compound_retention(get_cas_norm(get(rec, "cas_norm", "")); memory_library=memory_library)
        K = compound_bioaccumulation_factor(get_cas_norm(get(rec, "cas_norm", "")); memory_library=memory_library)
        push!(valid_compounds, (
            record = rec,
            axis = axis,
            effect_code = uppercase(strip(get(rec, "effect_code", ""))),
            cas_norm = get_cas_norm(get(rec, "cas_norm", "")),
            cas_hyphenated = get(rec, "cas_hyphenated", ""),
            chemical_name = get(rec, "chemical_name", "Unknown"),
            NOEC_median = noec,
            EC50_median = ec50,
            retention_rho_monthly = rho,
            bioaccumulation_factor = K
        ))
    catch
        continue
    end
end

# Sort deterministically
sort!(valid_compounds, by = x -> (string(x.axis), string(x.effect_code), string(x.cas_norm)))

selected_compounds = []
axis_counts = Dict{Symbol, Int}(:assimilation=>0, :maintenance=>0, :growth=>0, :reproduction=>0)
effect_codes_by_axis = Dict{Symbol, Set{String}}(:assimilation=>Set(), :maintenance=>Set(), :growth=>Set(), :reproduction=>Set())

# Try to select matching the ideal profile: >=3 M, >=2 G, >=2 R
for c in valid_compounds
    if c.axis == :maintenance && axis_counts[:maintenance] < 4
        push!(selected_compounds, c)
        axis_counts[:maintenance] += 1
        push!(effect_codes_by_axis[:maintenance], c.effect_code)
    elseif c.axis == :growth && axis_counts[:growth] < 3
        push!(selected_compounds, c)
        axis_counts[:growth] += 1
        push!(effect_codes_by_axis[:growth], c.effect_code)
    elseif c.axis == :reproduction && axis_counts[:reproduction] < 3
        push!(selected_compounds, c)
        axis_counts[:reproduction] += 1
        push!(effect_codes_by_axis[:reproduction], c.effect_code)
    elseif c.axis == :assimilation && axis_counts[:assimilation] < 2
        push!(selected_compounds, c)
        axis_counts[:assimilation] += 1
        push!(effect_codes_by_axis[:assimilation], c.effect_code)
    end
    if length(selected_compounds) >= req_n_compounds && axis_counts[:maintenance] >= 2 && length(effect_codes_by_axis[:maintenance]) >= 2
        break
    end
end

# Fill up to requested number if short
if length(selected_compounds) < req_n_compounds
    for c in valid_compounds
        if !(c in selected_compounds)
            push!(selected_compounds, c)
            axis_counts[c.axis] += 1
            push!(effect_codes_by_axis[c.axis], c.effect_code)
        end
        if length(selected_compounds) >= req_n_compounds
            break
        end
    end
end

if length(selected_compounds) < 10
    error("Failed to select at least 10 valid compounds.")
end

global has_same_axis_overlap = false
for (axis, count) in axis_counts
    if count >= 2
        global has_same_axis_overlap = true
    end
end

if !has_same_axis_overlap
    error("Failed to select at least two compounds on the same axis.")
end

df_compounds = DataFrame([
    (
        compound_index = i,
        cas_norm = c.cas_norm,
        cas_hyphenated = c.cas_hyphenated,
        chemical_name = c.chemical_name,
        effect_code = c.effect_code,
        deb_axis = string(c.axis),
        NOEC_median = c.NOEC_median,
        EC50_median = c.EC50_median,
        retention_rho_monthly = c.retention_rho_monthly,
        bioaccumulation_factor = c.bioaccumulation_factor
    ) for (i, c) in enumerate(selected_compounds)
])
CSV.write(joinpath(out_dir, "selected_compounds.csv"), df_compounds)
println("Selected $(length(selected_compounds)) compounds.")


# -----------------------------------------------------------------------------
# 3. Create Synthetic Grid Patterns
# -----------------------------------------------------------------------------
println("Creating Synthetic Concentration Grid...")

n_c = length(selected_compounds)
C_t = zeros(Float32, nx, ny, n_months, n_c)

# Helper functions for patterns
dist(x1, y1, x2, y2) = sqrt((x1 - x2)^2 + (y1 - y2)^2)

for (ic, comp) in enumerate(selected_compounds)
    pattern_type = mod(ic - 1, 4) + 1

    for m in 1:n_months
        for y in 1:ny
            for x in 1:nx
                I = 0.0

                if pattern_type == 1
                    # 1. River corridor / downstream plume (diagonal)
                    river_dist = abs(x - y * (nx / ny)) / sqrt(1 + (nx / ny)^2)
                    I = exp(-river_dist / 10.0) * (1.0 + 0.5 * sin(2 * pi * m / 12))

                elseif pattern_type == 2
                    # 2. Diffuse agricultural runoff (higher in top right, spring/summer pulse)
                    grad_x = x / nx
                    grad_y = y / ny
                    season_pulse = m in 4:8 ? 1.5 : 0.2
                    I = (grad_x * grad_y) * season_pulse

                elseif pattern_type == 3
                    # 3. Urban wastewater hotspot (center)
                    cx, cy = nx / 2, ny / 2
                    d = dist(x, y, cx, cy)
                    I = 2.0 * exp(-(d^2) / 200.0)

                elseif pattern_type == 4
                    # 4. Industrial point-source (episodic summer peak)
                    cx, cy = nx / 4, ny * 3 / 4
                    d = dist(x, y, cx, cy)
                    episodic = m in 6:8 ? 3.0 : 0.5
                    I = episodic * exp(-(d^2) / 100.0)
                end

                # Background baseline 0.1 * NOEC
                baseline = 0.1

                # Scale using NOEC and EC50
                # C = max(0, NOEC * baseline + I * (EC50 - NOEC))
                # So I=0 gives baseline, I=1 gives ~EC50, I>1 goes above EC50
                NOEC = comp.NOEC_median
                EC50 = comp.EC50_median
                val = max(0.0, NOEC * baseline + I * (EC50 - NOEC))

                C_t[x, y, m, ic] = val
            end
        end
    end
end

# Write Input NetCDF
nc_in_path = joinpath(out_dir, "synthetic_compound_concentration_grid.nc")
if isfile(nc_in_path) rm(nc_in_path) end

NCDataset(nc_in_path, "c") do ds
    defDim(ds, "x", nx)
    defDim(ds, "y", ny)
    defDim(ds, "month", n_months)
    defDim(ds, "compound", n_c)

    v = defVar(ds, "C_t", Float32, ("x", "y", "month", "compound"))
    v.attrib["units"] = "same concentration units as selected ECOTOX NOEC/EC50 records"

    ds.attrib["scenario"] = "synthetic_spatial_compound_grid"
    ds.attrib["generated_by"] = "examples/ecotox_amp_synthetic_grid_mixture_demo.jl"
    ds.attrib["synthetic"] = "true"
    ds.attrib["deterministic"] = "true"
    ds.attrib["dimension_order"] = "[x, y, month, compound]"

    v[:, :, :, :] = C_t
end
println("Wrote $nc_in_path")


# -----------------------------------------------------------------------------
# 4. Load from Disk & Run Model Computation
# -----------------------------------------------------------------------------
println("Loading Grids and Computing Reponses...")

C_t_loaded = NCDataset(nc_in_path, "r") do ds
    Array(ds["C_t"])
end

models = [
    "axis_toxic_unit_sum",
    "independent_action_axis_effects",
    "grouped_ca_then_ia_axis_effects"
]
n_models = length(models)
n_s = length(selected_species)

# Outputs
Q_out = zeros(Float32, nx, ny, n_months, n_s, n_models)
A_out = zeros(Float32, nx, ny, n_months, n_s, n_models)
lambda_out = zeros(Float32, nx, ny, n_months, n_s, n_models)
F_out = zeros(Float32, nx, ny, n_months, n_s, n_models)

# Axis outputs
E_A_out = zeros(Float32, nx, ny, n_months, n_models)
E_M_out = zeros(Float32, nx, ny, n_months, n_models)
E_G_out = zeros(Float32, nx, ny, n_months, n_models)
E_R_out = zeros(Float32, nx, ny, n_months, n_models)

X_A_out = zeros(Float32, nx, ny, n_months)
X_M_out = zeros(Float32, nx, ny, n_months)
X_G_out = zeros(Float32, nx, ny, n_months)
X_R_out = zeros(Float32, nx, ny, n_months)

active_cells = zeros(Int, n_months, 4)
multi_cells = zeros(Int, n_months, 4)

for y in 1:ny
    for x in 1:nx
        state = EcotoxExposureState()

        for m in 1:n_months
            concentrations = Dict{String, Float64}()
            for ic in 1:n_c
                concentrations[selected_compounds[ic].cas_norm] = C_t_loaded[x, y, m, ic]
            end

            # Since we need per-compound burdens for grouped mixture arithmetic,
            # we run update_internal_burden! and then ecotox_record_to_deb_burden manually

            compound_records_with_burdens = []

            for ic in 1:n_c
                comp = selected_compounds[ic]
                cas_norm = comp.cas_norm
                C_val = concentrations[cas_norm]

                # Update burden
                B_val = update_internal_burden!(state, cas_norm, C_val;
                            retention = comp.retention_rho_monthly,
                            bioaccumulation_factor = comp.bioaccumulation_factor)

                # Add to record list for mixture
                # Create a pseudo-record dict for the mixture aggregator
                rec_dict = Dict(
                    "cas_norm" => comp.cas_norm, "cas" => comp.cas_norm, "cas_hyphenated" => comp.cas_hyphenated,
                    "effect_code" => comp.effect_code,
                    "taxon_class" => "dummy", "n_NOEC" => 1, "n_EC50" => 1, "NOEC_median" => comp.NOEC_median,
                    "EC50_median" => comp.EC50_median
                )

                # To simulate ecotox_records_to_deb_burden_stateful! but keep the record level separation
                # we just pass a list of dicts to aggregate_axis_mixture_effects
                # Wait, aggregate_axis_mixture_effects expects records to have burden_assimilation etc.

                b_axes = ecotox_record_to_deb_burden(B_val, rec_dict)

                push!(compound_records_with_burdens, (
                    burden_assimilation = b_axes.assimilation,
                    burden_maintenance = b_axes.maintenance,
                    burden_growth = b_axes.growth,
                    burden_reproduction = b_axes.reproduction,
                    effect_code = comp.effect_code,
                    cas_norm = comp.cas_norm
                ))
            end

            # Now run for each model
            for (im, model_name) in enumerate(models)
                res = aggregate_axis_mixture_effects(compound_records_with_burdens; mixture_effect_model = model_name)

                E_A_out[x, y, m, im] = res.E_assimilation
                E_M_out[x, y, m, im] = res.E_maintenance
                E_G_out[x, y, m, im] = res.E_growth
                E_R_out[x, y, m, im] = res.E_reproduction

                if im == 1
                    X_A_out[x, y, m] = res.X_assimilation
                    X_M_out[x, y, m] = res.X_maintenance
                    X_G_out[x, y, m] = res.X_growth
                    X_R_out[x, y, m] = res.X_reproduction

                    if res.X_assimilation > 0 active_cells[m, 1] += 1 end
                    if res.X_maintenance > 0 active_cells[m, 2] += 1 end
                    if res.X_growth > 0 active_cells[m, 3] += 1 end
                    if res.X_reproduction > 0 active_cells[m, 4] += 1 end

                    # Estimate multi-contributor trivially
                    c_a = sum([r.burden_assimilation > 0 for r in compound_records_with_burdens])
                    c_m = sum([r.burden_maintenance > 0 for r in compound_records_with_burdens])
                    c_g = sum([r.burden_growth > 0 for r in compound_records_with_burdens])
                    c_r = sum([r.burden_reproduction > 0 for r in compound_records_with_burdens])

                    if c_a > 1 multi_cells[m, 1] += 1 end
                    if c_m > 1 multi_cells[m, 2] += 1 end
                    if c_g > 1 multi_cells[m, 3] += 1 end
                    if c_r > 1 multi_cells[m, 4] += 1 end
                end

                # Now species response
                for (is, sp) in enumerate(selected_species)
                    impairments = (
                        assimilation = res.E_assimilation,
                        maintenance = res.E_maintenance,
                        growth = res.E_growth,
                        reproduction = res.E_reproduction
                    )

                    # compute_adaptive_margin_response_from_impairment
                    sp_res = compute_adaptive_margin_response_from_impairment(impairments, sp.params;
                        X_axis = (assimilation=res.X_assimilation, maintenance=res.X_maintenance, growth=res.X_growth, reproduction=res.X_reproduction),
                        mixture_effect_model = model_name)

                    Q_out[x, y, m, is, im] = sp_res.Q_t
                    A_out[x, y, m, is, im] = sp_res.A_t
                    lambda_out[x, y, m, is, im] = sp_res.lambda_t
                    F_out[x, y, m, is, im] = sp_res.F_t
                end
            end
        end
    end
end

# -----------------------------------------------------------------------------
# 5. Write Outputs to NetCDF
# -----------------------------------------------------------------------------
println("Writing Output NetCDF...")
nc_out_path = joinpath(out_dir, "synthetic_grid_response_outputs.nc")
if isfile(nc_out_path) rm(nc_out_path) end

NCDataset(nc_out_path, "c") do ds
    defDim(ds, "x", nx)
    defDim(ds, "y", ny)
    defDim(ds, "month", n_months)
    defDim(ds, "species", n_s)
    defDim(ds, "mixture_effect_model", n_models)

    ds.attrib["dimension_order"] = "[x, y, month, species, mixture_effect_model]"

    v_Q = defVar(ds, "Q_t", Float32, ("x", "y", "month", "species", "mixture_effect_model"))
    v_A = defVar(ds, "A_t", Float32, ("x", "y", "month", "species", "mixture_effect_model"))
    v_lam = defVar(ds, "lambda_t", Float32, ("x", "y", "month", "species", "mixture_effect_model"))
    v_F = defVar(ds, "F_t", Float32, ("x", "y", "month", "species", "mixture_effect_model"))

    v_Q[:, :, :, :, :] = Q_out
    v_A[:, :, :, :, :] = A_out
    v_lam[:, :, :, :, :] = lambda_out
    v_F[:, :, :, :, :] = F_out
end
println("Wrote $nc_out_path")


# -----------------------------------------------------------------------------
# 6. CSV Summaries
# -----------------------------------------------------------------------------
println("Generating CSV Summaries...")

function safe_quantile(arr, q)
    vals = filter(isfinite, arr)
    if isempty(vals) return NaN end
    return quantile(vals, q)
end

function safe_mean(arr)
    vals = filter(isfinite, arr)
    if isempty(vals) return NaN end
    return mean(vals)
end

function safe_max(arr)
    vals = filter(isfinite, arr)
    if isempty(vals) return NaN end
    return maximum(vals)
end

function safe_min(arr)
    vals = filter(isfinite, arr)
    if isempty(vals) return NaN end
    return minimum(vals)
end

# Species summary
species_summary_rows = []
for m in 1:n_months
    for is in 1:n_s
        for im in 1:n_models
            Q_slice = Q_out[:, :, m, is, im]
            A_slice = A_out[:, :, m, is, im]
            F_slice = F_out[:, :, m, is, im]

            n_cells = length(Q_slice)

            push!(species_summary_rows, (
                species_key = selected_species[is].species_key,
                species_name = selected_species[is].species_name,
                mixture_effect_model = models[im],
                month = m,
                n_cells = n_cells,
                mean_Q_t = safe_mean(Q_slice),
                max_Q_t = safe_max(Q_slice),
                p95_Q_t = safe_quantile(Q_slice, 0.95),
                mean_A_t = safe_mean(A_slice),
                min_A_t = safe_min(A_slice),
                p05_A_t = safe_quantile(A_slice, 0.05),
                mean_F_t = safe_mean(F_slice),
                max_F_t = safe_max(F_slice),
                p95_F_t = safe_quantile(F_slice, 0.95),
                fraction_cells_F_gt_1_01 = sum(F_slice .> 1.01) / n_cells,
                fraction_cells_F_gt_1_10 = sum(F_slice .> 1.10) / n_cells,
                fraction_cells_Q_gt_0_50 = sum(Q_slice .> 0.50) / n_cells,
                fraction_cells_Q_gt_0_80 = sum(Q_slice .> 0.80) / n_cells
            ))
        end
    end
end
CSV.write(joinpath(out_dir, "synthetic_grid_species_summary.csv"), DataFrame(species_summary_rows))

# Axis summary
axis_summary_rows = []
axes_names = ["assimilation", "maintenance", "growth", "reproduction"]
axes_X_outs = [X_A_out, X_M_out, X_G_out, X_R_out]
axes_E_outs = [E_A_out, E_M_out, E_G_out, E_R_out]

for m in 1:n_months
    for (iax, ax_name) in enumerate(axes_names)
        for im in 1:n_models
            X_slice = axes_X_outs[iax][:, :, m]
            E_slice = axes_E_outs[iax][:, :, m, im]

            push!(axis_summary_rows, (
                mixture_effect_model = models[im],
                month = m,
                axis = ax_name,
                mean_X_axis = safe_mean(X_slice),
                max_X_axis = safe_max(X_slice),
                p95_X_axis = safe_quantile(X_slice, 0.95),
                mean_E_axis = safe_mean(E_slice),
                max_E_axis = safe_max(E_slice),
                p95_E_axis = safe_quantile(E_slice, 0.95),
                fraction_active_cells = active_cells[m, iax] / (nx*ny),
                fraction_multi_contributor_cells = multi_cells[m, iax] / (nx*ny)
            ))
        end
    end
end
CSV.write(joinpath(out_dir, "synthetic_grid_axis_summary.csv"), DataFrame(axis_summary_rows))

# Mixture Model Comparison
comparison_rows = []
for is in 1:n_s
    for m in 1:n_months
        Q_TU = Q_out[:, :, m, is, 1]
        Q_IA = Q_out[:, :, m, is, 2]
        Q_grouped = Q_out[:, :, m, is, 3]

        F_TU = F_out[:, :, m, is, 1]
        F_IA = F_out[:, :, m, is, 2]
        F_grouped = F_out[:, :, m, is, 3]

        push!(comparison_rows, (
            species_key = selected_species[is].species_key,
            species_name = selected_species[is].species_name,
            month = m,
            max_F_t_TU = safe_max(F_TU),
            max_F_t_IA = safe_max(F_IA),
            max_F_t_grouped = safe_max(F_grouped),
            delta_max_F_t_IA_minus_TU = safe_max(abs.(F_IA .- F_TU)),
            delta_max_F_t_grouped_minus_TU = safe_max(abs.(F_grouped .- F_TU)),
            delta_max_F_t_grouped_minus_IA = safe_max(abs.(F_grouped .- F_IA)),
            max_Q_t_TU = safe_max(Q_TU),
            max_Q_t_IA = safe_max(Q_IA),
            max_Q_t_grouped = safe_max(Q_grouped),
            delta_max_Q_t_IA_minus_TU = safe_max(abs.(Q_IA .- Q_TU)),
            delta_max_Q_t_grouped_minus_TU = safe_max(abs.(Q_grouped .- Q_TU)),
            delta_max_Q_t_grouped_minus_IA = safe_max(abs.(Q_grouped .- Q_IA))
        ))
    end
end
CSV.write(joinpath(out_dir, "synthetic_grid_mixture_model_comparison.csv"), DataFrame(comparison_rows))


# -----------------------------------------------------------------------------
# 7. PNG Maps
# -----------------------------------------------------------------------------
println("Generating Maps...")

# Helper to plot 2D grids as heatmaps with CairoMakie
# 1. Compound patterns
fig1 = Figure(size = (1200, 800))
# Plot first 4 compounds, month 8 (summer peak)
m_plot = 8
for i in 1:min(4, n_c)
    row = (i - 1) ÷ 2 + 1
    col = (i - 1) % 2 + 1
    ax = Axis(fig1[row, col], title = "Compound $(i): $(selected_compounds[i].chemical_name) (Month $m_plot)")
    hm = heatmap!(ax, C_t_loaded[:, :, m_plot, i], colormap=:viridis)
    #Colorbar(fig1[row, col+1], hm)
end
save(joinpath(out_dir, "synthetic_grid_compound_patterns.png"), fig1)


# 2. Axis impairment maps (Maintenance, Month 8)
fig2 = Figure(size = (1200, 400))
for im in 1:n_models
    ax = Axis(fig2[1, im], title = "E_maintenance ($(models[im]))")
    hm = heatmap!(ax, E_M_out[:, :, m_plot, im], colormap=:inferno, colorrange=(0, 1))
    #Colorbar(fig2[1, im+1], hm)
end
save(joinpath(out_dir, "synthetic_grid_axis_impairment_maps.png"), fig2)


# 3. Max F_t maps (Species 1, Month 8)
fig3 = Figure(size = (1200, 400))
sp_plot = 1
for im in 1:n_models
    ax = Axis(fig3[1, im], title = "F_t: $(selected_species[sp_plot].species_name) ($(models[im]))")
    hm = heatmap!(ax, F_out[:, :, m_plot, sp_plot, im], colormap=:magma)
    #Colorbar(fig3[1, im+1], hm)
end
save(joinpath(out_dir, "synthetic_grid_max_F_maps.png"), fig3)


# 4. Mixture model deltas
fig4 = Figure(size = (1200, 400))
ax1 = Axis(fig4[1, 1], title = "F_t_IA - F_t_TU")
hm1 = heatmap!(ax1, F_out[:, :, m_plot, sp_plot, 2] .- F_out[:, :, m_plot, sp_plot, 1], colormap=:bwr, colorrange=(-0.1, 0.1))
#Colorbar(fig4[1, 2], hm1)

ax2 = Axis(fig4[1, 3], title = "F_t_grouped - F_t_TU")
hm2 = heatmap!(ax2, F_out[:, :, m_plot, sp_plot, 3] .- F_out[:, :, m_plot, sp_plot, 1], colormap=:bwr, colorrange=(-0.1, 0.1))
#Colorbar(fig4[1, 4], hm2)

ax3 = Axis(fig4[1, 5], title = "F_t_grouped - F_t_IA")
hm3 = heatmap!(ax3, F_out[:, :, m_plot, sp_plot, 3] .- F_out[:, :, m_plot, sp_plot, 2], colormap=:bwr, colorrange=(-0.1, 0.1))
#Colorbar(fig4[1, 6], hm3)

save(joinpath(out_dir, "synthetic_grid_mixture_model_deltas.png"), fig4)

println("Done.")
