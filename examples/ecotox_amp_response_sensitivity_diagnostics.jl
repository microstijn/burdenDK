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

output_dir = joinpath(@__DIR__, "..", "output", "ecotox_amp_response_sensitivity_diagnostics")
mkpath(output_dir)

amp = load_amp_species_library()
ecotox = load_ecotox_library()
memory = load_compound_memory_library()

species_name = "Abatus cordatus"
base_params = amp_species_deb_params(amp, species_name)

cas_target = "7440-43-9" # Cadmium

# --- Select Compound Record ---
cas_target = "7440-43-9" # Cadmium

# Filter ecotox library for cadmium records
cadmium_records = filter(r -> r["cas_norm"] == TwoTimescaleResilience._normalize_cas_runtime(cas_target) || r["cas"] == cas_target, ecotox)

# Priority order: MOR, GRO, REP, then first valid
selected_record = nothing

# Wrap validate_ecotox_record in try/catch to safely ignore invalid records
function safe_validate(r)
    try
        return TwoTimescaleResilience.validate_ecotox_record(r)
    catch
        return false
    end
end

for code in ["MOR", "GRO", "REP"]
    recs = filter(r -> r["effect_code"] == code && safe_validate(r), cadmium_records)
    if !isempty(recs)
        global selected_record = recs[1]
        break
    end
end

if selected_record === nothing
    recs = filter(safe_validate, cadmium_records)
    if !isempty(recs)
        global selected_record = recs[1]
    else
        error("No valid ECOTOX record found for Cadmium (CAS: $(cas_target))")
    end
end

println("=== Selected ECOTOX Record ===")
println("CAS Hyphenated:  ", selected_record["cas_hyphenated"])
println("CAS Normalized:  ", selected_record["cas_norm"])
println("Taxon Class:     ", selected_record["taxon_class"])
println("Effect Code:     ", selected_record["effect_code"])
println("Routed Axis:     ", ecotox_effect_to_deb_axis(selected_record["effect_code"]))
println("NOEC_median:     ", selected_record["NOEC_median"])
println("EC50_median:     ", selected_record["EC50_median"])
println("n_NOEC:          ", selected_record["n_NOEC"])
println("n_EC50:          ", selected_record["n_EC50"])
println()

# --- Memory Parameters ---
rho = compound_retention(cas_target; memory_library=memory)
K = compound_bioaccumulation_factor(cas_target; memory_library=memory)

println("=== Compound Memory Parameters ===")
println("Retention (rho): ", rho)
println("Bioaccumulation Factor (K): ", K)
println()

# --- Monthly Exposure Series ---
EC50 = selected_record["EC50_median"]
NOEC = selected_record["NOEC_median"]

pulse = 10.0 * EC50
if !isfinite(pulse)
    error("Pulse value is not finite.")
end

C_t = [
    0.0, 0.0,
    pulse, pulse, pulse,
    0.0, 0.0, 0.0,
    pulse, pulse,
    0.0, 0.0
]

# The pulse is deliberately high to diagnose the response curve.
# This is a pedagogical stress scenario, not a claim about realistic ambient concentrations.

# --- Parameter Variants ---
function copy_deb_params_with(
    p;
    A0 = p.A0,
    alpha_axes = p.alpha_axes,
    lambda_min = p.lambda_min,
    lambda_max = p.lambda_max,
    KA = p.KA,
    recovery_axes = p.recovery_axes,
    use_axis_recovery_penalty = p.use_axis_recovery_penalty,
    use_buffer_recovery_factor = p.use_buffer_recovery_factor,
    beta_Z = p.beta_Z
)
    return DEBAxisParams(
        A0 = A0,
        alpha_axes = alpha_axes,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        KA = KA,
        recovery_axes = recovery_axes,
        use_axis_recovery_penalty = use_axis_recovery_penalty,
        use_buffer_recovery_factor = use_buffer_recovery_factor,
        beta_Z = beta_Z
    )
end

variants = Dict(
    "baseline" => base_params,
    "low_A0" => copy_deb_params_with(base_params, A0 = 100.0),
    "very_low_A0" => copy_deb_params_with(base_params, A0 = 50.0),
    "high_alpha" => copy_deb_params_with(base_params, alpha_axes = Tuple(x * 10.0 for x in base_params.alpha_axes)),
    "very_high_alpha" => copy_deb_params_with(base_params, alpha_axes = Tuple(x * 50.0 for x in base_params.alpha_axes)),
    "low_KA" => copy_deb_params_with(base_params, KA = min(base_params.KA, 50.0)),
    "very_low_KA" => copy_deb_params_with(base_params, KA = 10.0),
    "combined_sensitive" => copy_deb_params_with(
        base_params,
        A0 = 100.0,
        alpha_axes = Tuple(x * 10.0 for x in base_params.alpha_axes),
        KA = 50.0
    )
)

# Filter valid variants
valid_variants = Dict{String, DEBAxisParams}()
for (k, v) in variants
    if v.A0 > 0 && v.KA > 0 && v.lambda_max > v.lambda_min && all(isfinite, v.alpha_axes) && all(>=(0), v.alpha_axes)
        valid_variants[k] = v
    end
end

variant_order = [
    "baseline", "low_A0", "very_low_A0", "high_alpha", 
    "very_high_alpha", "low_KA", "very_low_KA", "combined_sensitive"
]
valid_ordered_names = filter(name -> haskey(valid_variants, name), variant_order)

println("=== Parameter Variants ===")
@printf("%-20s | %-6s | %-35s | %-6s | %-10s | %-10s | %-10s\n", "variant", "A0", "alpha_axes", "KA", "lambda_min", "lambda_max", "lambda0")
for name in valid_ordered_names
    p = valid_variants[name]
    l0 = restoring_force_from_margin(p.A0, p)
    alpha_str = string(round.(p.alpha_axes, digits=3))
    @printf("%-20s | %-6.1f | %-35s | %-6.1f | %-10.4f | %-10.4f | %-10.4f\n", name, p.A0, alpha_str, p.KA, p.lambda_min, p.lambda_max, l0)
end
println()

# --- Run Monthly Stateful Response for Each Variant ---

results = Dict{String, Any}()

println("=== Diagnostic Summaries ===")
@printf("%-20s | %-6s | %-6s | %-25s | %-8s | %-8s | %-10s | %-8s | %-8s | %-8s\n", "variant", "A0", "KA", "alpha_scale_description", "max_x", "min_A", "min_lambda", "max_F", "final_B", "final_F")

for name in valid_ordered_names
    p = valid_variants[name]
    
    state = EcotoxExposureState()
    
    track_month = Int[]
    track_C = Float64[]
    track_B = Float64[]
    track_x = Float64[]
    track_A = Float64[]
    track_lambda = Float64[]
    track_F = Float64[]
    
    for m in 1:length(C_t)
        C = C_t[m]
        concs = Dict(cas_target => C)
        burden = ecotox_records_to_deb_burden_stateful!(state, concs, [selected_record]; memory_library = memory)
        
        B = get_internal_burden(state, cas_target)
        x_val = ecotox_active_stress(B, NOEC, EC50)
        
        response = ecotox_burden_to_response(burden, p)
        
        push!(track_month, m)
        push!(track_C, C)
        push!(track_B, B)
        push!(track_x, x_val)
        push!(track_A, response.A)
        push!(track_lambda, response.lambda)
        push!(track_F, response.amplification)
    end
    
    max_B = maximum(track_B)
    max_x = maximum(track_x)
    min_A = minimum(track_A)
    min_lambda = minimum(track_lambda)
    max_F = maximum(track_F)
    final_B = track_B[end]
    final_F = track_F[end]
    
    results[name] = (
        month = track_month,
        C = track_C,
        B = track_B,
        x = track_x,
        A = track_A,
        lambda = track_lambda,
        F = track_F
    )
    
    alpha_desc = (name == "high_alpha" || name == "combined_sensitive") ? "10x baseline" : (name == "very_high_alpha" ? "50x baseline" : "baseline")
    
    @printf("%-20s | %-6.1f | %-6.1f | %-25s | %-8.2f | %-8.2f | %-10.4f | %-8.2f | %-8.2f | %-8.2f\n", name, p.A0, p.KA, alpha_desc, max_x, min_A, min_lambda, max_F, final_B, final_F)
end
println()

println("=== Amplification Bounds ===")
@printf("%-20s | %-10s | %-10s | %-20s\n", "variant", "lambda0", "lambda_min", "theoretical_F_if_A_zero")
for name in valid_ordered_names
    p = valid_variants[name]
    l0 = restoring_force_from_margin(p.A0, p)
    l_min_calc = restoring_force_from_margin(0.0, p)
    th_F = l0 / l_min_calc
    @printf("%-20s | %-10.4f | %-10.4f | %-20.2f\n", name, l0, p.lambda_min, th_F)
end
println()

# --- Plot 1: lambda(A) Response Curves ---
fig1 = Figure(size = (800, 600))
ax1 = Axis(fig1[1, 1], xlabel = "A / A0 (Normalized Adaptive Margin)", ylabel = "Restoring Force lambda(A)", title = "Response Curves by Variant")

colors = Makie.wong_colors()
for (i, name) in enumerate(valid_ordered_names)
    p = valid_variants[name]
    A_vals = range(0.0, p.A0, length=300)
    A_norm = A_vals ./ p.A0
    l_vals = [restoring_force_from_margin(A, p) for A in A_vals]
    lines!(ax1, A_norm, l_vals, label=name, color=colors[mod1(i, length(colors))], linewidth=2)
end
axislegend(ax1, position=:rb)
save(joinpath(output_dir, "lambda_vs_normalized_A.png"), fig1)

# --- Plot 2: Monthly Amplification by Variant ---
fig2 = Figure(size = (800, 600))
ax2 = Axis(fig2[1, 1], xlabel = "Month", ylabel = "Amplification (F)", title = "Monthly Amplification by Variant")

for (i, name) in enumerate(valid_ordered_names)
    F_vals = results[name].F
    months = results[name].month
    lines!(ax2, months, F_vals, label=name, color=colors[mod1(i, length(colors))], linewidth=2)
    scatter!(ax2, months, F_vals, color=colors[mod1(i, length(colors))], markersize=8)
end
axislegend(ax2, position=:lt)
save(joinpath(output_dir, "monthly_amplification_by_variant.png"), fig2)

# --- Plot 3: Monthly Internal Burden and Active Stress ---
fig3 = Figure(size = (800, 800))
ax3a = Axis(fig3[1, 1], xlabel = "Month", ylabel = "Concentration", title = "Baseline Stateful Run: C_t and B_t")
ax3b = Axis(fig3[2, 1], xlabel = "Month", ylabel = "Active Stress (x_stateful)", title = "Baseline Stateful Run: x_stateful")

baseline_res = results["baseline"]
months = baseline_res.month

lines!(ax3a, months, baseline_res.C, label="Ambient C", color=:blue, linewidth=2)
scatter!(ax3a, months, baseline_res.C, color=:blue, markersize=8)
lines!(ax3a, months, baseline_res.B, label="Internal Burden B", color=:red, linewidth=2)
scatter!(ax3a, months, baseline_res.B, color=:red, markersize=8)
axislegend(ax3a, position=:lt)

lines!(ax3b, months, baseline_res.x, color=:green, linewidth=2)
scatter!(ax3b, months, baseline_res.x, color=:green, markersize=8)

save(joinpath(output_dir, "monthly_B_and_x.png"), fig3)

# --- Plot 4: Monthly A by Variant ---
fig4 = Figure(size = (800, 600))
ax4 = Axis(fig4[1, 1], xlabel = "Month", ylabel = "Adaptive Margin (A)", title = "Monthly A by Variant")

for (i, name) in enumerate(valid_ordered_names)
    A_vals = results[name].A
    months = results[name].month
    lines!(ax4, months, A_vals, label=name, color=colors[mod1(i, length(colors))], linewidth=2)
    scatter!(ax4, months, A_vals, color=colors[mod1(i, length(colors))], markersize=8)
end
axislegend(ax4, position=:lb)
save(joinpath(output_dir, "monthly_A_by_variant.png"), fig4)

# --- Interpretation Block ---
println("=== Interpretation & Findings ===")
println("1. Chemical memory can produce sustained internal burden B after ambient concentration C returns to zero.")
println("2. Bioaccumulation factor K magnifies internal burden relative to ambient concentration.")
println("3. Large active stress x does not automatically imply large amplification F.")
println("4. F depends on how much A is depleted relative to A0 and how sensitive lambda(A) is over that range.")
println("5. If A remains large compared with KA, lambda(A) can remain near lambda_max and amplification stays near 1.")
println("6. Parameter variants with lower A0, higher alpha_axes, or smaller KA are diagnostic tools to reveal sensitivity. They are not claims that the AmP-derived baseline parameters are wrong.")
println("7. This script is for response sensitivity diagnostics, not calibration.")
println()

# --- Sanity Checks ---
println("=== Sanity Checks ===")

baseline_B = results["baseline"].B
if baseline_B[5] <= 0
    @warn "Baseline B should be greater than zero after pulse months."
end
if baseline_B[6] <= 0
    @warn "Baseline B should remain greater than zero immediately after pulse ends."
end

baseline_max_x = results["baseline"].x |> maximum
if baseline_max_x <= 1.0
    @warn "max_x_stateful for baseline should be > 1 under the high pulse scenario. Found: $baseline_max_x"
end

baseline_max_F = results["baseline"].F |> maximum
any_sensitive_greater = false
for name in valid_ordered_names
    if name != "baseline" && maximum(results[name].F) > baseline_max_F
        global any_sensitive_greater = true
        break
    end
end

if !any_sensitive_greater
    @warn "At least one sensitive variant should have max_F greater than baseline max_F."
end

for name in valid_ordered_names
    p = valid_variants[name]
    l0 = restoring_force_from_margin(p.A0, p)
    l_min = restoring_force_from_margin(0.0, p)
    if l0 < l_min
        @warn "lambda(A0) should be >= lambda(0.0) for variant $name. l0=$l0, l_min=$l_min"
    end
end

println("Response sensitivity diagnostics completed.")
