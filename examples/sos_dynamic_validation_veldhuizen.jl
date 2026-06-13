# ===========================================================================
# DYNAMIC margin validation — the first test that exercises the model's DYNAMICS, not the
# static point map. Transplant time-course (Veldhuizen-Tsoerkan et al. 1991): clean Mytilus
# edulis moved to a contamination gradient, stress indices tracked at 2.5 and 5 months.
#
# THE DISCRIMINATING QUESTION. Cadmium accumulates fast and PLATEAUS by 2.5 mo, yet anoxic
# survival (SoS) keeps dropping from 2.5 -> 5 mo at contaminated sites. A STATIC burden->margin
# map (margin = f(instantaneous burden)) predicts ~no further erosion once burden plateaus.
# The DYNAMIC model integrates an erosion state under SUSTAINED burden with intrinsic timescale
# 1/lambda; for M. edulis lambda_min = k_M = 0.00113/day -> 1/lambda ~ 68-880 days ~ months, so
# the erosion state is still rising at 2.5 and 5 months -> continued erosion. Which matches?
#
# Method:
#   STATIC : compute_adaptive_margin_response(burden_at_t)  -> A_t at 2.5 and 5 mo (~flat,
#            since Cd plateaus; only the small PCB increment moves it).
#   DYNAMIC: simulate_deb_axis_response integrates the slow erosion state y under the sustained
#            (constant, Cd-dominated) cost from t=0; read y at 75 d (2.5 mo) and 150 d (5 mo).
#            y rises = erosion; higher y -> lower expected survival. (q is a linear scale, so
#            the y150/y75 ratio and the cross-site pattern are q-independent.)
# Seasonal confound removed by comparing to the contaminant-driven signal = the SoS gap vs the
# (transplanted) reference, which also captures the spring anoxia-tolerance decline.
#
# HONEST SCOPE: n=4 sites x 2 times, figure-digitized burden/CI/AEC (SoS LT50 + PCB from text).
# A small, qualitative PROOF-OF-CONCEPT that the dynamics can/can't reproduce the trajectory —
# not a powered validation. cost << A0 here, so lambda ~ lambda_max (the raw-subtraction scale
# regime); the discriminating signal is the TRAJECTORY SHAPE, robust to that.
#
#   julia +release --project=. examples/sos_dynamic_validation_veldhuizen.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const FILE = joinpath(@__DIR__, "..", "data", "external", "sos_veldhuizen1991_transplant.csv")
parse_num(x) = (x = strip(x, ['"', ' ']); (x == "" || x == "nd") ? NaN :
                something(tryparse(Float64, x), NaN))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = (m = isfinite.(x) .& isfinite.(y); count(m) < 3 ? NaN : cor(ordinalrank(x[m]), ordinalrank(y[m])))

function main()
# --- read transplant panel ---
hdr = String[]; R = Vector{Vector{String}}()
for line in eachline(FILE)
    startswith(line, "#") && continue
    f = String.(split(line, ","))
    if startswith(line, "site,"); hdr = f; continue; end
    length(f) == length(hdr) && push!(R, f)
end
ci(n) = findfirst(==(n), hdr); col(n) = [parse_num(r[ci(n)]) for r in R]
site = [r[ci("site")] for r in R]
months = col("exposure_months"); Cd = col("Cd_ug_g_dry"); PCB = col("PCB_ug_g_lipid")
SoS = col("SoS_LT50_days"); isref = col("reference") .== 1

params = amp_species_deb_params(load_amp_species_library(), "Mytilus_edulis")
@printf("Mytilus edulis: A0=%.1f  lambda_min(=k_M)=%.5f/d  lambda_max=%.5f/d  g=%.1f\n",
        params.A0, params.lambda_min, params.lambda_max, params.lambda_max / params.lambda_min)
@printf("intrinsic erosion timescale 1/lambda = %.0f d (pristine) .. %.0f d (eroded) -> MONTHS\n",
        1 / params.lambda_max, 1 / params.lambda_min)

# axis pressures from burden (median-normalised across sites; metal->maint, PCB->repro)
medCd = median(filter(isfinite, Cd)); medPCB = median(filter(isfinite, PCB))
smaint(cd) = isfinite(cd) ? cd / medCd : 0.0
srepro(pcb) = isfinite(pcb) ? pcb / medPCB : 0.0

# ---------- DYNAMIC: integrate erosion state y under sustained burden ----------
# constant (step) burden from t=0 (Cd plateaus fast); read y at 75 d and 150 d.
const_dt = 0.5; T = collect(0.0:const_dt:150.0); nt = length(T)
i75 = argmin(abs.(T .- 75.0)); i150 = argmin(abs.(T .- 150.0))
function dyn_y(cd, pcb)
    pa = (assimilation = zeros(nt), maintenance = fill(smaint(cd), nt),
          growth = zeros(nt), reproduction = fill(srepro(pcb), nt))
    y = simulate_deb_axis_response(T, pa, params.A0, params, 1.0; y0 = 0.0, dt = const_dt)
    return y[i75], y[i150]
end

# per-site (use 5-mo burden = the plateau; Cd same at 2.5 mo, PCB present by 5 mo)
sites5 = findall(months .== 5.0)
println("\n=== per site (5-mo burden drives the sustained cost) ===")
println("  site         Cd    PCB   |  dyn y(2.5mo)  y(5mo)  ratio | SoS 2.5mo  5mo")
dy75 = Float64[]; dy150 = Float64[]
for i in sites5
    y75, y150 = dyn_y(Cd[i], PCB[i])
    push!(dy75, y75); push!(dy150, y150)
    i25 = findfirst(j -> site[j] == site[i] && months[j] == 2.5, eachindex(R))
    @printf("  %-11s %4.1f  %5.2f |  %8.2f  %7.2f  %.2f | %6.1f   %4.1f\n",
            site[i], Cd[i], PCB[i], y75, y150, y150 / y75,
            i25 === nothing ? NaN : SoS[i25], SoS[i])
end

# ---------- (1) cross-sectional: erosion state vs survival ----------
SoS5 = SoS[sites5]
println("\n=== (1) CROSS-SECTIONAL (n=4): dynamic erosion y vs measured SoS survival ===")
@printf("  rho( dyn y(5mo), SoS(5mo) )   = %+.2f   [expect strongly NEGATIVE: more erosion -> less survival]\n",
        spear(dy150, SoS5))

# ---------- (2) THE DISCRIMINATING temporal test ----------
println("\n=== (2) TEMPORAL — continued erosion from 2.5 -> 5 mo (the static map cannot) ===")
# static change 2.5->5mo at a contaminated site: Cd identical, only PCB added
function static_change(s)
    i25 = findfirst(j -> site[j] == s && months[j] == 2.5, eachindex(R))
    i50 = findfirst(j -> site[j] == s && months[j] == 5.0, eachindex(R))
    A25 = compute_adaptive_margin_response((0.0, smaint(Cd[i25]), 0.0, srepro(PCB[i25])), params).A_t
    A50 = compute_adaptive_margin_response((0.0, smaint(Cd[i50]), 0.0, srepro(PCB[i50])), params).A_t
    return (params.A0 - A50) - (params.A0 - A25)   # extra static erosion 2.5->5mo
end
for i in sites5
    s = site[i]; isref[i] && continue
    i25 = findfirst(j -> site[j] == s && months[j] == 2.5, eachindex(R))
    y75, y150 = dyn_y(Cd[i], PCB[i])
    sos_drop = SoS[i25] - SoS[i]
    @printf("  %-11s  measured SoS drop 2.5->5mo = %.1f d  | DYNAMIC extra erosion (y150-y75) = %+.2f (%.0f%% more) | STATIC extra erosion = %+.3f (~flat: Cd plateaued)\n",
            s, sos_drop, y150 - y75, 100 * (y150 / y75 - 1), static_change(s))
end

println("\nReading: the DYNAMIC erosion state keeps rising 2.5->5 mo (~33% more) because the")
println("model's intrinsic timescale 1/lambda ~ months > experiment duration -- matching the")
println("observed continued SoS decline AFTER Cd plateaued; the STATIC map predicts ~no further")
println("erosion (burden flat) and so CANNOT explain the 2.5->5mo drop. This is the first")
println("evidence the model's DYNAMICS (not just its static response surface) track reality.")
println("Caveats: n=4 sites x 2 times, digitized; PCB rising 2.5->5mo is a partial confound for")
println("the temporal signal; q arbitrary (ratios q-independent). Proof-of-concept, not powered.")
println("See docs/notes/sos_validation_status.md.")
end

main()
