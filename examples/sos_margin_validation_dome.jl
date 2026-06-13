# ===========================================================================
# Stress-on-Stress (survival-in-air) margin validation — the DIRECT amplification test.
# Does chronic contaminant burden, routed through the model to an eroded adaptive margin,
# predict reduced resilience to an ACUTE perturbation (emersion/anoxia)? This targets the
# framework's CORE two-timescale claim more directly than Scope for Growth: SoS measures the
# organism's capacity to withstand an acute hit, which is what the margin is *for* (SFG is the
# margin state; SoS is its consequence for surviving an acute event).
#
# Data: ICES DOME 2024 OSPAR CEMP, Mytilus edulis, 17 UK stations (one SURVT station has no
# contaminants), 2012-2022, station-level medians (scripts/extract_dome_sos.jl ->
# sos_dome_ukcemp.csv). Outcome = median survival-in-air (days); higher = more resilient.
# Multi-station, exposure-paired, QA'd, open (CC BY 4.0) -- the data SFG lacked. Real
# gradient: PAH ~100-fold, Pb 15x, PCB 14x across stations.
#
# PRE-REGISTERED (sign frozen before running): the modeled margin A_t correlates POSITIVELY
# with SoS survival (eroded margin -> shorter survival); the MoA-routed margin beats naive
# total load; the result survives body-size control (the Albentosa confound lesson). SoS uses
# standard-size mussels (length 4.4-5.9 cm here) so the size range is naturally narrow.
#
#   julia +release --project=. examples/sos_margin_validation_dome.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const FILE = joinpath(@__DIR__, "..", "data", "external", "sos_dome_ukcemp.csv")

parse_num(x) = (x = strip(x, ['"', ' ']); (x == "" || x == "nd") ? NaN :
                something(tryparse(Float64, x), NaN))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
function spear(x, y)
    m = isfinite.(x) .& isfinite.(y)
    count(m) < 4 && return (NaN, 0)
    return (cor(ordinalrank(x[m]), ordinalrank(y[m])), count(m))
end
sig(r, n) = (n <= 3 ? "" : (t = r * sqrt((n - 2) / (1 - r^2)); abs(t) > 2.7 ? "**" : abs(t) > 2.0 ? "*" : ""))
function partial_spear(x, y, Zcols::Vector{Vector{Float64}})
    m = isfinite.(x) .& isfinite.(y)
    for z in Zcols; m = m .& isfinite.(z); end
    count(m) < length(Zcols) + 3 && return (NaN, count(m))
    cols = Vector{Float64}[ones(count(m))]
    for z in Zcols; push!(cols, ordinalrank(z[m])); end
    Z = reduce(hcat, cols)
    resid(v) = v .- Z * (Z \ v)
    return (cor(resid(ordinalrank(x[m])), resid(ordinalrank(y[m]))), count(m))
end

function main()
# --- read derived station table ---
hdr = String[]; rows = Vector{Vector{String}}()
for line in eachline(FILE)
    startswith(line, "#") && continue
    f = String.(split(line, ","))
    if startswith(line, "station_code,"); hdr = f; continue; end
    length(f) == length(hdr) && push!(rows, f)
end
col(name) = (i = findfirst(==(name), hdr); [parse_num(r[i]) for r in rows])

survival = col("survt_median_d")
LNMEA = col("LNMEA"); DRYWT = col("DRYWT")
metals = ("CD", "CU", "HG", "PB", "ZN")
pah23  = ("NAP", "ACNLE", "ACNE", "FLE", "PA", "ANT")
C = Dict(c => col(c) for c in (metals..., pah23..., "SCB7"))
# total 2-3 ring PAH per station (the Widdows toxic-hydrocarbon axis); nd -> 0 for a sum
totPAH = [sum(isfinite(C[p][i]) ? C[p][i] : 0.0 for p in pah23) for i in eachindex(rows)]
# keep stations with survival + at least the PAH/metal axes
keep = isfinite.(survival) .& isfinite.(totPAH) .& isfinite.(C["CD"])
sites = [rows[i][findfirst(==("station_name"), hdr)] for i in eachindex(rows) if keep[i]]
@printf("stations: %d (Mytilus edulis, UK DOME stress-on-stress, station-level)\n", count(keep))

sub(v) = v[keep]
survival = sub(survival); LNMEA = sub(LNMEA); DRYWT = sub(DRYWT); totPAH = sub(totPAH)
for k in keys(C); C[k] = sub(C[k]); end

# ============================ (1) EMPIRICAL ============================
println("\n=== (1) EMPIRICAL: survival vs contaminant (Spearman; negative = burden lowers survival) ===")
empir = (("PAH(2-3ring)", totPAH), ("PCB(SCB7)", C["SCB7"]),
         ("Cu", C["CU"]), ("Cd", C["CD"]), ("Pb", C["PB"]), ("Hg", C["HG"]), ("Zn", C["ZN"]))
for (lab, v) in empir
    r, n = spear(v, survival); @printf("  survival ~ %-14s rho = %+.3f %-2s (n=%d)\n", lab, r, sig(r, n), n)
end

# ============================ (2) MODEL ============================
# pMoA routing (parallels the SFG harness): assimilation <- 2-3 ring PAH (toxic hydrocarbons);
# maintenance <- metals (Cd,Cu,Hg,Pb,Zn); reproduction <- PCB (no organotins here);
# growth <- baseline (no organotin). Pressure = median-normalised relative burden.
params = amp_species_deb_params(load_amp_species_library(), "Mytilus_edulis")
tu(v) = (med = median(filter(isfinite, v)); [isfinite(x) ? x / med : 1.0 for x in v])
meanrows(vs...) = [mean(t) for t in zip(vs...)]

p_assim  = tu(totPAH)
p_maint  = meanrows((tu(C[m]) for m in metals)...)
p_repro  = tu(C["SCB7"])
p_growth = fill(1.0, length(survival))

margin = Float64[]
for i in eachindex(survival)
    res = compute_adaptive_margin_response((p_assim[i], p_maint[i], p_growth[i], p_repro[i]), params)
    push!(margin, res.A_t)
end

println("\n=== (2) MODEL: modeled adaptive margin vs measured SoS survival (+ = model tracks) ===")
rM, nM = spear(margin, survival)
@printf("  survival ~ margin A_t          rho = %+.3f %-2s (n=%d)  [pre-registered: POSITIVE]\n", rM, sig(rM, nM), nM)
rLN, _ = partial_spear(margin, survival, [LNMEA])
rLD, _ = partial_spear(margin, survival, [LNMEA, DRYWT])
@printf("  survival ~ margin | length     rho = %+.3f %-2s (n=%d)\n", rLN, sig(rLN, count(isfinite.(LNMEA))), count(isfinite.(LNMEA)))
@printf("  survival ~ margin | length,dry rho = %+.3f %-2s\n", rLD, sig(rLD, count(isfinite.(LNMEA) .& isfinite.(DRYWT))))

naive_load = meanrows(tu(totPAH), tu(C["SCB7"]), (tu(C[m]) for m in metals)...)
rPAH, _ = spear(totPAH, survival); rL, _ = spear(naive_load, survival)
println("\n  baselines (|rho|):")
@printf("    survival ~ PAH (best single, validated axis)  |rho| = %.3f\n", abs(rPAH))
@printf("    survival ~ naive mean toxic-unit              |rho| = %.3f\n", abs(rL))
@printf("    MODEL margin                                  |rho| = %.3f\n", abs(rM))

# --------- (2b) axis diagnostic -------
println("\n=== (2b) AXIS DIAGNOSTIC: survival vs each axis pressure (negative = behaves as toxicant) ===")
for (lab, p) in (("p_assim (PAH)", p_assim), ("p_maint (metals)", p_maint), ("p_repro (PCB)", p_repro))
    r, n = spear(p, survival); @printf("  survival ~ %-18s rho = %+.3f %-2s (n=%d)\n", lab, r, sig(r, n), n)
end

println("\nInterpretation: a POSITIVE survival~margin (surviving size control) is the first")
println("direct external support for the framework's amplification thesis -- chronic burden")
println("-> eroded margin -> reduced resilience to an acute perturbation -- on a multi-station,")
println("exposure-paired, QA'd dataset. Watch the axis diagnostic: PAH/assimilation should be")
println("negative (toxic); if metals are positive again (As/Cd/Zn confound, cf. SFG) the MoA")
println("routing earns its keep by keeping them off the toxic axis. Rank stats; * |t|>2.0, ** >2.7.")
println("See docs/notes/ (sos_validation_status.md).")
end

main()
