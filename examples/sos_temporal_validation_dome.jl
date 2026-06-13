# ===========================================================================
# Stress-on-Stress TEMPORAL margin validation (power-boost + erosion-over-time).
# Companion to the station-level test (examples/sos_margin_validation_dome.jl, n=17, +0.39).
# Uses the station-YEAR panel (scripts/extract_dome_sos_yearly.jl) to ask the dynamic
# question the station-level test cannot: WITHIN a station, do years with higher modeled
# margin (lower contaminant burden) have higher survival-in-air? That is the erosion
# mechanism tested OVER TIME, and a within-station fixed-effects design removes ALL fixed
# between-station differences (the cleanest possible confound control).
#
# Three readouts:
#   (A) WITHIN-STATION fixed effects (THE temporal test): de-mean margin & survival by
#       station, correlate the year-to-year deviations. Removes every station-level confound.
#   (B) Pooled station-years (more obs, but pseudoreplicated between-station => significance
#       is OVERstated; reported for context only).
#   (C) Station-level recap (the honest primary unit, n=17).
#
# PRE-REGISTERED sign: positive (higher margin -> higher survival), within and pooled.
#   julia +release --project=. examples/sos_temporal_validation_dome.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const FILE = joinpath(@__DIR__, "..", "data", "external", "sos_dome_ukcemp_yearly.csv")

parse_num(x) = (x = strip(x, ['"', ' ']); (x == "" || x == "nd") ? NaN :
                something(tryparse(Float64, x), NaN))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
function spear(x, y)
    m = isfinite.(x) .& isfinite.(y); count(m) < 4 && return (NaN, 0)
    return (cor(ordinalrank(x[m]), ordinalrank(y[m])), count(m))
end
sig(r, n) = (n <= 3 ? "" : (t = r * sqrt((n - 2) / (1 - r^2)); abs(t) > 2.7 ? "**" : abs(t) > 2.0 ? "*" : ""))

function main()
hdr = String[]; R = Vector{Vector{String}}()
for line in eachline(FILE)
    startswith(line, "#") && continue
    f = String.(split(line, ","))
    if startswith(line, "station_code,"); hdr = f; continue; end
    length(f) == length(hdr) && push!(R, f)
end
ci(name) = findfirst(==(name), hdr)
getc(name) = [parse_num(r[ci(name)]) for r in R]
station = [r[ci("station_code")] for r in R]
survival = getc("survt_median_d"); LNMEA = getc("LNMEA"); DRYWT = getc("DRYWT")
metals = ("CD","CU","HG","PB","ZN"); pah23 = ("NAP","ACNLE","ACNE","FLE","PA","ANT")
C = Dict(c => getc(c) for c in (metals..., pah23..., "SCB7"))
totPAH = [sum(isfinite(C[p][i]) ? C[p][i] : 0.0 for p in pah23) for i in eachindex(R)]

# QC: keep station-years with a real PAH+metal match; drop the Cu>50 mg/kg outlier (7903/2013)
keep = isfinite.(survival) .& (totPAH .> 0) .& isfinite.(C["CD"]) .& (C["CU"] .< 50000)
idx = findall(keep)
@printf("usable station-years: %d (dropped %d: no-match / QC outlier)\n", length(idx), length(R) - length(idx))

sub(v) = v[idx]
station = station[idx]; survival = sub(survival); LNMEA = sub(LNMEA); DRYWT = sub(DRYWT)
totPAH = sub(totPAH); for k in keys(C); C[k] = sub(C[k]); end

# margin per station-year (same MoA routing as the station-level test)
params = amp_species_deb_params(load_amp_species_library(), "Mytilus_edulis")
tu(v) = (med = median(filter(isfinite, v)); [isfinite(x) ? x / med : 1.0 for x in v])
meanrows(vs...) = [mean(t) for t in zip(vs...)]
p_assim = tu(totPAH); p_maint = meanrows((tu(C[m]) for m in metals)...)
p_repro = tu(C["SCB7"]); p_growth = fill(1.0, length(survival))
margin = [compute_adaptive_margin_response((p_assim[i], p_maint[i], p_growth[i], p_repro[i]), params).A_t
          for i in eachindex(survival)]

# ---------- (A) WITHIN-STATION fixed effects (the temporal test) ----------
# de-mean margin & survival within each station; correlate year-to-year deviations.
stations = unique(station)
multi = [s for s in stations if count(==(s), station) >= 2]
dm = Float64[]; ds = Float64[]
for s in multi
    rows = findall(==(s), station)
    mbar = mean(margin[rows]); sbar = mean(survival[rows])
    append!(dm, margin[rows] .- mbar); append!(ds, survival[rows] .- sbar)
end
nstat = length(multi); nobs = length(dm); dfw = nobs - nstat
rW = cor(ordinalrank(dm), ordinalrank(ds))
tW = rW * sqrt(dfw / (1 - rW^2))
println("\n=== (A) WITHIN-STATION fixed effects — erosion OVER TIME (the clean temporal test) ===")
@printf("  %d multi-year stations, %d station-years, df=%d\n", nstat, nobs, dfw)
@printf("  within-station  survival ~ margin   rho = %+.3f  (t=%.2f, %s)\n",
        rW, tW, abs(tW) > 2.7 ? "p<0.01" : abs(tW) > 2.0 ? "p<0.05" : "n.s.")
@printf("  [+ = years with higher modeled margin have higher acute-stress survival, within station]\n")

# ---------- (B) pooled station-years (context; pseudoreplicated) ----------
println("\n=== (B) pooled station-years (more obs, but between-station pseudoreplicated) ===")
rP, nP = spear(margin, survival)
@printf("  survival ~ margin   rho = %+.3f %-2s (n=%d)   [significance OVERstated — see (A)/(C)]\n", rP, sig(rP, nP), nP)
for (lab, p) in (("p_assim (PAH)", p_assim), ("p_maint (metals)", p_maint), ("p_repro (PCB)", p_repro))
    r, n = spear(p, survival); @printf("    survival ~ %-18s rho = %+.3f %-2s (n=%d)\n", lab, r, sig(r, n), n)
end

# ---------- (C) station-level recap (the honest primary unit) ----------
println("\n=== (C) station-level recap (one obs per station = the honest primary unit) ===")
smarg = Float64[]; ssurv = Float64[]
for s in stations
    rows = findall(==(s), station)
    push!(smarg, mean(margin[rows])); push!(ssurv, mean(survival[rows]))
end
rS, nS = spear(smarg, ssurv)
@printf("  survival ~ margin   rho = %+.3f %-2s (n=%d)\n", rS, sig(rS, nS), nS)

println("\nInterpretation: (A) is the key new result — a within-station fixed-effects design")
println("removes every fixed between-station confound (size regime, population, local")
println("hydrography), so a positive within-station margin~survival is the cleanest evidence")
println("that the SAME mussels' acute-stress resilience erodes as modeled margin falls year to")
println("year. Rank stats; pre-registered positive. Caveat: thin panel (mostly 2 yr/station).")
println("See docs/notes/sos_validation_status.md.")
end

main()
