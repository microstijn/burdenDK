# ===========================================================================
# Response-curve-form sensitivity — ranking-stability robustness check.
#
# The licensed use of the framework is RELATIVE (mechanistically-structured rankings
# of sites/scenarios), so the defensibility question is: does the corroboration AND
# the site ranking survive a change in the per-axis impairment curve? We already
# sensitivity-test the MoA routing (Discussion); this adds the response-curve form.
#
# Method (exact, reuses the real engine). The engine applies E(x)=x/(1+x) per axis.
# To reproduce an ALTERNATIVE impairment E_alt through the same engine (same AmP
# weights w_i, same A0, same lambda), feed a transformed pressure
#     x' = E_alt(x) / (1 - E_alt(x)).
# Then E_engine(x') = x'/(1+x') = E_alt(x) exactly. So:
#   - Hill exponent h:  E_alt = x^h/(1+x^h)  =>  x' = x^h     (h=1 is baseline)
#   - saturating-exp:   E_alt = 1 - exp(-x)  =>  x' = exp(x)-1
# Every curve here is monotone, bounded in [0,1), threshold-free, E(0)=0, and
# HALF-SATURATES AT THE REFERENCE x=1 by construction -- a SHAPE probe, NOT a
# reintroduced half-saturation location knob (cf. the removed K_A). No tuning knob.
#
# Anchors: Widdows 1995 SFG (n=36) and ICES DOME stress-on-stress (n=17). We report
#   (a) margin<->outcome Spearman rho under each curve (does corroboration hold?), and
#   (b) DOME ranking stability: Spearman between the baseline station ranking and each
#       alternative's ranking (does the licensed relative ordering hold?).
#
#   julia +release --project=. examples/response_curve_sensitivity.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const EXT = joinpath(@__DIR__, "..", "data", "external")

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
function spear(x, y)
    m = isfinite.(x) .& isfinite.(y)
    count(m) < 4 && return (NaN, 0)
    return (cor(ordinalrank(x[m]), ordinalrank(y[m])), count(m))
end
sig(r, n) = (n <= 3 ? "" : (t = r * sqrt((n - 2) / (1 - r^2)); abs(t) > 2.7 ? "**" : abs(t) > 2.0 ? "*" : ""))

# --- the response-curve family (transform on per-axis pressure x) ---
const CURVES = [
    ("x/(1+x)  [baseline]",   x -> x),
    ("Hill h=2  (steeper)",   x -> x^2),
    ("Hill h=0.5 (gentler)",  x -> sqrt(max(x, 0.0))),
    ("1-exp(-x) (sat-exp)",   x -> exp(x) - 1),
]

# margin per row under a given curve transform T
function margins(pressures::NTuple{4,Vector{Float64}}, T, params)
    p = map(v -> T.(v), pressures)
    [compute_adaptive_margin_response((p[1][i], p[2][i], p[3][i], p[4][i]), params).A_t
     for i in eachindex(p[1])]
end

tu(v) = (med = median(filter(isfinite, v)); [isfinite(x) ? x / med : 1.0 for x in v])
meanrows(vs...) = [mean(t) for t in zip(vs...)]

# =========================== DOME stress-on-stress ===========================
function load_dome()
    parse_num(x) = (x = strip(x, ['"', ' ']); (x == "" || x == "nd") ? NaN :
                    something(tryparse(Float64, x), NaN))
    hdr = String[]; rows = Vector{Vector{String}}()
    for line in eachline(joinpath(EXT, "sos_dome_ukcemp.csv"))
        startswith(line, "#") && continue
        f = String.(split(line, ","))
        if startswith(line, "station_code,"); hdr = f; continue; end
        length(f) == length(hdr) && push!(rows, f)
    end
    col(name) = (i = findfirst(==(name), hdr); [parse_num(r[i]) for r in rows])
    survival = col("survt_median_d")
    metals = ("CD", "CU", "HG", "PB", "ZN"); pah23 = ("NAP","ACNLE","ACNE","FLE","PA","ANT")
    C = Dict(c => col(c) for c in (metals..., pah23..., "SCB7"))
    totPAH = [sum(isfinite(C[p][i]) ? C[p][i] : 0.0 for p in pah23) for i in eachindex(rows)]
    keep = isfinite.(survival) .& isfinite.(totPAH) .& isfinite.(C["CD"])
    sub(v) = v[keep]
    survival = sub(survival); totPAH = sub(totPAH); for k in keys(C); C[k] = sub(C[k]); end
    p_assim = tu(totPAH); p_maint = meanrows((tu(C[m]) for m in metals)...)
    p_repro = tu(C["SCB7"]); p_growth = fill(1.0, length(survival))
    return (p_assim, p_maint, p_growth, p_repro), survival
end

# =========================== Widdows 1995 SFG ===========================
function load_sfg()
    norm_site(s) = lowercase(replace(strip(s), r"[^A-Za-z0-9]" => ""))
    parse_conc(x) = (x = strip(x); (x == "" || x == "nd") ? NaN :
                     (startswith(x, "<") ? parse(Float64, x[2:end]) : parse(Float64, x)))
    sfg = Dict{String,Float64}()
    for line in eachline(joinpath(EXT, "sfg_widdows1995_northsea.csv"))
        (startswith(line, "#") || startswith(line, "site,")) && continue
        f = split(line, ","); length(f) >= 6 || continue
        sfg[norm_site(f[1])] = parse(Float64, f[6])
    end
    chead = String[]; cont = Dict{String,Dict{String,Float64}}()
    for line in eachline(joinpath(EXT, "sfg_widdows1995_contaminants.csv"))
        startswith(line, "#") && continue
        f = split(line, ",")
        if startswith(line, "site,"); chead = String.(f); continue; end
        length(f) == length(chead) || continue
        d = Dict{String,Float64}(); for j in 2:length(f); d[chead[j]] = parse_conc(f[j]); end
        cont[norm_site(f[1])] = d
    end
    findsite(k, ks) = (k in ks && return k; for kk in ks; (startswith(kk, k) || startswith(k, kk)) && return kk; end; nothing)
    sk = collect(keys(sfg)); SFG = Float64[]; C = Dict{String,Vector{Float64}}()
    for ck in keys(cont)
        s = findsite(ck, sk); s === nothing && continue
        push!(SFG, sfg[s]); for (name, v) in cont[ck]; push!(get!(C, name, Float64[]), v); end
    end
    p_assim = tu(C["total_toxic_HC"]); p_maint = meanrows(tu(C["Cd"]),tu(C["Cu"]),tu(C["Hg"]),tu(C["Pb"]),tu(C["Zn"]))
    p_repro = meanrows(tu(C["TBT"]), tu(C["PCB"])); p_growth = tu(C["DBT"])
    return (p_assim, p_maint, p_growth, p_repro), SFG
end

function report(name, pressures, outcome, params)
    println("\n=== $name (n=$(length(outcome))) ===")
    @printf("  %-22s  %10s   %s\n", "response curve", "rho(margin,outcome)", "ranking-stability vs baseline")
    base = margins(pressures, CURVES[1][2], params)
    for (lab, T) in CURVES
        m = margins(pressures, T, params)
        r, n = spear(m, outcome)
        rstab, _ = spear(base, m)                      # Spearman between rankings
        stab = lab == CURVES[1][1] ? "  (—)" : @sprintf("  rho_rank = %+.3f", rstab)
        @printf("  %-22s  %+0.3f %-2s          %s\n", lab, r, sig(r, n), stab)
    end
end

function main()
    params = amp_species_deb_params(load_amp_species_library(), "Mytilus_edulis")
    dome_p, dome_y = load_dome()
    sfg_p,  sfg_y  = load_sfg()
    println("Response-curve-form sensitivity (ranking stability). All curves are threshold-free,")
    println("bounded, monotone, and half-saturate at the reference x=1 (shape probes, no knob).")
    report("Widdows 1995 SFG (margin <-> SFG)", sfg_p, sfg_y, params)
    report("ICES DOME stress-on-stress (margin <-> survival-in-air)", dome_p, dome_y, params)
    println("\nReading: if rho(margin,outcome) stays the same sign/strength and ranking-stability")
    println("(rho_rank) is near +1 across curves, then neither the corroboration nor the licensed")
    println("relative ranking depends on the specific impairment form E=x/(1+x). * |t|>2.0, ** >2.7.")
end

main()
