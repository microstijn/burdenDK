# ===========================================================================
# GlobTherm coherence / general-resilience probe (external validation #5).
#
# PRE-REGISTERED (before computing): the adaptive-margin model does NOT predict
# thermal limits directly, so this is NOT a test of the margin's erosion mechanism
# (that is the Scope-for-Growth follow-up). It is a COHERENCE probe of one
# question, with a directional prior:
#
#   H (general resilience): species with greater recovery CAPACITY (k_M, the
#   COMADRE-validated recovery floor; lambda_max, the fast bound) tolerate BROADER
#   thermal ranges (CTmax - CTmin), beyond body size and latitude (Janzen) confounds.
#
#   Null is informative: no signal => recovery capacity is stressor-SPECIFIC
#   (thermal tolerance is a separate physiological axis), which honestly bounds the
#   model's reach. A positive => the capacity axis generalises to thermal resilience.
#
# Circularity guard: AmP rates are temperature-corrected and carry the species'
# thermal regime, which also drives CTmax. So (a) control body size + |latitude|,
# and (b) re-run on ECTOTHERMS only (Class not Mammalia/Aves) where T_typical is
# environmental, additionally controlling T_typical. Endotherm T_typical is the AmP
# default body temperature (~310.65 K) and cannot de-confound regime.
#
# READ-ONLY. Input: data/external/globtherm_amp_matched.csv
#   julia +release --project=. examples/globtherm_validation.jl
# ===========================================================================

using Statistics, Printf

const TBL = joinpath(@__DIR__, "..", "data", "external", "globtherm_amp_matched.csv")

pf(x) = (x == "" ? NaN : parse(Float64, x))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = cor(ordinalrank(x), ordinalrank(y))
sig(r, df) = (df <= 1 ? "" : (t = r * sqrt(df / (1 - r^2)); abs(t) > 2.6 ? "**" : abs(t) > 1.97 ? "*" : ""))
# partial rank correlation of x,y controlling covariates (ranked)
function prank(x, y, ctrls::Vector)
    n = length(x); Z = hcat(ones(n), (ordinalrank(c) for c in ctrls)...)
    res(v) = (vr = ordinalrank(v); vr .- Z * (Z \ vr)); cor(res(x), res(y))
end

# columns: species,k_M,E_m,lambda_max,g,T_A,T_typical,Ww_i,CTmax,CTmin,lat,class
rows = NamedTuple[]
for (i, line) in enumerate(eachline(TBL))
    i == 1 && continue
    f = split(line, ","); length(f) >= 12 || continue
    push!(rows, (kM = pf(f[2]), Em = pf(f[3]), lmax = pf(f[4]), g = pf(f[5]),
                 Ttyp = pf(f[7]), Ww = pf(f[8]), CTmax = pf(f[9]), CTmin = pf(f[10]),
                 lat = pf(f[11]), class = String(f[12])))
end

function arrs(rs)
    (kM = [r.kM for r in rs], lmax = [r.lmax for r in rs], g = [r.g for r in rs],
     Em = [r.Em for r in rs], Ww = [r.Ww for r in rs], Ttyp = [r.Ttyp for r in rs],
     CTmax = [r.CTmax for r in rs], CTmin = [r.CTmin for r in rs],
     breadth = [r.CTmax - r.CTmin for r in rs], lat = [abs(r.lat) for r in rs])
end

function report(rs, title; ectotherm = false)
    a = arrs(rs)
    println("\n=== ", title, " (n=", length(rs), ") ===")
    preds = (("k_M", log.(a.kM)), ("lambda_max", log.(a.lmax)), ("g", log.(a.g)), ("E_m (A0)", log.(a.Em)))
    targets = (("CTmax", a.CTmax), ("CTmin", a.CTmin), ("breadth", a.breadth))
    @printf("%-12s %14s %14s %14s   (raw Spearman)\n", "model \\ thermal", "CTmax", "CTmin", "breadth")
    for (pl, px) in preds
        cells = String[]
        for (_, ty) in targets
            m = isfinite.(px) .& isfinite.(ty); r = spear(px[m], ty[m])
            push!(cells, @sprintf("%+.3f %-2s", r, sig(r, count(m) - 2)))
        end
        @printf("%-12s %14s %14s %14s\n", pl, cells[1], cells[2], cells[3])
    end
    # pre-registered headline: capacity -> breadth, controlling |lat| + log size (+T_typ if ectotherm)
    println("-- partial: predictor -> breadth, controlling |lat|, log Ww_i", ectotherm ? ", T_typical --" : " --")
    for (pl, px) in preds
        ctrls = ectotherm ? Any[a.lat, log.(a.Ww), a.Ttyp] : Any[a.lat, log.(a.Ww)]
        m = isfinite.(px) .& isfinite.(a.breadth)
        for c in ctrls; m = m .& isfinite.(c); end
        r = prank(px[m], a.breadth[m], [c[m] for c in ctrls])
        @printf("   %-12s rho(->breadth | controls) = %+.3f %-2s  (n=%d)\n",
                pl, r, sig(r, count(m) - 2 - length(ctrls)), count(m))
    end
end

report(rows, "ALL taxa")
report(filter(r -> !(r.class in ("Mammalia", "Aves")), rows), "ECTOTHERMS only (Class != Mammalia/Aves)"; ectotherm = true)

println("\n* p<0.05, ** p<0.01. PRE-REGISTERED coherence probe (see header): a clear")
println("positive capacity->breadth partial = recovery capacity generalises to thermal")
println("resilience; a null = thermal tolerance is a separate axis (honest, bounding).")
