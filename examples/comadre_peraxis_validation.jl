# ===========================================================================
# Per-axis validation (Idea B): do different DEB process rates predict different
# COMADRE demographic-resilience components? Tests the pre-registered mapping in
# docs/notes/comadre_peraxis_prereg.md:
#     k_M (maintenance)   -> recovery   (damping ratio)
#     R_i (reproduction)  -> compensation (reactivity, max col sum of A/lambda1)
#     r_B (growth)        -> recovery
#   resistance = attenuation (min col sum of A/lambda1).
#
# READ-ONLY. Inputs (all committed):
#   data/external/comadre_recovery.csv         (+reactivity,+attenuation cols)
#   data/external/comadre_amp_namemap.csv       (harmonised names)
#   data/external/amp_reproduction_rates.csv    (k_M,R_i,r_B,kap_R,Ww_i)
# Reports the rate x component rank-correlation matrix, raw and controlling
# generation time, and evaluates the pre-registered diagonal predictions.
#   julia +release --project=. examples/comadre_peraxis_validation.jl
# ===========================================================================

using Statistics, Printf

const DIR = joinpath(@__DIR__, "..", "data", "external")

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = cor(ordinalrank(x), ordinalrank(y))
pcor(a, b, c) = (rab = cor(a, b); rac = cor(a, c); rbc = cor(b, c);
                 (rab - rac * rbc) / sqrt((1 - rac^2) * (1 - rbc^2)))
pspear(x, y, z) = pcor(ordinalrank(x), ordinalrank(y), ordinalrank(z))
sig(r, df) = (df <= 1 ? "" : (t = r * sqrt(df / (1 - r^2)); abs(t) > 2.6 ? "**" : abs(t) > 1.97 ? "*" : ""))

pf(x) = x == "" ? NaN : parse(Float64, x)

# COMADRE: species => (logdamp, gen, reactivity, attenuation)
function read_comadre()
    d = Dict{String, NTuple{4, Float64}}()
    for line in Iterators.drop(eachline(joinpath(DIR, "comadre_recovery.csv")), 1)
        f = split(line, ","); length(f) >= 9 || continue
        d[String(f[1])] = (pf(f[3]), pf(f[4]), pf(f[8]), pf(f[9]))
    end
    return d
end

read_namemap() = (m = Dict{String, String}();
    for line in Iterators.drop(eachline(joinpath(DIR, "comadre_amp_namemap.csv")), 1)
        f = split(line, ","); length(f) >= 2 && f[2] != "" && (m[String(f[1])] = String(f[2]))
    end; m)

# AmP rates: amp_key => (kM, Ri, rB, Wwi)
function read_rates()
    d = Dict{String, NTuple{4, Float64}}()
    for line in Iterators.drop(eachline(joinpath(DIR, "amp_reproduction_rates.csv")), 1)
        f = split(line, ","); length(f) >= 7 || continue
        d[String(f[1])] = (pf(f[2]), pf(f[3]), pf(f[4]), pf(f[7]))
    end
    return d
end

function main()
    com = read_comadre(); namemap = read_namemap(); rates = read_rates()

    # group COMADRE rows by amp_key (collapse synonyms -> average, like the scalar analysis)
    grouped = Dict{String, Vector{NTuple{4, Float64}}}()
    for (sp, v) in com
        key = get(namemap, sp, nothing); key === nothing && continue
        push!(get!(grouped, key, NTuple{4, Float64}[]), v)
    end

    keys_used = String[]
    kM = Float64[]; Ri = Float64[]; rB = Float64[]; Wwi = Float64[]
    recov = Float64[]; comp = Float64[]; resist = Float64[]; gen = Float64[]
    for (key, vs) in grouped
        haskey(rates, key) || continue
        r = rates[key]
        all(isfinite, (r[1], r[2], r[3])) || continue
        push!(keys_used, key)
        push!(kM, r[1]); push!(Ri, r[2]); push!(rB, r[3]); push!(Wwi, r[4])
        push!(recov, mean(v[1] for v in vs))
        push!(gen, mean(x for x in (v[2] for v in vs)))   # may be NaN if all missing
        push!(comp, mean(v[3] for v in vs))
        push!(resist, mean(v[4] for v in vs))
    end
    @printf("matched species (model rates + COMADRE components): %d\n", length(keys_used))

    models = (("k_M  (maint)", log.(kM)), ("R_i  (repro)", log.(Ri)), ("r_B  (growth)", log.(rB)))
    comps  = (("resistance", resist), ("compensat.", comp), ("recovery", recov))
    # finite mask per component (reactivity/attenuation/gen can be missing)
    finite_all(xs...) = [all(isfinite, t) for t in zip(xs...)]

    println("\n=== RAW rank correlation (rate x component) ===")
    @printf("%-14s %12s %12s %12s\n", "rate \\ comp", "resistance", "compensat.", "recovery")
    for (mlab, x) in models
        row = String[]
        for (_, y) in comps
            m = finite_all(x, y); r = spear(x[m], y[m])
            push!(row, @sprintf("%+.3f %-2s", r, sig(r, count(m) - 2)))
        end
        @printf("%-14s %12s %12s %12s\n", mlab, row[1], row[2], row[3])
    end

    println("\n=== PARTIAL controlling generation time (the real test) ===")
    @printf("%-14s %12s %12s %12s\n", "rate \\ comp", "resistance", "compensat.", "recovery")
    parmat = zeros(3, 3)
    for (i, (mlab, x)) in enumerate(models)
        row = String[]
        for (j, (_, y)) in enumerate(comps)
            m = finite_all(x, y, gen); r = pspear(x[m], y[m], gen[m])
            parmat[i, j] = r
            push!(row, @sprintf("%+.3f %-2s", r, sig(r, count(m) - 3)))
        end
        @printf("%-14s %12s %12s %12s\n", mlab, row[1], row[2], row[3])
    end

    # pre-registered prediction checks (on the gen-controlled matrix)
    rownames = ["k_M", "R_i", "r_B"]; colnames = ["resistance", "compensation", "recovery"]
    println("\n=== Pre-registered predictions (gen-controlled partials) ===")
    # P1: R_i's strongest component is compensation
    ri = 2; bestcol = argmax(abs.(parmat[ri, :]))
    @printf("P1  R_i strongest component = %-13s (pred: compensation)  %s\n",
            colnames[bestcol], bestcol == 2 ? "PASS" : "fail")
    # P1b: compensation's strongest rate is R_i
    cc = 2; bestrow = argmax(abs.(parmat[:, cc]))
    @printf("P1b compensation strongest rate = %-6s (pred: R_i)         %s\n",
            rownames[bestrow], bestrow == 2 ? "PASS" : "fail")
    # P2: k_M's strongest component is recovery
    km = 1; bestcol2 = argmax(abs.(parmat[km, :]))
    @printf("P2  k_M strongest component = %-13s (pred: recovery)     %s\n",
            colnames[bestcol2], bestcol2 == 3 ? "PASS" : "fail")
    # P2b: recovery's strongest rate
    rc = 3; bestrow2 = argmax(abs.(parmat[:, rc]))
    @printf("P2b recovery strongest rate = %-6s (pred: k_M or r_B)    %s\n",
            rownames[bestrow2], bestrow2 in (1, 3) ? "PASS" : "fail")

    # secondary: also control body mass (Ww_i) alongside generation time, via
    # double residualisation on ranks (control gen first, then mass).
    println("\n=== R_i -> compensation, controlling gen AND body mass (Ww_i) ===")
    m = finite_all(log.(Ri), comp, gen, log.(max.(Wwi, eps())))
    if count(m) > 5
        ri_, co_, ge_, ww_ = log.(Ri)[m], comp[m], gen[m], log.(Wwi[m])
        # residual of each on (gen, mass) by rank, then correlate residuals
        function resid_on(v, c1, c2)
            X = hcat(ones(length(v)), ordinalrank(c1), ordinalrank(c2))
            v .- X * (X \ v)
        end
        rr = cor(resid_on(ordinalrank(ri_), ge_, ww_), resid_on(ordinalrank(co_), ge_, ww_))
        @printf("partial rank rho(R_i, compensation | gen, Ww_i) = %+.3f %s  (n=%d)\n",
                rr, sig(rr, count(m) - 4), count(m))
    end

    println("\n* p<0.05, ** p<0.01. compensation=reactivity, resistance=attenuation,")
    println("recovery=damping ratio (all log10). Rates log-transformed. See prereg note.")
end

main()
