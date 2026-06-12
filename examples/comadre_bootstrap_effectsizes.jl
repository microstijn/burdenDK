# ===========================================================================
# Effect sizes + multiple testing (roadmap Part 3) for the COMADRE external
# validation. Puts bootstrap (resample-over-species) 95% CIs on the headline
# partial-rank correlations and applies a Benjamini-Hochberg correction across
# the family of model quantities tested, so the manuscript can report intervals
# and an honest multiplicity adjustment rather than bare point estimates.
#
# READ-ONLY. Inputs (all committed):
#   data/external/comadre_recovery.csv         (recovery, gen, reactivity, attenuation)
#   data/external/comadre_amp_namemap.csv        (harmonised names)
#   data/external/amp_reproduction_rates.csv     (R_i, r_B, a_p, Ww_i)
#   data/AmP_Species_Library.json                (k_M, g via DEBAxisParams)
#   julia +release --project=. examples/comadre_bootstrap_effectsizes.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf, Random

const DIR = joinpath(@__DIR__, "..", "data", "external")

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
# partial rank correlation of x,y controlling any number of covariates (ranks).
function prank(x, y, ctrls::Vector)
    n = length(x)
    Z = hcat(ones(n), (ordinalrank(c) for c in ctrls)...)
    res(v) = (vr = ordinalrank(v); vr .- Z * (Z \ vr))
    rx, ry = res(x), res(y)
    return cor(rx, ry)
end
# two-sided p from a correlation via the t approximation (df = n - 2 - #controls).
function pval(r, df)
    df <= 0 && return NaN
    t = r * sqrt(df / (1 - r^2))
    # survival of |t| under Student-t via the regularized incomplete beta
    x = df / (df + t^2)
    return betainc_reg(df / 2, 0.5, x)   # = 2 * P(T > |t|)
end
# regularized incomplete beta I_x(a,b) via continued fraction (Numerical Recipes).
function betainc_reg(a, b, x)
    (x <= 0) && return 0.0
    (x >= 1) && return 1.0
    lbeta = lgamma(a) + lgamma(b) - lgamma(a + b)
    front = exp(a * log(x) + b * log(1 - x) - lbeta) / a
    f = 1.0; c = 1.0; d = 0.0
    for i in 0:200
        m = div(i, 2)
        if i == 0
            num = 1.0
        elseif i % 2 == 0
            num = (m * (b - m) * x) / ((a + 2m - 1) * (a + 2m))
        else
            num = -((a + m) * (a + b + m) * x) / ((a + 2m) * (a + 2m + 1))
        end
        d = 1.0 + num * d; abs(d) < 1e-30 && (d = 1e-30); d = 1 / d
        c = 1.0 + num / c; abs(c) < 1e-30 && (c = 1e-30)
        f *= d * c
        abs(1 - d * c) < 1e-10 && break
    end
    return front * (f - 1)
end
# Lanczos lgamma
function lgamma(z)
    g = 7.0
    c = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
    z < 0.5 && return log(pi / sin(pi * z)) - lgamma(1 - z)
    z -= 1
    a = c[1]; t = z + g + 0.5
    for i in 2:length(c); a += c[i] / (z + i - 1); end
    return 0.5 * log(2pi) + (z + 0.5) * log(t) - t + log(a)
end
bh(ps) = (n = length(ps); o = sortperm(ps); adj = similar(ps);
          m = 1.0; for k in n:-1:1; m = min(m, ps[o[k]] * n / k); adj[o[k]] = m; end; adj)

pf(x) = x == "" ? NaN : parse(Float64, x)

function load()
    com = Dict{String, NTuple{4, Float64}}()   # sp => (recovery, gen, reactivity, attenuation)
    for line in Iterators.drop(eachline(joinpath(DIR, "comadre_recovery.csv")), 1)
        f = split(line, ","); length(f) >= 9 || continue
        com[String(f[1])] = (pf(f[3]), pf(f[4]), pf(f[8]), pf(f[9]))
    end
    namemap = Dict{String, String}()
    for line in Iterators.drop(eachline(joinpath(DIR, "comadre_amp_namemap.csv")), 1)
        f = split(line, ","); length(f) >= 2 && f[2] != "" && (namemap[String(f[1])] = String(f[2]))
    end
    rates = Dict{String, NTuple{4, Float64}}()  # key => (R_i, r_B, Ww_i, a_p)
    for line in Iterators.drop(eachline(joinpath(DIR, "amp_reproduction_rates.csv")), 1)
        f = split(line, ","); length(f) >= 8 || continue
        rates[String(f[1])] = (pf(f[3]), pf(f[4]), pf(f[7]), pf(f[8]))
    end
    return com, namemap, rates
end

function main()
    com, namemap, rates = load()
    lib = load_amp_species_library()

    grouped = Dict{String, Vector{NTuple{4, Float64}}}()
    for (sp, v) in com
        k = get(namemap, sp, nothing); k === nothing && continue
        push!(get!(grouped, k, NTuple{4, Float64}[]), v)
    end

    kM = Float64[]; g = Float64[]; Ri = Float64[]; rB = Float64[]; ap = Float64[]; Ww = Float64[]
    recov = Float64[]; comp = Float64[]; resist = Float64[]; gen = Float64[]
    for (key, vs) in grouped
        haskey(lib, key) && haskey(rates, key) || continue
        local p; try; p = amp_species_deb_params(lib, key); catch; continue; end
        (0 < p.lambda_min <= p.lambda_max && isfinite(p.A0)) || continue
        r = rates[key]; all(isfinite, (r[1], r[2])) || continue
        push!(kM, p.lambda_min); push!(g, p.lambda_max / p.lambda_min)
        push!(Ri, r[1]); push!(rB, r[2]); push!(Ww, r[3]); push!(ap, r[4])
        push!(recov, mean(x[1] for x in vs)); push!(gen, mean(x[2] for x in vs))
        push!(comp, mean(x[3] for x in vs)); push!(resist, mean(x[4] for x in vs))
    end
    N = length(kM)
    @printf("species: %d\n\n", N)

    # Family of tests: (label, x, y, controls, n_controls) — gen-controlled partials.
    lg(v) = log.(v)
    tests = [
        ("k_M -> recovery        | gen",        lg(kM), recov,  [gen]),
        ("g   -> recovery        | gen (null)", lg(g),  recov,  [gen]),
        ("R_i -> compensation    | gen",        lg(Ri), comp,   [gen]),
        ("R_i -> compensation    | gen,mass",   lg(Ri), comp,   [gen, lg(Ww)]),
        ("a_p -> compensation    | gen,R_i",    lg(ap), comp,   [gen, lg(Ri)]),
        ("k_M -> resistance      | gen",        lg(kM), resist, [gen]),
        ("r_B -> recovery        | gen",        lg(rB), recov,  [gen]),
    ]

    # finite mask per test, then point estimate, analytic p, bootstrap CI.
    Random.seed!(20260612); B = 5000
    rows = NamedTuple[]
    for (lab, x, y, ctrls) in tests
        m = [all(isfinite, t) for t in zip(x, y, ctrls...)]
        xs, ys = x[m], y[m]; cs = [c[m] for c in ctrls]; n = count(m)
        r = prank(xs, ys, cs)
        p = pval(r, n - 2 - length(cs))
        boot = Float64[]
        idx = collect(1:n)
        for _ in 1:B
            s = rand(idx, n)
            push!(boot, prank(xs[s], ys[s], [c[s] for c in cs]))
        end
        sort!(boot)
        lo, hi = boot[max(1, round(Int, 0.025B))], boot[round(Int, 0.975B)]
        push!(rows, (lab = lab, r = r, lo = lo, hi = hi, p = p, n = n))
    end
    padj = bh([row.p for row in rows])

    @printf("%-34s %7s %18s %9s %9s %4s\n", "test (partial rank rho)", "rho", "95% CI (bootstrap)", "p", "p_BH", "n")
    for (row, pa) in zip(rows, padj)
        star = pa < 0.01 ? "**" : pa < 0.05 ? "*" : ""
        @printf("%-34s %+6.3f   [%+6.3f, %+6.3f] %9.4f %9.4f %4d %s\n",
                row.lab, row.r, row.lo, row.hi, row.p, pa, row.n, star)
    end
    println("\nBootstrap: $B resamples over species (seed 20260612), 95% percentile CI.")
    println("p_BH = Benjamini-Hochberg across the $(length(rows)) tests. * p_BH<0.05, ** <0.01.")
    println("A CI excluding 0 = the effect is robust to which species are sampled.")
end

main()
