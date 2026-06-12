# ===========================================================================
# External validation (#3): does the model predict COMADRE demographic recovery
# BEYOND raw pace-of-life AND beyond phylogeny (taxonomic clade)?
#
# READ-ONLY. Reads the COMADRE-derived recovery table
# (data/external/comadre_recovery.csv, produced by
# scripts/extract_comadre_recovery.jl), matches species to AmP via the
# harmonised name map (data/external/comadre_amp_namemap.csv, produced by
# scripts/resolve_comadre_amp_names.jl — recovers synonyms/typos beyond exact
# string match), and reports:
#   - raw rank correlation,
#   - partial controlling COMADRE generation time (pace-of-life),
#   - Order-controlled partial: within-clade, group-mean-centered, ALSO
#     controlling generation time (a tree-free phylogenetic proxy for PGLS).
#
# Rationale + caveats: docs/notes/comadre_partial_validation.md
#   julia +release --project=. examples/comadre_partial_validation.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics
using Printf

const CSV = joinpath(@__DIR__, "..", "data", "external", "comadre_recovery.csv")
const MAP = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_namemap.csv")

function read_comadre(path)
    sp = String[]; logdamp = Float64[]; gen = Float64[]; order = String[]
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, ","); length(f) >= 6 || continue
        push!(sp, f[1]); push!(logdamp, parse(Float64, f[3]))
        push!(gen, f[4] == "" ? NaN : parse(Float64, f[4])); push!(order, f[6])
    end
    return sp, logdamp, gen, order
end

# comadre_species => amp_key, resolved rows only (empty amp_key skipped).
function read_namemap(path)
    m = Dict{String, String}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, ","); length(f) >= 2 || continue
        f[2] == "" || (m[String(f[1])] = String(f[2]))
    end
    return m
end

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = cor(ordinalrank(x), ordinalrank(y))
pcor(a, b, c) = (rab = cor(a, b); rac = cor(a, c); rbc = cor(b, c);
                 (rab - rac * rbc) / sqrt((1 - rac^2) * (1 - rbc^2)))
pspear(x, y, z) = pcor(ordinalrank(x), ordinalrank(y), ordinalrank(z))
sig(r, df) = (t = r * sqrt(df / (1 - r^2)); abs(t) > 2.6 ? "**" : abs(t) > 1.97 ? "*" : "")

# within-group (Order) mean-centered partial correlation, controlling z.
function order_controlled(x, y, z, grp; min_group = 2)
    cnt = Dict{String, Int}(); for g in grp; cnt[g] = get(cnt, g, 0) + 1; end
    keep = [cnt[grp[i]] >= min_group for i in eachindex(grp)]
    xs, ys, zs, gs = x[keep], y[keep], z[keep], grp[keep]
    function ctr(v)
        r = ordinalrank(v); m = Dict(g => mean(r[gs .== g]) for g in unique(gs))
        return [r[i] - m[gs[i]] for i in eachindex(r)]
    end
    return pcor(ctr(xs), ctr(ys), ctr(zs)), count(keep), length(unique(gs))
end

function main()
    sp, logdamp, gen, order = read_comadre(CSV)
    namemap = read_namemap(MAP)
    lib = load_amp_species_library()

    # Group COMADRE rows by their resolved AmP key. Multiple COMADRE species can
    # resolve to one AmP species (e.g. "Chen caerulescens" + "Anser caerulescens"
    # snow goose); averaging their COMADRE quantities avoids pseudoreplication.
    grouped = Dict{String, Vector{Tuple{Float64, Float64, String}}}()
    for k in eachindex(sp)
        key = get(namemap, sp[k], nothing)
        key === nothing && continue
        push!(get!(grouped, key, Tuple{Float64, Float64, String}[]),
              (logdamp[k], gen[k], order[k]))
    end

    lamA0 = Float64[]; lammin = Float64[]; gval = Float64[]
    rec = Float64[]; cgen = Float64[]; ord = String[]
    for (key, rows) in grouped
        haskey(lib, key) || continue
        local p
        try; p = amp_species_deb_params(lib, key); catch; continue; end
        (0 < p.lambda_min <= p.lambda_max && isfinite(p.A0)) || continue
        gens = [r[2] for r in rows if isfinite(r[2])]
        push!(lamA0, restoring_force_from_margin(p.A0, p))
        push!(lammin, p.lambda_min); push!(gval, p.lambda_max / p.lambda_min)
        push!(rec, mean(r[1] for r in rows))
        push!(cgen, isempty(gens) ? NaN : mean(gens)); push!(ord, rows[1][3])
    end
    m = isfinite.(cgen); n = count(m)
    @printf("AmP-COMADRE matched species: %d  (with generation time: %d)\n\n", length(rec), n)

    println("Predicting COMADRE demographic recovery (log damping ratio):")
    @printf("%-24s %10s %16s %22s\n", "model quantity", "raw", "| gen.time", "| gen.time + Order")
    for (lab, x) in (("lambda(A0) recov. rate", lamA0), ("lambda_min = k_M", lammin),
                     ("g (amplification)", gval))
        raw = spear(x, rec)
        par = pspear(x[m], rec[m], cgen[m])
        oc, neff, nord = order_controlled(x[m], rec[m], cgen[m], ord[m])
        @printf("%-24s %7.3f %2s %12.3f %2s %16.3f %2s\n",
                lab, raw, sig(raw, length(rec) - 2), par, sig(par, n - 3),
                oc, sig(oc, neff - nord - 2))
    end
    oc, neff, nord = order_controlled(lammin[m], rec[m], cgen[m], ord[m])
    @printf("\nOrder control: %d species in %d multi-species Orders (effective n)\n", neff, nord)
    @printf("reference  rho(recovery, generation time) = %+.3f\n", spear(rec[m], cgen[m]))
    println("\n* p<0.05, ** p<0.01. The rightmost column is the strongest test: within")
    println("taxonomic Orders AND controlling generation time. Signal there = the model")
    println("predicts recovery beyond both pace-of-life and coarse phylogeny.")
end

main()
