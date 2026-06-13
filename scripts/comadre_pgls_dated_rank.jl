# ===========================================================================
# RANK-based dated-tree PGLS (phylogenetic Spearman). Companion to
# scripts/comadre_pgls_dated.jl, which tests the LOG-LINEAR form and finds k_M
# nulls under a generation-time covariate. The surviving COMADRE signal was always
# the RANK partial (Spearman +0.26**, within-Order +0.19*), which a parametric
# PGLS does not test. This script rank-transforms y, the predictor and generation
# time, then runs the SAME Pagel's-lambda PGLS on the ranks — i.e. a phylogenetic
# Spearman. The standardized-rank beta is the phylogeny-corrected rank coefficient.
#
# STANDALONE — needs Distributions (throwaway env):
#   julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add("Distributions")'
#   julia +release --project=<that-env> scripts/comadre_pgls_dated_rank.jl
#
# Reads:  data/external/comadre_amp_matched.csv, data/external/comadre_amp_dated_tree.nwk
# ===========================================================================

using LinearAlgebra, Statistics, Printf, Distributions

const DIR = joinpath(@__DIR__, "..", "data", "external")
const NWK = joinpath(DIR, "comadre_amp_dated_tree.nwk")
const TBL = joinpath(DIR, "comadre_amp_matched.csv")

function read_table()
    rows = Dict{String, NamedTuple}()
    for (i, line) in enumerate(eachline(TBL))
        i == 1 && continue
        f = split(line, ","); length(f) >= 7 || continue
        rows[String(f[1])] = (kM = parse(Float64, f[3]), lamA0 = parse(Float64, f[4]),
                              g = parse(Float64, f[5]), y = parse(Float64, f[6]),
                              gen = f[7] == "" ? NaN : parse(Float64, f[7]))
    end
    return rows
end

isfile(NWK) || error("Dated tree not found at $NWK — run scripts/comadre_pgls_dated.jl first.")

# --- Newick parser WITH branch lengths (identical to comadre_pgls_dated.jl) ---
mutable struct Node; label::String; len::Float64; children::Vector{Int}; parent::Int; end
function parse_newick(s)
    s = strip(replace(s, r"\s" => "")); endswith(s, ";") && (s = s[1:end-1])
    nodes = Node[]; i = Ref(1)
    function newnode(parent)
        push!(nodes, Node("", 0.0, Int[], parent)); id = length(nodes)
        parent > 0 && push!(nodes[parent].children, id); return id
    end
    stop = (',', ')', '(', ';', ':')
    function readtok()
        a = i[]; while i[] <= lastindex(s) && !(s[i[]] in stop); i[] = nextind(s, i[]); end
        return String(s[a:prevind(s, i[])])
    end
    function clade(parent)
        id = 0
        if s[i[]] == '('
            id = newnode(parent); i[] = nextind(s, i[])
            while true
                clade(id); c = s[i[]]; i[] = nextind(s, i[]); c == ')' && break
            end
            nodes[id].label = readtok()
        else
            id = newnode(parent); nodes[id].label = readtok()
        end
        if i[] <= lastindex(s) && s[i[]] == ':'
            i[] = nextind(s, i[]); nodes[id].len = something(tryparse(Float64, readtok()), 0.0)
        end
        return id
    end
    return nodes, clade(0)
end

function phylo_vcv(nodes, root)
    n = length(nodes)
    rootdist = fill(0.0, n); order = Int[]; st = [root]
    while !isempty(st); v = pop!(st); push!(order, v); append!(st, nodes[v].children); end
    for v in order; v == root || (rootdist[v] = rootdist[nodes[v].parent] + nodes[v].len); end
    tips = [v for v in 1:n if isempty(nodes[v].children)]
    chain = Dict(t => (c = Int[]; v = t; while v != 0; push!(c, v); v = nodes[v].parent; end; c) for t in tips)
    m = length(tips); C = Matrix{Float64}(undef, m, m)
    for a in 1:m
        S = Set(chain[tips[a]])
        for b in a:m
            mrca = tips[a]; for v in chain[tips[b]]; (v in S) && (mrca = v; break); end
            C[a, b] = C[b, a] = rootdist[mrca]
        end
    end
    norm(l) = (t = split(replace(l, " " => "_"), "_"); length(t) >= 2 ? string(t[1], "_", t[2]) : l)
    idx = Dict{String, Int}(); for k in 1:m; idx[norm(nodes[tips[k]].label)] = k; end
    return C, idx
end

function pgls(y, X, C)
    n = length(y); p = size(X, 2)
    fit(lam) = (V = lam .* C; for i in 1:n; V[i, i] = C[i, i]; end;
                Vi = inv(Symmetric(V)); XtVi = X' * Vi; beta = (XtVi * X) \ (XtVi * y);
                r = y - X * beta; s2 = (r'*Vi*r)/n;
                (-0.5*(n*log(2π)+logdet(Symmetric(V))+n*log(s2)+n), beta, s2, Vi, XtVi))
    gr = (sqrt(5)-1)/2; a, b = 0.0, 1.0; c = b-gr*(b-a); dd = a+gr*(b-a)
    fc = fit(c)[1]; fd = fit(dd)[1]
    for _ in 1:60
        if fc < fd; a = c; c = dd; fc = fd; dd = a+gr*(b-a); fd = fit(dd)[1]
        else; b = dd; dd = c; fd = fc; c = b-gr*(b-a); fc = fit(c)[1]; end
    end
    lam = (a+b)/2; _, beta, s2, Vi, XtVi = fit(lam)
    covb = s2 * inv(XtVi*X); se = sqrt.(diag(covb)); t = beta./se
    (beta, se, t, 2 .* ccdf.(TDist(n-p), abs.(t)), lam)
end
z(v) = (v .- mean(v)) ./ std(v)

# average (tied) ranks — rank is monotone so rank(x) == rank(log x)
function tiedrank(v)
    n = length(v); p = sortperm(v); r = zeros(Float64, n); i = 1
    while i <= n
        j = i; while j < n && v[p[j+1]] == v[p[i]]; j += 1; end
        avg = (i + j) / 2
        for k in i:j; r[p[k]] = avg; end
        i = j + 1
    end
    return r
end

function main()
    rows = read_table()
    nodes, root = parse_newick(read(NWK, String))
    C, idx = phylo_vcv(nodes, root)
    keys_in = String[]; cidx = Int[]
    for (k, _) in rows
        haskey(idx, k) && (push!(keys_in, k); push!(cidx, idx[k]))
    end
    isempty(keys_in) && error("no tree tips matched the matched-table species")
    Csub = C[cidx, cidx]
    kM = [rows[k].kM for k in keys_in]; lamA0 = [rows[k].lamA0 for k in keys_in]
    g = [rows[k].g for k in keys_in]; y = [rows[k].y for k in keys_in]; gen = [rows[k].gen for k in keys_in]
    @printf("RANK-based dated-tree PGLS (phylogenetic Spearman, Pagel's lambda ML).\n")
    @printf("matched species on dated tree: %d\n\n", length(y))
    quantities = (("lambda(A0)", lamA0), ("k_M", kM), ("g", g))
    ry = z(tiedrank(y))
    println("=== rank(comadre_log_damping) ~ rank(predictor) ===")
    @printf("%-12s %9s %8s %8s %8s %8s\n", "predictor", "beta*", "SE", "t", "p", "Pagel-l")
    for (lab, x) in quantities
        b, se, t, p, lam = pgls(ry, hcat(ones(length(y)), z(tiedrank(x))), Csub)
        @printf("%-12s %9.3f %8.3f %8.2f %8.4f %8.3f\n", lab, b[2], se[2], t[2], p[2], lam)
    end
    keep = findall(isfinite, gen); Ck = Csub[keep, keep]
    ryk = z(tiedrank(y[keep])); genk_r = z(tiedrank(gen[keep]))
    println("\n=== rank(comadre_log_damping) ~ rank(predictor) + rank(generation_time) (n=", length(keep), ") ===")
    @printf("%-12s %9s %8s %8s %8s %8s\n", "predictor", "beta*", "SE", "t", "p", "Pagel-l")
    for (lab, x) in quantities
        b, se, t, p, lam = pgls(ryk, hcat(ones(length(keep)), z(tiedrank(x[keep])), genk_r), Ck)
        @printf("%-12s %9.3f %8.3f %8.2f %8.4f %8.3f\n", lab, b[2], se[2], t[2], p[2], lam)
    end
    println("\nbeta* on standardized ranks = the phylogeny-corrected partial-Spearman coefficient.")
    println("Compare with the non-phylogenetic Spearman partials (k_M +0.264** on gen, +0.190* within-Order).")
end

main()
