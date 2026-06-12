# ===========================================================================
# Phylogenetic generalized least squares (PGLS) validation — Idea A. Tests
# whether the model's recovery quantities predict COMADRE demographic recovery
# beyond pace-of-life AND a REAL phylogeny (vs. the taxonomic-Order proxy used
# in examples/comadre_partial_validation.jl).
#
# Pure-Julia PGLS (no phylo-package dependency): parse the OTL induced-subtree
# Newick, assign Grafen (1989) ultrametric branch lengths, build the phylogenetic
# variance-covariance matrix, and fit GLS with Pagel's lambda estimated by ML.
# Subsetting the full-tree VCV to the generation-time subset is exact under the
# BM model (marginalization consistency), so the tree is not re-pruned.
#
# STANDALONE — needs Distributions (for the t p-value) in a throwaway env:
#   julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add("Distributions")'
#   julia +release --project=<that-env> scripts/comadre_pgls.jl
#
# Reads:  data/external/comadre_amp_matched.csv  (model + COMADRE + taxonomy)
#         data/external/comadre_amp_tree.nwk      (OTL topology)
#         data/external/comadre_amp_tree_map.csv  (amp_key -> ott_id)
#
# OTL gives TOPOLOGY ONLY; Grafen branch lengths are a coarse ultrametric proxy.
# Pagel's lambda (estimated) absorbs how much residual phylogenetic signal the
# branch lengths imply. A dated vertebrate tree (VertLife/TimeTree) is the
# stronger follow-up (needs a manual download — see the validation roadmap).
# ===========================================================================

using LinearAlgebra, Statistics, Printf, Distributions

const DIR = joinpath(@__DIR__, "..", "data", "external")

# --- Newick parser (no branch lengths, internal node labels present) ------
mutable struct Node; label::String; children::Vector{Int}; parent::Int; end

function parse_newick(s)
    s = strip(replace(s, r"\s" => ""))
    endswith(s, ";") && (s = s[1:end-1])
    nodes = Node[]; i = Ref(1)
    function newnode(label, parent)
        push!(nodes, Node(label, Int[], parent))
        id = length(nodes)
        parent > 0 && push!(nodes[parent].children, id)
        return id
    end
    function readlabel()
        a = i[]
        while i[] <= lastindex(s) && !(s[i[]] in (',', ')', '(', ';')); i[] = nextind(s, i[]); end
        return String(s[a:prevind(s, i[])])
    end
    function clade(parent)
        if s[i[]] == '('
            id = newnode("", parent); i[] = nextind(s, i[])
            while true
                clade(id)
                c = s[i[]]; i[] = nextind(s, i[])
                c == ')' && break
            end
            nodes[id].label = readlabel()
            return id
        else
            return newnode(readlabel(), parent)
        end
    end
    root = clade(0)
    return nodes, root
end

# --- Grafen branch lengths -> root-distance per node ----------------------
# kappa = descendant-tip count; height a = kappa-1; H = (a/a_root)^rho (rho=1);
# distance from root = 1 - H (tips at 1, root at 0). Ultrametric.
function grafen_rootdist(nodes, root; rho = 1.0)
    n = length(nodes)
    istip = [isempty(nodes[k].children) for k in 1:n]
    kappa = zeros(Int, n)
    order = Int[]; stack = [root]            # post-order
    while !isempty(stack); v = pop!(stack); push!(order, v); append!(stack, nodes[v].children); end
    for v in reverse(order)
        kappa[v] = istip[v] ? 1 : sum(kappa[c] for c in nodes[v].children)
    end
    aroot = kappa[root] - 1
    H = [(kappa[v] - 1) / aroot for v in 1:n] .^ rho
    return 1 .- H, istip
end

# --- phylogenetic VCV over tips -------------------------------------------
function phylo_vcv(nodes, root)
    rootdist, istip = grafen_rootdist(nodes, root)
    tips = findall(istip)
    # ancestor chain (tip -> root) per tip
    chain = Dict{Int, Vector{Int}}()
    for t in tips
        c = Int[]; v = t
        while v != 0; push!(c, v); v = nodes[v].parent; end
        chain[t] = c
    end
    m = length(tips)
    C = Matrix{Float64}(undef, m, m)
    for a in 1:m
        ta = tips[a]; S = Set(chain[ta])
        for b in a:m
            tb = tips[b]
            mrca = ta
            for v in chain[tb]; if v in S; mrca = v; break; end; end
            C[a, b] = C[b, a] = rootdist[mrca]
        end
    end
    # tip label -> matrix index
    idx = Dict(nodes[tips[k]].label => k for k in 1:m)
    return C, idx
end

# --- I/O ------------------------------------------------------------------
function read_matched(path)
    rows = Dict{String, NamedTuple}()
    hdr = true
    for line in eachline(path)
        if hdr; hdr = false; continue; end
        f = split(line, ","); length(f) >= 7 || continue
        rows[String(f[1])] = (kM = parse(Float64, f[3]), lamA0 = parse(Float64, f[4]),
                              g = parse(Float64, f[5]), y = parse(Float64, f[6]),
                              gen = f[7] == "" ? NaN : parse(Float64, f[7]))
    end
    return rows
end

# amp_key -> ott_id (placed only)
function read_treemap(path)
    m = Dict{String, String}()
    hdr = true
    for line in eachline(path)
        if hdr; hdr = false; continue; end
        f = split(line, ","); length(f) >= 3 || continue
        f[2] != "" && f[3] == "true" && (m[String(f[1])] = String(f[2]))
    end
    return m
end

# --- PGLS with Pagel's lambda (ML) ----------------------------------------
# Returns (beta, se, tval, pval, lambda, n) for each column of X (intercept first).
function pgls(y, X, C)
    n = length(y); p = size(X, 2)
    function fit(lam)
        V = lam .* C; for i in 1:n; V[i, i] = C[i, i]; end   # scale off-diagonals only
        Vi = inv(Symmetric(V))
        XtVi = X' * Vi
        beta = (XtVi * X) \ (XtVi * y)
        r = y - X * beta
        rss = (r' * Vi * r)
        s2 = rss / n
        ll = -0.5 * (n * log(2π) + logdet(Symmetric(V)) + n * log(s2) + n)
        return ll, beta, s2, Vi, XtVi
    end
    # golden-section maximize ll over lambda in [0,1]
    gr = (sqrt(5) - 1) / 2; a, b = 0.0, 1.0
    c = b - gr * (b - a); d = a + gr * (b - a)
    fc = fit(c)[1]; fd = fit(d)[1]
    for _ in 1:60
        if fc < fd; a = c; c = d; fc = fd; d = a + gr * (b - a); fd = fit(d)[1]
        else; b = d; d = c; fd = fc; c = b - gr * (b - a); fc = fit(c)[1]; end
    end
    lam = (a + b) / 2
    ll, beta, s2, Vi, XtVi = fit(lam)
    covb = s2 * inv(XtVi * X)
    se = sqrt.(diag(covb))
    tval = beta ./ se
    pval = 2 .* ccdf.(TDist(n - p), abs.(tval))
    return beta, se, tval, pval, lam, n
end

function loglik_at(y, X, C, lam)
    n = length(y); V = lam .* C; for i in 1:n; V[i, i] = C[i, i]; end
    Vi = inv(Symmetric(V)); XtVi = X' * Vi
    beta = (XtVi * X) \ (XtVi * y); r = y - X * beta
    s2 = (r' * Vi * r) / n
    return -0.5 * (n * log(2π) + logdet(Symmetric(V)) + n * log(s2) + n)
end

z(v) = (v .- mean(v)) ./ std(v)   # standardize predictor for comparable betas

function main()
    rows = read_matched(joinpath(DIR, "comadre_amp_matched.csv"))
    treemap = read_treemap(joinpath(DIR, "comadre_amp_tree_map.csv"))
    nwk = read(joinpath(DIR, "comadre_amp_tree.nwk"), String)
    nodes, root = parse_newick(nwk)
    C, tipidx = phylo_vcv(nodes, root)

    # build aligned arrays: only amp_keys present in matched table AND on tree
    # tip labels look like Genus_species_ott12345; map via ott_id suffix.
    ott_to_key = Dict(v => k for (k, v) in treemap)
    keys_in_order = String[]; cidx = Int[]
    for (label, j) in tipidx
        mm = match(r"^(.*)_ott(\d+)$", label)
        mm === nothing && continue
        binom, ott = String(mm.captures[1]), String(mm.captures[2])
        # primary: match by OTT id (principled); fallback: OTL kept the binomial
        # name but remapped the ott -> match the Genus_species directly to amp_key.
        k = get(ott_to_key, ott, haskey(rows, binom) ? binom : nothing)
        (k === nothing || !haskey(rows, k)) && continue
        push!(keys_in_order, k); push!(cidx, j)
    end
    Csub = C[cidx, cidx]
    kM   = [rows[k].kM    for k in keys_in_order]
    lamA0 = [rows[k].lamA0 for k in keys_in_order]
    gg   = [rows[k].g     for k in keys_in_order]
    y    = [rows[k].y     for k in keys_in_order]
    gen  = [rows[k].gen   for k in keys_in_order]

    @printf("PGLS on OTL tree (Grafen branch lengths, Pagel's lambda ML).\n")
    @printf("species on tree & in table: %d\n\n", length(y))

    quantities = (("lambda(A0)", lamA0), ("k_M", kM), ("g", gg))

    println("=== Model: comadre_log_damping ~ predictor (no pace-of-life control) ===")
    @printf("%-12s %9s %8s %8s %8s %8s\n", "predictor", "beta*", "SE", "t", "p", "Pagel-l")
    for (lab, x) in quantities
        X = hcat(ones(length(y)), z(log.(x)))
        b, se, t, p, lam, _ = pgls(y, X, Csub)
        @printf("%-12s %9.3f %8.3f %8.2f %8.4f %8.3f %s\n", lab, b[2], se[2], t[2], p[2], lam, stars(p[2]))
    end

    # generation-time-controlled (drop NaN gen; subset VCV — exact under BM)
    keep = findall(isfinite, gen)
    Ck = Csub[keep, keep]; yk = y[keep]; genk = gen[keep]
    println("\n=== Model: comadre_log_damping ~ predictor + generation_time ===")
    @printf("(n = %d with finite generation time)\n", length(keep))
    @printf("%-12s %9s %8s %8s %8s %8s\n", "predictor", "beta*", "SE", "t", "p", "Pagel-l")
    for (lab, x) in quantities
        xk = x[keep]
        X = hcat(ones(length(keep)), z(log.(xk)), z(log.(genk)))
        b, se, t, p, lam, _ = pgls(yk, X, Ck)
        @printf("%-12s %9.3f %8.3f %8.2f %8.4f %8.3f %s\n", lab, b[2], se[2], t[2], p[2], lam, stars(p[2]))
    end
    println("\nbeta* = standardized partial slope (per SD of log predictor).")
    println("* p<0.05, ** p<0.01. Pagel-l = ML phylogenetic signal in residuals (0=none,1=BM).")

    # Verify the lambda~0 estimate is a real ML optimum (not an optimizer artifact)
    # and that the OTL+Grafen tree simply carries little signal: profile logL.
    println("\n=== Pagel's lambda profile (logL of k_M + gen model) — sanity check ===")
    Xkm = hcat(ones(length(keep)), z(log.(kM[keep])), z(log.(genk)))
    @printf("%-8s %12s\n", "lambda", "logL")
    for lam in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
        @printf("%-8.2f %12.3f\n", lam, loglik_at(yk, Xkm, Ck, lam))
    end
    println("logL flat/declining in lambda => OTL+Grafen tree adds ~no phylogenetic")
    println("structure; this PGLS reduces to OLS. A DATED tree is the stronger test.")
end

stars(p) = p < 0.01 ? "**" : p < 0.05 ? "*" : ""

main()
