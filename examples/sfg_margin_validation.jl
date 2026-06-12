# ===========================================================================
# Scope-for-Growth margin validation (external validation #6 — the real margin test).
# Does the model's adaptive margin, driven by per-site contaminant pressure, track the
# measured energetic margin (SFG) of Mytilus edulis along the North Sea pollution
# gradient (Widdows et al. 1995, MEPS 127)? Same organisational level (individual
# energetics) => NO scale bridge; SFG is an independent measurement.
#
# Two layers:
#   (1) EMPIRICAL — which contaminants predict SFG (replicate Widdows' finding that
#       hydrocarbons dominate); establishes the gradient + the naive baseline.
#   (2) MODEL — for M. edulis, route contaminants to DEB axes (documented mode-of-action
#       below), build median-normalised toxic-unit pressures, run the margin-first point
#       API (compute_adaptive_margin_response), and correlate modeled margin with SFG.
#       The question that matters: does the MoA-routed margin track SFG, and how does it
#       compare to the single best raw contaminant predictor?
#
# HONEST CAVEATS (single-species design):
#   - With ONE species, the AmP capacity (A0, k_M, weights) is constant across sites, so
#     this validates the EROSION mechanism + MoA/mixture aggregation, NOT the capacity
#     weighting (that needs across-species SFG).
#   - Tissue body burden vs exposure-based EC50 is a known mismatch (critical body
#     residue); we use median-normalised relative burden as a threshold-free pressure
#     proxy, so potency differences between contaminants are NOT encoded.
#   - The contaminant->axis routing is a documented pMoA assignment, approximate.
#
#   julia +release --project=. examples/sfg_margin_validation.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const DIR = joinpath(@__DIR__, "..", "data", "external")
norm_site(s) = lowercase(replace(strip(s), r"[^A-Za-z0-9]" => ""))
parse_conc(x) = (x = strip(x); (x == "" || x == "nd") ? NaN :
                 (startswith(x, "<") ? parse(Float64, x[2:end]) : parse(Float64, x)))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
function spear(x, y)
    m = isfinite.(x) .& isfinite.(y)
    count(m) < 4 && return (NaN, 0)
    return (cor(ordinalrank(x[m]), ordinalrank(y[m])), count(m))
end
sig(r, n) = (n <= 3 ? "" : (t = r * sqrt((n - 2) / (1 - r^2)); abs(t) > 2.7 ? "**" : abs(t) > 2.0 ? "*" : ""))

function main()
# --- read SFG (skip comment lines) ---
sfg = Dict{String, Float64}()
for line in eachline(joinpath(DIR, "sfg_widdows1995_northsea.csv"))
    (startswith(line, "#") || startswith(line, "site,")) && continue
    f = split(line, ","); length(f) >= 6 || continue
    sfg[norm_site(f[1])] = parse(Float64, f[6])    # sfg_J_g_h
end

# --- read contaminants ---
chead = String[]; cont = Dict{String, Dict{String, Float64}}()
for line in eachline(joinpath(DIR, "sfg_widdows1995_contaminants.csv"))
    startswith(line, "#") && continue
    f = split(line, ",")
    if startswith(line, "site,"); chead = String.(f); continue; end
    length(f) == length(chead) || continue
    d = Dict{String, Float64}()
    for j in 2:length(f); d[chead[j]] = parse_conc(f[j]); end
    cont[norm_site(f[1])] = d
end

# --- join (exact normalised, else prefix) ---
function findsite(k, keys)
    k in keys && return k
    for kk in keys; (startswith(kk, k) || startswith(k, kk)) && return kk; end
    return nothing
end
sk = collect(keys(sfg))
sites = String[]; SFG = Float64[]; C = Dict{String, Vector{Float64}}()
for ck in keys(cont)
    s = findsite(ck, sk); s === nothing && continue
    push!(sites, ck); push!(SFG, sfg[s])
    for (name, v) in cont[ck]; push!(get!(C, name, Float64[]), v); end
end
@printf("joined sites: %d (Mytilus edulis, North Sea)\n", length(sites))

# ============================ (1) EMPIRICAL ============================
println("\n=== (1) EMPIRICAL: SFG vs contaminant (Spearman; negative = stress lowers SFG) ===")
for c in ("total_toxic_HC", "PAH_2_3ring", "THC_UVF", "TBT", "PCB", "Cu", "Cd", "Pb", "Zn", "Hg")
    haskey(C, c) || continue
    r, n = spear(C[c], SFG)
    @printf("  SFG ~ %-16s rho = %+.3f %-2s (n=%d)\n", c, r, sig(r, n), n)
end

# ============================ (2) MODEL ============================
# Documented pMoA routing (mussel toxicology):
#   assimilation <- hydrocarbons (PAH/THC suppress clearance/feeding; Widdows' mechanism)
#   maintenance  <- metals (Cd,Cu,Hg,Pb,Zn: metabolic/gill cost)
#   reproduction <- organotins + organochlorines (TBT imposex; PCB/OC reproductive tox)
#   growth       <- (minor) organotin DBT
# Pressure per contaminant = median-normalised relative burden (threshold-free toxic-unit
# proxy): tu = conc / median(conc). Axis pressure = mean tu over routed contaminants.
params = amp_species_deb_params(load_amp_species_library(), "Mytilus_edulis")
tu(c) = (v = C[c]; med = median(filter(isfinite, v)); [isfinite(x) ? x / med : 1.0 for x in v])
meanrows(vs...) = [mean(t) for t in zip(vs...)]

p_assim = tu("total_toxic_HC")
p_maint = meanrows(tu("Cd"), tu("Cu"), tu("Hg"), tu("Pb"), tu("Zn"))
p_repro = meanrows(tu("TBT"), tu("PCB"))
p_growth = tu("DBT")

margin = Float64[]; Qt = Float64[]
for i in eachindex(sites)
    res = compute_adaptive_margin_response((p_assim[i], p_maint[i], p_growth[i], p_repro[i]), params)
    push!(margin, res.A_t); push!(Qt, res.Q_t)
end

println("\n=== (2) MODEL: modeled adaptive margin vs measured SFG ===")
rM, nM = spear(margin, SFG)
@printf("  SFG ~ modeled margin A_t      rho = %+.3f %-2s (n=%d)  [+ = model tracks SFG]\n", rM, sig(rM, nM), nM)
# naive baselines for comparison
rHC, nHC = spear(C["total_toxic_HC"], SFG)
naive_load = meanrows(tu("total_toxic_HC"), tu("TBT"), tu("PCB"), tu("Cu"), tu("Cd"), tu("Pb"), tu("Zn"), tu("Hg"))
rL, nL = spear(naive_load, SFG)
@printf("  (baseline) SFG ~ total_toxic_HC          |rho| = %.3f\n", abs(rHC))
@printf("  (baseline) SFG ~ naive mean toxic-unit   |rho| = %.3f\n", abs(rL))
@printf("  MODEL margin                              |rho| = %.3f\n", abs(rM))

println("\nInterpretation: a strong positive SFG~margin (|rho| comparable to or above the")
println("best raw predictor) = the margin model reproduces the measured energetic-margin")
println("gradient through its pressure->axis->margin mechanism. With one species this")
println("tests the EROSION mechanism + MoA aggregation, not the capacity weighting.")
println("* |t|>2.0, ** |t|>2.7. See docs/notes/sfg_validation_status.md for caveats.")
end

main()
