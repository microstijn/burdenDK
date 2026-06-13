# ===========================================================================
# Scope-for-Growth margin validation — THIRD study (Widdows et al. 2002, Irish Sea).
# Mytilus edulis, eastern Irish Sea (Phase I, 1996), 23 sites. Same group/method/species as
# the validated Widdows 1995 (examples/sfg_margin_validation.jl); a second, larger,
# independent gradient of the same kind. Same organisational level => NO scale bridge.
#
# WHY: third point in the across-species SFG programme (Widdows 1995 = +0.41 validated;
# Albentosa 2012 = -0.11 bounding/confounded). This is the screened "exposure-driven" mussel
# gradient the post-Albentosa rule asks for: hydrocarbon/industrial inputs, not a food
# gradient. The paper pre-establishes the mechanism: SFG = 26.8 - 9.4*log(toxic HC), i.e.
# hydrocarbons (our assimilation route) dominate; metals are non-toxic (As is POSITIVELY
# related to SFG, like Albentosa's Zn). Question: does the MoA-routed margin track SFG here?
#
# DATA CAVEAT: SFG is DIGITIZED from Fig 2A (no numeric table in the paper), ~+/-1 J/g/h,
# anchored to 4 Table-2 values -> RANK statistics only. Contaminant metals/TBT/PAH are
# reliable; sparse organochlorines (minor repro axis) lower-confidence. See the CSV headers.
#
#   julia +release --project=. examples/sfg_margin_validation_widdows2002.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const DIR = joinpath(@__DIR__, "..", "data", "external")
const SFG_FILE  = joinpath(DIR, "sfg_widdows2002_irishsea.csv")
const CONT_FILE = joinpath(DIR, "sfg_widdows2002_contaminants.csv")

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

function read_table(path)
    head = String[]; rows = Vector{Tuple{String, Vector{String}}}()
    for line in eachline(path)
        startswith(line, "#") && continue
        f = String.(split(line, ","))
        if startswith(line, "site,"); head = f; continue; end
        length(f) == length(head) || continue
        push!(rows, (f[1], f))
    end
    return head, rows
end

function main()
# --- read SFG ---
shead, srows = read_table(SFG_FILE)
isfg = findfirst(==("sfg_J_g_h"), shead)
sfg = Dict{String, Float64}()
for (site, f) in srows; sfg[norm_site(site)] = parse(Float64, f[isfg]); end

# --- read contaminants ---
chead, crows = read_table(CONT_FILE)
cont = Dict{String, Dict{String, Float64}}()
for (site, f) in crows
    d = Dict{String, Float64}()
    for j in 2:length(f); d[chead[j]] = parse_conc(f[j]); end
    cont[norm_site(site)] = d
end

# --- join ---
sites = String[]; SFG = Float64[]; C = Dict{String, Vector{Float64}}()
for ck in sort(collect(keys(cont)))
    haskey(sfg, ck) || continue
    push!(sites, ck); push!(SFG, sfg[ck])
    for (name, v) in cont[ck]; push!(get!(C, name, Float64[]), v); end
end
@printf("joined sites: %d (Mytilus edulis, eastern Irish Sea, Phase I 1996)\n", length(sites))

# ============================ (1) EMPIRICAL ============================
println("\n=== (1) EMPIRICAL: SFG vs contaminant (Spearman; negative = stress lowers SFG) ===")
for c in ("PAH", "TBT", "PCB", "DDE", "Cu", "Cd", "Pb", "Hg", "Zn", "As", "Se")
    haskey(C, c) || continue
    r, n = spear(C[c], SFG)
    @printf("  SFG ~ %-6s rho = %+.3f %-2s (n=%d)\n", c, r, sig(r, n), n)
end

# ============================ (2) MODEL ============================
# Documented pMoA routing (parallels Widdows 1995 / Albentosa):
#   assimilation <- PAH (2+3-ring toxic hydrocarbons; the validated dominant driver)
#   maintenance  <- metals (Cd,Cu,Hg,Pb,Zn)
#   reproduction <- organotin + organochlorines (TBT, PCB, DDE)
#   growth       <- DBT
# Pressure per contaminant = median-normalised relative burden tu = conc/median(conc).
params = amp_species_deb_params(load_amp_species_library(), "Mytilus_edulis")
tu(c) = (v = C[c]; med = median(filter(isfinite, v)); [isfinite(x) ? x / med : 1.0 for x in v])
meanrows(vs...) = [mean(t) for t in zip(vs...)]

p_assim  = tu("PAH")
p_maint  = meanrows(tu("Cd"), tu("Cu"), tu("Hg"), tu("Pb"), tu("Zn"))
p_repro  = meanrows(tu("TBT"), tu("PCB"), tu("DDE"))
p_growth = tu("DBT")

margin = Float64[]
for i in eachindex(sites)
    res = compute_adaptive_margin_response((p_assim[i], p_maint[i], p_growth[i], p_repro[i]), params)
    push!(margin, res.A_t)
end

println("\n=== (2) MODEL: modeled adaptive margin vs measured SFG (+ = model tracks SFG) ===")
rM, nM = spear(margin, SFG)
@printf("  SFG ~ modeled margin A_t   rho = %+.3f %-2s (n=%d)\n", rM, sig(rM, nM), nM)
naive_load = meanrows(tu("PAH"), tu("TBT"), tu("PCB"), tu("DDE"),
                      tu("Cu"), tu("Cd"), tu("Pb"), tu("Hg"), tu("Zn"))
rPAH, _ = spear(C["PAH"], SFG); rL, _ = spear(naive_load, SFG)
@printf("\n  baselines (|rho|):\n")
@printf("    SFG ~ PAH (best single, the validated axis)   |rho| = %.3f\n", abs(rPAH))
@printf("    SFG ~ naive mean toxic-unit                   |rho| = %.3f\n", abs(rL))
@printf("    MODEL margin                                  |rho| = %.3f\n", abs(rM))

# --------- (2b) axis diagnostic: which axis carries the (anti)signal? -------
println("\n=== (2b) AXIS DIAGNOSTIC: SFG vs each axis pressure (negative = behaves as toxicant) ===")
for (lab, p) in (("p_assim (PAH)", p_assim), ("p_maint (metals)", p_maint), ("p_repro (TBT+OC)", p_repro))
    r, n = spear(p, SFG); @printf("  SFG ~ %-18s rho = %+.3f %-2s (n=%d)\n", lab, r, sig(r, n), n)
end

println("\nInterpretation: a positive SFG~margin replicates the Widdows-1995 result in a second,")
println("larger, same-method M. edulis gradient. Expect the PAH/assimilation axis to carry the")
println("signal (negative) and metals to be a weak/positive confound (the paper's own finding:")
println("SFG = 26.8 - 9.4*log(toxic HC); As positively related). Caveat: single-marker")
println("correlations weaken over the large Irish-Sea scale (per the paper). SFG digitized from")
println("Fig 2A -> rank only. * |t|>2.0, ** |t|>2.7. See docs/notes/sfg_validation_status.md.")
end

main()
