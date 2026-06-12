# ===========================================================================
# Scope-for-Growth margin validation — SECOND species (Mytilus galloprovincialis).
# Replicates the Widdows/M. edulis margin test (examples/sfg_margin_validation.jl) on the
# Iberian SMP survey of Albentosa et al. 2012 (STOTEN 435-436:430-445; 2007 & 2008).
# Same organisational level (individual energetics) => NO scale bridge; SFG is independent.
#
# WHY a second species: the single-species Widdows test validated the EROSION mechanism +
# MoA aggregation but holds AmP capacity constant. Adding M. galloprovincialis is step 1
# toward the across-species CAPACITY test (route (a), chemical). This script's immediate
# job is REPLICATION: does the MoA-routed margin track measured SFG in a second species/
# region, as it did in the North Sea?
#
# OUTCOME (see docs/notes/sfg_validation_status.md): it does NOT replicate — and the
# AUTHORS explain why. Albentosa et al. show SFG here is dominated by CONDITION INDEX
# (SFG~CI r=-0.617***, R^2=51.7%) and AGE/shell thickness (r=-0.465**, R^2=26.4%);
# chemicals add only 16.95% of variance, and only the organochlorines (DDTs, chlordanes)
# act as toxicants — Zn is POSITIVELY related to SFG (a confound). So tissue burden here
# indexes food/condition, not toxic exposure, and a burden-driven margin anti-tracks SFG.
# This is a BOUNDING result (when burden!=exposure), not a tuning failure. We do NOT
# re-route to chase a negative (that would be p-hacking against an invariant).
#
# DESIGN NOTES specific to this dataset:
#   - Two surveys (2007, 2008). At repeated sites SFG is systematically LOWER in 2008 while
#     CPI is ~unchanged => a survey batch effect (the authors trace it to higher CI in
#     2008). So we report POOLED and PER-SURVEY.
#   - CPI (chemical pollution index) = the authors' own aggregate; published naive baseline.
#   - No organotins in this survey, so the growth axis has no contaminant route here
#     (Widdows routed DBT -> growth); growth pressure = baseline 1.0.
#   - CONFOUND CONTROL: for the 8 stations sampled in both years (16 records, Table 5) the
#     paper tabulates CI and ST, so we partial them out and ask whether the margin~SFG
#     relationship is rescued once the authors' confounders are removed.
#
#   julia +release --project=. examples/sfg_margin_validation_albentosa2012.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const DIR = joinpath(@__DIR__, "..", "data", "external")
const SFG_FILE  = joinpath(DIR, "sfg_albentosa2012_iberia.csv")
const CONT_FILE = joinpath(DIR, "sfg_albentosa2012_contaminants.csv")
const BIOM_FILE = joinpath(DIR, "sfg_albentosa2012_biometric_repeated.csv")

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

# Spearman partial correlation of x,y controlling for the columns in Zcols (rank-residuals).
function partial_spear(x, y, Zcols::Vector{Vector{Float64}})
    m = isfinite.(x) .& isfinite.(y)
    for z in Zcols; m = m .& isfinite.(z); end
    count(m) < length(Zcols) + 3 && return (NaN, count(m))
    cols = Vector{Float64}[ones(count(m))]
    for z in Zcols; push!(cols, ordinalrank(z[m])); end
    Z = reduce(hcat, cols)
    resid(v) = v .- Z * (Z \ v)
    ex = resid(ordinalrank(x[m])); ey = resid(ordinalrank(y[m]))
    return (cor(ex, ey), count(m))
end

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
# --- read SFG + survey ---
shead, srows = read_table(SFG_FILE)
isfg = findfirst(==("sfg_J_g_h"), shead)
isurv = findfirst(==("survey"), shead)
sfg = Dict{String, Float64}(); survey = Dict{String, String}()
for (site, f) in srows
    sfg[norm_site(site)] = parse(Float64, f[isfg])
    survey[norm_site(site)] = f[isurv]
end

# --- read contaminants ---
chead, crows = read_table(CONT_FILE)
cont = Dict{String, Dict{String, Float64}}()
for (site, f) in crows
    d = Dict{String, Float64}()
    for j in 2:length(f); d[chead[j]] = parse_conc(f[j]); end
    cont[norm_site(site)] = d
end

# --- join on exact survey-tagged site key ---
sites = String[]; SFG = Float64[]; SURV = String[]; C = Dict{String, Vector{Float64}}()
for ck in sort(collect(keys(cont)))
    haskey(sfg, ck) || continue
    push!(sites, ck); push!(SFG, sfg[ck]); push!(SURV, survey[ck])
    for (name, v) in cont[ck]; push!(get!(C, name, Float64[]), v); end
end
n07 = count(==("2007"), SURV); n08 = count(==("2008"), SURV)
@printf("joined sites: %d (Mytilus galloprovincialis, Iberian SMP)  [2007: %d, 2008: %d]\n",
        length(sites), n07, n08)

# helper: spearman over a survey subset (or all)
function spear_sub(x, y, mask)
    xs = x[mask]; ys = y[mask]; spear(xs, ys)
end
m_all = trues(length(sites)); m07 = SURV .== "2007"; m08 = SURV .== "2008"
report(label, x, y) = begin
    ra, na = spear_sub(x, y, m_all); r7, n7 = spear_sub(x, y, m07); r8, n8 = spear_sub(x, y, m08)
    @printf("  %-22s pooled %+.3f %-2s (n=%d) | 2007 %+.3f %-2s (n=%d) | 2008 %+.3f %-2s (n=%d)\n",
            label, ra, sig(ra, na), na, r7, sig(r7, n7), n7, r8, sig(r8, n8), n8)
end

# ============================ (1) EMPIRICAL ============================
println("\n=== (1) EMPIRICAL: SFG vs contaminant (Spearman; negative = stress lowers SFG) ===")
for c in ("PAH13", "PCB7", "DDTs", "chlordanes", "Cu", "Cd", "Pb", "Hg", "Zn", "Ar", "CPI")
    haskey(C, c) || continue
    report("SFG ~ " * c, C[c], SFG)
end

# ============================ (2) MODEL ============================
# Documented pMoA routing (parallels the Widdows harness):
#   assimilation <- PAHs (PAH13: hydrocarbons suppress clearance/feeding — Widdows' mechanism)
#   maintenance  <- metals (Cd,Cu,Hg,Pb,Zn: metabolic/gill cost)
#   reproduction <- organochlorines (DDTs, PCB7, chlordanes: reproductive tox)
#   growth       <- (no organotin in this survey -> baseline)
# Pressure per contaminant = median-normalised relative burden tu = conc/median(conc).
# Median is taken over the full pooled set so the pressure scale is consistent across
# surveys; per-survey correlations subset these same pressures.
params = amp_species_deb_params(load_amp_species_library(), "Mytilus_galloprovincialis")
tu(c) = (v = C[c]; med = median(filter(isfinite, v)); [isfinite(x) ? x / med : 1.0 for x in v])
meanrows(vs...) = [mean(t) for t in zip(vs...)]

p_assim  = tu("PAH13")
p_maint  = meanrows(tu("Cd"), tu("Cu"), tu("Hg"), tu("Pb"), tu("Zn"))
p_repro  = meanrows(tu("DDTs"), tu("PCB7"), tu("chlordanes"))
p_growth = fill(1.0, length(sites))

margin = Float64[]; Qt = Float64[]
for i in eachindex(sites)
    res = compute_adaptive_margin_response((p_assim[i], p_maint[i], p_growth[i], p_repro[i]), params)
    push!(margin, res.A_t); push!(Qt, res.Q_t)
end

println("\n=== (2) MODEL: modeled adaptive margin vs measured SFG (+ = model tracks SFG) ===")
report("margin A_t", margin, SFG)
println("\n  baselines (|rho|, pooled):")
naive_load = meanrows(tu("PAH13"), tu("PCB7"), tu("DDTs"), tu("chlordanes"),
                      tu("Cu"), tu("Cd"), tu("Pb"), tu("Hg"), tu("Zn"))
rBest, _ = spear(C["PAH13"], SFG)
rL, _ = spear(naive_load, SFG)
rCPI, _ = spear(C["CPI"], SFG)
rM, _ = spear(margin, SFG)
@printf("    SFG ~ PAH13 (best single)        |rho| = %.3f\n", abs(rBest))
@printf("    SFG ~ naive mean toxic-unit      |rho| = %.3f\n", abs(rL))
@printf("    SFG ~ CPI (published index)      |rho| = %.3f\n", abs(rCPI))
@printf("    MODEL margin                     |rho| = %.3f\n", abs(rM))

# --------- (2b) per-axis pressure diagnostics: which axis carries the (anti)signal? -------
# (signed; pressure is higher=worse, so a TOXIC axis should be NEGATIVE with SFG.)
println("\n=== (2b) AXIS DIAGNOSTIC: SFG vs each axis pressure (negative = behaves as toxicant) ===")
report("SFG ~ p_maint (metals)", p_maint, SFG)
report("SFG ~ p_repro (organochl)", p_repro, SFG)
report("SFG ~ p_assim (PAHs)", p_assim, SFG)
println("  -> metals POSITIVE (confound: burden indexes food/condition, not stress);")
println("     organochlorines NEGATIVE (the only genuine toxic signal — matches the paper).")

# --------- (3) CONFOUND CONTROL on the 8 repeated stations (Table 5: CI, ST) -------------
bhead, brows = read_table(BIOM_FILE)
ici = findfirst(==("CI"), bhead); ist = findfirst(==("ST_g_cm2"), bhead)
CImap = Dict{String, Float64}(); STmap = Dict{String, Float64}()
for (site, f) in brows
    CImap[norm_site(site)] = parse(Float64, f[ici]); STmap[norm_site(site)] = parse(Float64, f[ist])
end
CI = [get(CImap, s, NaN) for s in sites]; ST = [get(STmap, s, NaN) for s in sites]

println("\n=== (3) CONFOUND CONTROL on the 8 repeated stations (n=16; CI + ST from Table 5) ===")
rci, nci = spear(CI, SFG); rst, nst = spear(ST, SFG)
@printf("  SFG ~ CI (condition)          rho = %+.3f %-2s (n=%d)  [authors: -0.617***, R^2 51.7%%]\n", rci, sig(rci, nci), nci)
@printf("  SFG ~ ST (age proxy)          rho = %+.3f %-2s (n=%d)  [authors: -0.465**,  R^2 26.4%%]\n", rst, sig(rst, nst), nst)
# restrict the raw correlation to rows with CI/ST present, for a like-for-like comparison
mask16 = isfinite.(CI) .& isfinite.(ST)
rm16r, _ = spear(margin[mask16], SFG[mask16])
rp, np = partial_spear(margin, SFG, [CI, ST])
@printf("  SFG ~ margin (raw, n=%d)       rho = %+.3f %-2s\n", count(mask16), rm16r, sig(rm16r, count(mask16)))
@printf("  SFG ~ margin | CI, ST          rho = %+.3f %-2s (n=%d)  [does controlling confounds rescue it?]\n", rp, sig(rp, np), np)

println("\nInterpretation: the margin does NOT replicate here, and the authors show why —")
println("SFG is condition/age-driven; tissue burden indexes food, not toxic exposure, so")
println("the metal route is positively confounded and the burden-margin anti-tracks SFG.")
println("Controlling CI+ST does not rescue it (the confound is in the pressure proxy, not")
println("just additive). The only genuine toxic axis (organochlorines) IS negative, as the")
println("paper found. BOUNDING result: the margin replicates where burden indexes exposure")
println("(Widdows, +0.41) and fails where burden is food-confounded (here). NOT a tuning")
println("failure — routing was held fixed. Capacity weighting remains untested (needs")
println("gradients where burden->exposure holds across species).")
println("* |t|>2.0, ** |t|>2.7. See docs/notes/sfg_validation_status.md.")
end

main()
