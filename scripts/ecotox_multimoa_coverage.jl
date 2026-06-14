# ecotox_multimoa_coverage.jl -- FEASIBILITY scan for the across-axis capacity-weighting test.
#
# The model's distinctive content is a per-species kappa-driven axis weighting
#   w = (1/2, kappa/4, kappa/4, (1-kappa)/2)   over (assim, maint, growth, repro),
# with kappa = alpha_axes[3] (the DEB allocation fraction). Assimilation is FIXED at 1/2 for every
# species, so it carries NO cross-species weighting signal; the discriminating contrast is
# maintenance (proportional to kappa) vs reproduction (proportional to 1-kappa). Growth has no
# independent chemical support (diflubenzuron only), so the test is maintenance-vs-reproduction.
#
# The test needs the SAME AmP-resident species recurring across maintenance AND reproduction,
# AND those species must SPREAD in kappa. This script measures that from
# data/external/ecotox_multimoa_extract.csv (scripts/extract_ecotox_multimoa.awk), split by the
# panel TIER (core = defensible pMoA; stratum = AChE/DDT/PAH robustness, contested pMoA).
#
# Run:  julia +release --project=. scripts/ecotox_multimoa_coverage.jl
import JSON
using Statistics, Printf

const EXTRACT = "data/external/ecotox_multimoa_extract.csv"
const PANEL   = "data/ecotox_multimoa_panel.csv"
const POTENCY = Set(["LC50", "EC50", "IC50"])

# --- AmP: species key -> (kappa, L_m structural) ---
lib = JSON.parsefile("data/AmP_Species_Library.json")
amp = Dict{String,NTuple{2,Float64}}()
for (sp, v) in lib
    (v isa AbstractDict && haskey(v, "alpha_axes")) || continue
    aa = v["alpha_axes"]; length(aa) == 4 || continue
    kappa = float(aa[3])
    (isfinite(kappa) && 0.0 < kappa < 1.0) || continue
    Lm = NaN
    a = get(v, "auxiliary_metrics", Dict())
    (a isa AbstractDict && haskey(a, "L_m")) && (Lm = float(a["L_m"]))
    amp[sp] = (kappa, Lm)
end
println("AmP species with usable kappa: ", length(amp))

# --- panel: cas -> tier ---
tier_of = Dict{String,String}()
for ln in eachline(PANEL)
    f = split(ln, ","); (length(f) < 5 || f[1] == "cas") && continue
    tier_of[String(f[1])] = String(f[5])
end

genus_species(l) = (t = split(l, "_"); length(t) >= 2 ? string(t[1], "_", t[2]) : String(l))

# --- read extract once: keep (species, axis, cas, tier) potency rows matched to AmP ---
rows = NamedTuple[]   # (sp, axis, cas, tier)
n_rows = 0
for ln in eachline(EXTRACT)
    global n_rows
    startswith(ln, "cas,") && continue
    f = split(ln, ","); length(f) >= 11 || continue
    n_rows += 1
    String(f[7]) in POTENCY || continue
    cas = String(f[1]); ax = String(f[3]); key = genus_species(String(f[5]))
    haskey(amp, key) || continue
    push!(rows, (sp = key, axis = ax, cas = cas, tier = get(tier_of, cas, "?")))
end
@printf("extract rows: %d ; potency rows matched to AmP: %d\n", n_rows, length(rows))

# --- analysis for a given allowed tier set ---
function analyse(allowed::Set{String}, label::String)
    axes_of = Dict{String,Dict{String,Set{String}}}()   # sp -> axis -> {cas}
    for r in rows
        r.tier in allowed || continue
        push!(get!(get!(axes_of, r.sp, Dict{String,Set{String}}()), r.axis, Set{String}()), r.cas)
    end
    axsp(ax) = Set(sp for (sp, d) in axes_of if haskey(d, ax))
    M = axsp("maintenance"); R = axsp("reproduction")
    mr = sort(collect(intersect(M, R)); by = sp -> -amp[sp][1])
    println("\n", "="^70, "\n[$label]  tiers = ", join(sort(collect(allowed)), "+"))
    @printf("  maintenance species: %d   reproduction species: %d   M∩R: %d\n",
            length(M), length(R), length(mr))
    if !isempty(mr)
        ks = [amp[sp][1] for sp in mr]
        @printf("  kappa over M∩R:  min=%.3f  median=%.3f  max=%.3f  range=%.3f\n",
                minimum(ks), median(ks), maximum(ks), maximum(ks) - minimum(ks))
    end
    return axes_of, mr
end

core_axes, core_mr = analyse(Set(["core"]), "CORE (defensible pMoA)")
analyse(Set(["core", "stratum"]), "CORE + STRATUM (incl. AChE/DDT/PAH)")

# --- core candidate matrix (the honest test set) ---
println("\n=== CORE candidate species (in maintenance AND reproduction): kappa, L_m, #chem ===")
@printf("  %-30s %5s %9s  %s\n", "species", "kappa", "L_m(cm)", "n_chem(maint/repro)")
nchem(sp, ax) = haskey(core_axes[sp], ax) ? length(core_axes[sp][ax]) : 0
for sp in core_mr
    k, Lm = amp[sp]
    @printf("  %-30s %5.3f %9.3g  %d/%d\n", sp, k, Lm, nchem(sp, "maintenance"), nchem(sp, "reproduction"))
end

open("data/external/ecotox_multimoa_core_candidates.csv", "w") do io
    println(io, "species,kappa,L_m,n_maint_core,n_repro_core")
    for sp in core_mr
        k, Lm = amp[sp]
        println(io, join([sp, round(k, digits = 4), round(Lm, sigdigits = 5),
            nchem(sp, "maintenance"), nchem(sp, "reproduction")], ","))
    end
end
println("\nwrote data/external/ecotox_multimoa_core_candidates.csv")
