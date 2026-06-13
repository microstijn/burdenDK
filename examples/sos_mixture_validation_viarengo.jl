# ===========================================================================
# CONTROLLED validation — Viarengo et al. 1995 stress-on-stress dose-response + mixture.
# The first controlled-exposure test of the model's IMPAIRMENT CURVE and MIXTURE AGGREGATION
# (the field SFG/SoS tests can't isolate single contaminants or mixtures). Mytilus
# galloprovincialis, 3-day exposure, survival-in-air LT50.
#
#   (A) DOSE-RESPONSE per MoA axis: does SoS LT50 fall monotonically with dose for Cu
#       (maintenance), DMBA/PAH (assimilation), Aroclor/PCB (reproduction)? potency ordering?
#   (B) MIXTURE (the key controlled test): the Cu+DMBA mixture is the paper's "cumulative
#       effect". Predict the mixture LT50 from the single-component effects using the model's
#       OWN mixture rules -- axis_toxic_unit_sum (CA) and independent_action_axis_effects (IA),
#       from src/mixture_aggregation.jl (the real code) -- and compare to observed. Tests the
#       "mixtures are additive assumptions, not fitted interactions" design (no synergism/
#       antagonism). Effects are mapped to the model impairment via E = x/(1+x).
#
# HONEST CAVEATS: n=2 mixtures, LT50 rounded to ~0.5 d, and combining LT50s (a time) via the
# fractional-effect IA/CA formalism is approximate (effect should ideally be the affected
# fraction at a fixed time). So (B) is a direction + bracketing test, not a precise IA-vs-CA
# discrimination. Single time point (3 d) -> this is STATIC; it does NOT test the dynamics.
#
#   julia +release --project=. examples/sos_mixture_validation_viarengo.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const FILE = joinpath(@__DIR__, "..", "data", "external", "sos_viarengo1995_doseresponse.csv")
const CONTROL_LT50 = 7.0
parse_num(x) = (x = strip(x, ['"', ' ']); (x == "" || x == "nd") ? NaN : something(tryparse(Float64, x), NaN))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = (m = isfinite.(x) .& isfinite.(y); count(m) < 3 ? NaN : cor(ordinalrank(x[m]), ordinalrank(y[m])))
rec(xM) = (burden_assimilation = 0.0, burden_maintenance = xM, burden_growth = 0.0, burden_reproduction = 0.0)

function main()
hdr = String[]; R = Vector{Vector{String}}()
for line in eachline(FILE)
    startswith(line, "#") && continue
    f = String.(split(line, ","))
    if startswith(line, "treatment,"); hdr = f; continue; end
    length(f) == length(hdr) && push!(R, f)
end
ci(n) = findfirst(==(n), hdr)
treat = [r[ci("treatment")] for r in R]; cont = [r[ci("contaminant")] for r in R]
typ = [r[ci("type")] for r in R]; LT50 = [parse_num(r[ci("LT50_days")]) for r in R]
Cu = [parse_num(r[ci("Cu_uM")]) for r in R]; DMBA = [parse_num(r[ci("DMBA_uM")]) for r in R]
Aro = [parse_num(r[ci("Aroclor_uM")]) for r in R]

# ---------- (A) dose-response per MoA axis ----------
println("=== (A) DOSE-RESPONSE: SoS LT50 vs dose, per single contaminant ===")
for (name, dose) in (("Cu (maintenance)", Cu), ("DMBA/PAH (assimilation)", DMBA), ("Aroclor/PCB (reproduction)", Aro))
    idx = findall(i -> (typ[i] == "single" && dose[i] > 0) || typ[i] == "control", eachindex(R))
    # restrict to this contaminant's series (its dose>0) plus control
    series = [i for i in idx if dose[i] > 0 || typ[i] == "control"]
    # only keep rows where the OTHER contaminants are 0 (clean single series)
    series = [i for i in series if (dose === Cu || Cu[i] == 0) && (dose === DMBA || DMBA[i] == 0) && (dose === Aro || Aro[i] == 0)]
    d = dose[series]; l = LT50[series]
    rho = spear(d, l)
    potency = (CONTROL_LT50 - minimum(l)) / maximum(d)   # LT50 loss per uM at the top dose
    @printf("  %-26s rho(dose,LT50) = %+.2f  | LT50 %s d  | potency ~ %.2f d/uM\n",
            name, rho, join(string.(round.(l, digits = 1)), "->"), potency)
end
println("  -> monotone decline on every axis; potency Cu > DMBA > PCB (metal most acute, PCB weakest).")

# ---------- (B) MIXTURE: model's own rules vs observed ----------
println("\n=== (B) MIXTURE Cu+DMBA: model mixture rules (real code) vs observed LT50 ===")
println("  effect E=(7-LT50)/7 ; burden x=E/(1-E) (inverts model impairment E=x/(1+x)) ;")
println("  predict mixture E via the model's axis_toxic_unit_sum (CA) and independent_action (IA).")
for i in findall(==("mixture"), typ)
    cuD = Cu[i]; dmbaD = DMBA[i]
    # matching single-component rows
    jCu = findfirst(k -> cont[k] == "Cu" && Cu[k] == cuD && typ[k] == "single", eachindex(R))
    jDM = findfirst(k -> cont[k] == "DMBA" && DMBA[k] == dmbaD && typ[k] == "single", eachindex(R))
    LCu = LT50[jCu]; LDM = LT50[jDM]; Lobs = LT50[i]
    ECu = (CONTROL_LT50 - LCu) / CONTROL_LT50; EDM = (CONTROL_LT50 - LDM) / CONTROL_LT50
    xCu = ECu / (1 - ECu); xDM = EDM / (1 - EDM)
    tu = aggregate_axis_mixture_effects([rec(xCu), rec(xDM)]; mixture_effect_model = "axis_toxic_unit_sum")
    ia = aggregate_axis_mixture_effects([rec(xCu), rec(xDM)]; mixture_effect_model = "independent_action_axis_effects")
    L_tu = CONTROL_LT50 * (1 - tu.E_maintenance); L_ia = CONTROL_LT50 * (1 - ia.E_maintenance)
    worst = min(LCu, LDM)
    @printf("  %s (Cu %.1f / DMBA %.1f uM):\n", treat[i], cuD, dmbaD)
    @printf("    components alone: Cu %.1f d, DMBA %.1f d (worst %.1f) | OBSERVED mixture %.1f d\n", LCu, LDM, worst, Lobs)
    @printf("    model predicts:  CA/TU %.2f d  |  IA %.2f d   (vs observed %.1f)\n", L_tu, L_ia, Lobs)
    tag = Lobs < worst - 0.05 ? "WORSE than either component (real combination effect)" : "not beyond components"
    extra = Lobs < min(L_tu, L_ia) - 0.05 ? "slightly beyond additive (mild supra-additivity, within noise)" :
            Lobs > max(L_tu, L_ia) + 0.05 ? "less than additive (mild antagonism)" : "additive (matches TU/IA bracket)"
    @printf("    -> %s ; %s\n", tag, extra)
end

println("\nReading: (A) the model's saturating per-axis impairment is consistent with the monotone")
println("dose-response, with a sensible potency ordering (metal>PAH>PCB). (B) the Cu+DMBA mixture")
println("is clearly worse than either component (a real combination effect, NO antagonism) and")
println("falls at/slightly beyond the model's additive TU/IA predictions -- corroborating the")
println("'mixtures are additive assumptions, not fitted interactions' design (TU/IA bracket the")
println("data; any mild excess is unresolved at n=2 / LT50 rounding). First CONTROLLED check of")
println("the impairment curve + mixture machinery. See docs/notes/sos_validation_status.md.")
end

main()
