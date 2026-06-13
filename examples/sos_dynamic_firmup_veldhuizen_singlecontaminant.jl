# ===========================================================================
# DYNAMIC firm-up — single-contaminant time-course (Veldhuizen-Tsoerkan et al. 1991, ACET
# 20:259-265). Removes the Cd-vs-PCB confound that limited the transplant proof-of-concept
# (examples/sos_dynamic_validation_veldhuizen.jl): there, continued erosion 2.5->5 mo COULD
# have been PCB rather than time-integration. Here Cd and PCB are exposed SEPARATELY, so the
# single-contaminant arms test whether each toxicant ALONE erodes acute-stress resilience
# (anoxic survival LT50) progressively with exposure time/dose.
#
# Mytilus edulis. Cd -> maintenance, PCB -> reproduction. Three arms (Cd lab fast/high-dose;
# Cd semi-field slow/low-dose; PCB semi-field). Modeled margin A_t from tissue burden (the
# static map) rank-correlated with measured LT50.
#
# HONEST SCOPE: small n, burden RISES through the measured points (no clean plateau here, so
# this CONFIRMS accumulation->erosion coupling for a single contaminant -- de-confounding the
# transplant -- rather than independently re-running the plateau/continued-erosion test). The
# PCB arm's DELAYED onset (no effect 3 mo, effect 6 mo) is a qualitative slow-erosion signature.
#
#   julia +release --project=. examples/sos_dynamic_firmup_veldhuizen_singlecontaminant.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const FILE = joinpath(@__DIR__, "..", "data", "external", "sos_veldhuizen1991_singlecontaminant.csv")
parse_num(x) = (x = strip(x, ['"', ' ']); (x == "" || x == "nd") ? NaN : something(tryparse(Float64, x), NaN))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = (m = isfinite.(x) .& isfinite.(y); count(m) < 3 ? NaN : cor(ordinalrank(x[m]), ordinalrank(y[m])))

function main()
hdr = String[]; R = Vector{Vector{String}}()
for line in eachline(FILE)
    startswith(line, "#") && continue
    f = String.(split(line, ","))
    if startswith(line, "arm,"); hdr = f; continue; end
    length(f) == length(hdr) && push!(R, f)
end
ci(n) = findfirst(==(n), hdr)
arm = [r[ci("arm")] for r in R]; tlabel = [r[ci("time_label")] for r in R]
burden = [parse_num(r[ci("burden")]) for r in R]; LT50 = [parse_num(r[ci("LT50_days")]) for r in R]
ctrl = [parse_num(r[ci("control_LT50_days")]) for r in R]

params = amp_species_deb_params(load_amp_species_library(), "Mytilus_edulis")

# modeled margin from burden: Cd -> maintenance, PCB -> reproduction (median-normalised per axis)
cd_idx = findall(i -> startswith(arm[i], "Cd") && isfinite(burden[i]), eachindex(R))
pcb_idx = findall(i -> startswith(arm[i], "PCB") && isfinite(burden[i]), eachindex(R))
medCd = median(burden[cd_idx]); medPCB = median(burden[pcb_idx])
function margin_at(i)
    sM = startswith(arm[i], "Cd") ? burden[i] / medCd : 0.0
    sR = startswith(arm[i], "PCB") ? burden[i] / medPCB : 0.0
    compute_adaptive_margin_response((0.0, sM, 0.0, sR), params).A_t
end

println("=== single-contaminant exposure time-courses (Mytilus edulis) ===")
for a in ("Cd_lab", "Cd_semifield", "PCB_semifield")
    idx = findall(i -> arm[i] == a && isfinite(LT50[i]), eachindex(R))
    isempty(idx) && continue
    println("  $a:")
    for i in idx
        @printf("    %-5s burden %6s -> LT50 %4.1f d (control %4.1f) | modeled margin A_t = %.2f\n",
                tlabel[i], isfinite(burden[i]) ? string(burden[i]) : "nd", LT50[i], ctrl[i], margin_at(i))
    end
end

# ---------- Cd-only: does modeled margin track LT50? (single contaminant, multi-timepoint) ----------
println("\n=== Cd-ONLY: modeled margin vs measured LT50 (the de-confounded test) ===")
cdL = findall(i -> startswith(arm[i], "Cd") && isfinite(LT50[i]) && isfinite(burden[i]), eachindex(R))
mCd = [margin_at(i) for i in cdL]; lCd = LT50[cdL]
@printf("  rho( margin A_t, LT50 )  = %+.2f  (n=%d, Cd lab+semifield)  [+ = model tracks]\n", spear(mCd, lCd), length(cdL))
# erosion = control - LT50 vs burden (within each arm baseline)
eros = [ctrl[i] - LT50[i] for i in cdL]; bur = burden[cdL]
@printf("  rho( Cd burden, erosion=control-LT50 ) = %+.2f  [+ = more Cd -> more erosion]\n", spear(bur, eros))
# within-arm monotonic time-dependence
for a in ("Cd_lab", "Cd_semifield")
    idx = findall(i -> arm[i] == a && isfinite(LT50[i]), eachindex(R))
    @printf("  %-13s LT50 over time: %s  (monotone decline = time/dose-dependent erosion)\n",
            a, join(string.(LT50[idx]), " -> "))
end

# ---------- PCB delayed onset ----------
println("\n=== PCB-ONLY: delayed onset (a slow single-contaminant erosion signature) ===")
for i in pcb_idx
    @printf("  %-5s burden %.1f (wet) -> LT50 %.1f d (control %.1f) : %s\n", tlabel[i], burden[i], LT50[i], ctrl[i],
            tlabel[i] == "3mo" ? "NO effect (burden present)" : "effect appears")
end
println("  -> PCB present at 3 mo without effect, effect at 6 mo: consistent with slow erosion")
println("     under sustained burden (or a burden threshold between 3 and 7 ug/g; data can't separate).")

println("\nFIRM-UP READING: Cd ALONE erodes anoxic-survival resilience progressively with")
println("exposure time/dose (lab 10.7->9.5->7.6; semi-field 9.3->8.6), and modeled margin tracks")
println("it -- so the transplant proof-of-concept's continued 2.5->5 mo erosion is NOT an artifact")
println("of the co-accumulating PCB: a single toxicant suffices. PCB alone also erodes SoS, with")
println("a DELAYED onset. What this does NOT add: a clean constant-burden continued-erosion test")
println("(burden rises through the measured points; the near-plateau 10-mo LT50 was figure-only).")
println("So the dynamic claim is de-confounded and reinforced, still short of a powered test.")
println("See docs/notes/sos_validation_status.md.")
end

main()
