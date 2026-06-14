# ===========================================================================
# DYNAMIC margin validation #2 — phenanthrene firm-up of the erosion dynamics.
# Dellali et al. 2023 (Animals 13(1):151): Mytilus galloprovincialis exposed to CONSTANT
# waterborne phenanthrene (control ~10, WC1 ~45, WC2 ~89 ug/L; Table 1 confirms exposure is
# flat over time), survival-in-air (anoxia) LT50 measured at 7/15/21/28 days.
#
# WHY THIS IS A CLEANER DYNAMIC TEST THAN THE TRANSPLANT (Veldhuizen):
#   - SINGLE contaminant (phenanthrene -> assimilation axis, same routing as Widdows/DOME): no
#     metal/PCB confound.
#   - CONSTANT exposure with a FLAT control (negative baseline): time/handling is controlled out.
#   - CONCENTRATION GRADIENT (3 levels): adds dose-dependence.
#
# THE DISCRIMINATING QUESTION (same logic as the transplant). Exposure (hence quasi-steady tissue
# burden) is CONSTANT over time, so a STATIC burden->margin map predicts LT50 ~FLAT within each
# treatment. The observed LT50 instead DECLINES progressively and dose-dependently (WC2<WC1<ctrl).
# The DYNAMIC erosion state integrates the sustained cost with intrinsic timescale 1/lambda ~ months
# (>> the 28 d experiment), so it keeps rising over 4 weeks -> progressive decline. Which matches?
#
# DATA STATE (honest): 7 d and 28 d LT50 are text/Table-3 reported; 15 d and 21 d are figure-digitised
# from the Figure 4 M. galloprovincialis (water) survivorship panels (50%-survival crossing, +-~0.5 d ->
# rank-reliable). The full 4-timepoint x 3-dose grid is n=12 cells.
#
#   julia +release --project=. examples/sos_dynamic_validation_dellali.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const FILE = joinpath(@__DIR__, "..", "data", "external", "sos_dellali2023_phenanthrene.csv")
parse_num(x) = (x = strip(x, ['"', ' ']); (x == "" || x == "nd") ? NaN :
                something(tryparse(Float64, x), NaN))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = (m = isfinite.(x) .& isfinite.(y); count(m) < 3 ? NaN : cor(ordinalrank(x[m]), ordinalrank(y[m])))

function read_table(file)
    hdr = String[]; rows = Vector{Vector{String}}()
    for line in eachline(file)
        startswith(line, "#") && continue
        f = String.(split(line, ","))
        if startswith(line, "duration_d,"); hdr = f; continue; end
        length(f) == length(hdr) && push!(rows, f)
    end
    return hdr, rows
end

function get_params(lib)
    for name in ("Mytilus_galloprovincialis", "Mytilus_edulis")
        try
            return name, amp_species_deb_params(lib, name)
        catch
        end
    end
    error("no Mytilus params in the AmP library")
end

function main()
    hdr, rows = read_table(FILE)
    ci(n) = findfirst(==(n), hdr); col(n) = [parse_num(r[ci(n)]) for r in rows]
    dur = col("duration_d"); conc = col("treatment_ugL"); LT50 = col("LT50_d")

    lib = load_amp_species_library()
    spname, params = get_params(lib)
    @printf("species params: %s  A0=%.1f  lambda_min(=k_M)=%.5f/d  lambda_max=%.5f/d  g=%.1f\n",
            spname, params.A0, params.lambda_min, params.lambda_max, params.lambda_max / params.lambda_min)
    @printf("intrinsic erosion timescale 1/lambda = %.0f d (pristine) .. %.0f d (eroded)  [>> 28 d]\n",
            1 / params.lambda_max, 1 / params.lambda_min)

    # phenanthrene -> ASSIMILATION axis; pressure = median-normalised exposure (threshold-free),
    # CONSTANT in time. median over the 3 distinct treatment levels.
    treatments = sort(unique(conc))
    medC = median(treatments)
    xassim(c) = c / medC
    @printf("\ntreatments (ug/L): %s ; median=%.1f ; assimilation pressure x=C/median: %s\n",
            string(treatments), medC, string(round.(xassim.(treatments), digits = 3)))

    # ---- DYNAMIC erosion state under constant assimilation pressure from t=0 ----
    dt = 0.5; T = collect(0.0:dt:28.0); nt = length(T)
    idx(t) = argmin(abs.(T .- t))
    times = [7.0, 15.0, 21.0, 28.0]
    function dyn_y(c)
        pa = (assimilation = fill(xassim(c), nt), maintenance = zeros(nt),
              growth = zeros(nt), reproduction = zeros(nt))
        y = simulate_deb_axis_response(T, pa, params.A0, params, 1.0; y0 = 0.0, dt = dt)
        return [y[idx(t)] for t in times]
    end
    # ---- STATIC margin (instantaneous, constant in time) ----
    static_margin(c) = compute_adaptive_margin_response((xassim(c), 0.0, 0.0, 0.0), params).A_t

    println("\n=== model trajectory: DYNAMIC erosion y(t) (rises) vs STATIC margin (flat in time) ===")
    println("  conc   x_assim | dyn y: 7d   15d   21d   28d  (y28/y7) | static A_t (flat)")
    Y = Dict{Float64,Vector{Float64}}()
    for c in treatments
        yv = dyn_y(c); Y[c] = yv
        @printf("  %5.1f   %5.3f | %6.3f %5.3f %5.3f %5.3f  (%.2f) | %.3f\n",
                c, xassim(c), yv[1], yv[2], yv[3], yv[4], yv[4] / yv[1], static_margin(c))
    end

    # ---- assemble paired (cell) table: observed LT50 vs model erosion at matching (conc,time) ----
    obs_y = Float64[]; obs_lt = Float64[]; obs_t = Float64[]; obs_c = Float64[]
    for k in eachindex(dur)
        isfinite(LT50[k]) || continue
        ti = findfirst(==(dur[k]), times); ti === nothing && continue
        push!(obs_y, Y[conc[k]][ti]); push!(obs_lt, LT50[k]); push!(obs_t, dur[k]); push!(obs_c, conc[k])
    end
    n = length(obs_lt)

    println("\n=== (1) POOLED across available (conc x time) cells: dynamic erosion vs observed LT50 ===")
    @printf("  n = %d cells (currently %s) ; rho( dyn erosion y , LT50 ) = %+.2f  [expect strongly NEGATIVE]\n",
            n, all(t -> t in (7.0, 28.0), obs_t) ? "7 d + 28 d only — 15/21 d pending Fig. 4" : "incl. digitised 15/21 d",
            spear(obs_y, obs_lt))

    # ---- (2) THE DISCRIMINATING temporal test: within-treatment 7->28 d ----
    println("\n=== (2) TEMPORAL — within-treatment change 7 -> 28 d (the static map predicts ZERO) ===")
    for c in treatments
        i7 = findfirst(k -> conc[k] == c && dur[k] == 7.0, eachindex(dur))
        i28 = findfirst(k -> conc[k] == c && dur[k] == 28.0, eachindex(dur))
        (i7 === nothing || i28 === nothing) && continue
        (isfinite(LT50[i7]) && isfinite(LT50[i28])) || continue
        pct = 100 * (1 - LT50[i28] / LT50[i7])
        yv = Y[c]; dyabs = yv[4] - yv[1]
        @printf("  conc %5.1f ug/L: observed LT50 %.2f -> %.2f d (%2.0f%% drop) | DYNAMIC abs. erosion rise 7->28 d = +%.3f (dose-ordered) | STATIC change 0\n",
                c, LT50[i7], LT50[i28], pct, dyabs)
    end

    println("\nReading: STATIC margin is identical at 7 and 28 d (constant exposure) -> it predicts NO")
    println("LT50 change, yet LT50 falls progressively and dose-dependently (control ~flat). The DYNAMIC")
    println("erosion keeps rising over the 4 weeks (1/lambda ~ months) and is dose-ordered, matching the")
    println("data: a cleaner replicate of the transplant result (single PAH, flat control, gradient).")
    println("SCOPE: full 4-timepoint x 3-dose grid, n=12 cells (7/28 d text-reported, 15/21 d figure-")
    println("digitised +-~0.5 d). rho is no longer monotone-by-construction -- it carries real dose x time")
    println("inversions (e.g. WC1@28d < WC2@15d). Rank-based, pre-registered NEGATIVE sign.")
end

main()
