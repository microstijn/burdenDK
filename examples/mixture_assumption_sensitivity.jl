# ===========================================================================
# Mixture-assumption sensitivity (with deliberate same-axis overlap)
#
# READ-ONLY. The large calibrated demo never has more than one compound on an
# axis, so TU / IA / grouped-CA-then-IA collapse to identical numbers and the
# mixture layer is effectively untested. This builds the missing overlap case:
# several compounds loading the SAME axis, and quantifies
#   (1) that the mixture models genuinely differ once there is overlap, and
#   (2) how much the margin-state output (Q_t, A_t, F) depends on the mixture
#       assumption -- i.e. how much it matters that the assumption is a proxy.
#
#   julia +release --project=. examples/mixture_assumption_sensitivity.jl
# ===========================================================================

using TwoTimescaleResilience
using Printf

# k compounds, each loading a single axis with active stress x, all others 0.
function axis_overlap(axis::Symbol, k::Int, x::Float64; distinct_codes::Bool = false)
    f(a) = a == axis ? x : 0.0
    return [(
        chemical_name = "C$i",
        effect_code = distinct_codes ? "MOR$i" : "MOR",   # same vs distinct effect codes
        burden_assimilation = f(:assimilation),
        burden_maintenance = f(:maintenance),
        burden_growth = f(:growth),
        burden_reproduction = f(:reproduction),
    ) for i in 1:k]
end

maintenance_overlap(k, x; kw...) = axis_overlap(:maintenance, k, x; kw...)
E_M(burdens, model) = aggregate_axis_mixture_effects(burdens; mixture_effect_model = model).E_maintenance
E_axis(burdens, model, axis) = getproperty(
    aggregate_axis_mixture_effects(burdens; mixture_effect_model = model),
    Symbol("E_", axis))

function main()
    println("=== (1) Do the mixture models differ under same-axis overlap? ===")
    println("    k compounds each at maintenance stress x; reporting E_maintenance.\n")
    @printf("%4s %6s %8s %8s %10s %10s %10s\n",
            "k", "x", "TU", "IA", "grp(same)", "grp(dist)", "|TU-IA|")
    for k in (2, 3, 5), x in (0.25, 0.5, 1.0, 2.0)
        b_same = maintenance_overlap(k, x; distinct_codes = false)
        b_dist = maintenance_overlap(k, x; distinct_codes = true)
        tu  = E_M(b_same, "axis_toxic_unit_sum")
        ia  = E_M(b_same, "independent_action_axis_effects")
        gs  = E_M(b_same, "grouped_ca_then_ia_axis_effects")
        gd  = E_M(b_dist, "grouped_ca_then_ia_axis_effects")
        @printf("%4d %6.2f %8.4f %8.4f %10.4f %10.4f %10.4f\n", k, x, tu, ia, gs, gd, abs(tu - ia))
    end
    println("\n  grp(same) == TU (concentration addition within one effect-code group);")
    println("  grp(dist) == IA (independent action across distinct codes). TU is always the")
    println("  most conservative (largest E); the spread |TU-IA| is the mixture-assumption")
    println("  uncertainty on the per-axis impairment.")

    println("\n=== (2) How much does the mixture choice move the margin state? ===")
    println("    The same E_a spread propagates very differently depending on the axis,")
    println("    because the kappa-rule weights are assimilation-led.\n")
    lib = load_amp_species_library()
    p = amp_species_deb_params(lib, "Daphnia_magna")
    w = axis_weights_for_species(p)
    wmap = (assimilation = w.w_assimilation, maintenance = w.w_maintenance)
    @printf("    Daphnia_magna: w_assimilation = %.3f (high), w_maintenance = %.3f (low)\n\n",
            w.w_assimilation, w.w_maintenance)

    @printf("%-13s %4s %6s | %20s | %18s\n", "overlap axis", "k", "x", "Q_t (TU / IA)", "F (TU / IA)")
    for axis in (:assimilation, :maintenance), k in (5,), x in (0.5, 1.0, 2.0)
        b = axis_overlap(axis, k, x)
        wa = getproperty(wmap, axis)
        QtTU = wa * E_axis(b, "axis_toxic_unit_sum", axis)
        QtIA = wa * E_axis(b, "independent_action_axis_effects", axis)
        FTU = amplification_from_margin(p.A0 * (1 - QtTU), p)
        FIA = amplification_from_margin(p.A0 * (1 - QtIA), p)
        @printf("%-13s %4d %6.2f | %8.4f / %8.4f | %7.4f / %7.4f\n",
                string(axis), k, x, QtTU, QtIA, FTU, FIA)
    end

    println("\nReading:")
    println("  - With same-axis overlap the mixture models genuinely differ on E_a (part 1):")
    println("    TU is the conservative upper bound, IA the lower, grouped CA-then-IA between")
    println("    them depending on shared effect codes. WITHOUT overlap (one compound per axis)")
    println("    all three are identical -- which is why the big calibrated demo shows zero")
    println("    mixture sensitivity.")
    println("  - But the *propagated* sensitivity is axis-dependent: a given E_a spread moves")
    println("    Q_t (and F) ~3x more on assimilation (w~0.5) than on maintenance (w~0.15),")
    println("    because the kappa-rule weights are assimilation-led. So the mixture-assumption")
    println("    proxy matters most for assimilation-targeting compound overlaps.")
end

main()
