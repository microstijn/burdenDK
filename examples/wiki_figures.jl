# ===========================================================================
# Generate the figures embedded in docs/wiki/ with CairoMakie.
#
#   julia +release --project=. examples/wiki_figures.jl
#
# Writes:
#   docs/wiki/figures/restoring_force_amplification.png
#   docs/wiki/figures/two_timescale_concept.png
#   docs/wiki/figures/kappa_collapse.png
#
# These illustrate the model mechanism and the known kappa-collapse. Edit this
# script (not the images) when the figures need updating.
# ===========================================================================

using TwoTimescaleResilience
using CairoMakie
using Statistics

const FIG_DIR = normpath(joinpath(@__DIR__, "..", "docs", "wiki", "figures"))
mkpath(FIG_DIR)

# Illustrative normalized restoring-force parameters (A0 = 1).
# kappa = 0.3 -> lambda_max/lambda_min = 1/kappa = 3.33, so the curve is visibly bent.
const A0 = 1.0
const LMIN = 0.30
const LMAX = 1.00
const KA  = 0.30 * A0

lambda(A) = LMIN + (LMAX - LMIN) * max(A, 0) / (KA + max(A, 0))
Famp(A_t) = lambda(A0) / lambda(A_t)

# --------------------------------------------------------------------------
# Figure 1: restoring force lambda(A) and amplification F(A_t)
# --------------------------------------------------------------------------
function figure_restoring_force()
    A = range(0, A0; length = 400)
    Q = range(0, 0.999; length = 400)           # erosion fraction; A_t = (1-Q)*A0
    F = [Famp((1 - q) * A0) for q in Q]

    fig = Figure(size = (920, 380))

    ax1 = Axis(fig[1, 1];
        title = "Restoring force λ(A)",
        xlabel = "adaptive margin  A / A₀", ylabel = "λ(A)")
    lines!(ax1, A ./ A0, lambda.(A); color = :steelblue, linewidth = 3)
    scatter!(ax1, [1.0], [lambda(A0)]; color = :seagreen, markersize = 12)
    text!(ax1, 1.0, lambda(A0); text = "  pristine  λ(A₀)", align = (:right, :bottom))
    scatter!(ax1, [0.25], [lambda(0.25 * A0)]; color = :firebrick, markersize = 12)
    text!(ax1, 0.25, lambda(0.25 * A0); text = "eroded  ", align = (:right, :top))
    hlines!(ax1, [LMIN]; color = :gray, linestyle = :dash)
    text!(ax1, 0.02, LMIN; text = "λ_min", align = (:left, :bottom), color = :gray)

    ax2 = Axis(fig[1, 2];
        title = "Amplification  F = λ(A₀)/λ(A_t)",
        xlabel = "margin eroded  Q  (0 = pristine, 1 = spent)", ylabel = "F")
    lines!(ax2, Q, F; color = :firebrick, linewidth = 3)
    hlines!(ax2, [1.0]; color = :gray, linestyle = :dash)
    text!(ax2, 0.02, 1.0; text = "no amplification", align = (:left, :bottom), color = :gray)

    save(joinpath(FIG_DIR, "restoring_force_amplification.png"), fig; px_per_unit = 2)
    println("wrote restoring_force_amplification.png")
end

# --------------------------------------------------------------------------
# Figure 2: the two timescales (slow margin erosion + fast event response)
# --------------------------------------------------------------------------
function figure_two_timescale()
    months = range(0, 120; length = 300)
    # slow erosion of the margin from A0 toward a depleted level
    A_slow = @. A0 * (0.4 + 0.6 * exp(-months / 45))
    t_event = 90.0
    A_at_event = A0 * (0.4 + 0.6 * exp(-t_event / 45))

    fig = Figure(size = (920, 640))

    ax1 = Axis(fig[1, 1];
        title = "Slow timescale — chronic erosion of adaptive margin",
        xlabel = "time (months)", ylabel = "A(t) / A₀")
    lines!(ax1, months, A_slow ./ A0; color = :steelblue, linewidth = 3)
    vlines!(ax1, [t_event]; color = :black, linestyle = :dot)
    scatter!(ax1, [t_event], [A_at_event / A0]; color = :firebrick, markersize = 12)
    text!(ax1, t_event, A_at_event / A0; text = "  acute event hits here", align = (:left, :center))

    # fast: impulse response y(tau) = exp(-lambda * tau); burden = area = 1/lambda
    tau = range(0, 60; length = 400)
    lam_pristine = lambda(A0)
    lam_eroded   = lambda(A_at_event)
    y_pristine = exp.(-lam_pristine .* tau)
    y_eroded   = exp.(-lam_eroded .* tau)

    ax2 = Axis(fig[2, 1];
        title = "Fast timescale — same shock, two recovery rates",
        xlabel = "event time τ", ylabel = "perturbation burden y(τ)")
    band!(ax2, tau, zeros(length(tau)), y_eroded; color = (:firebrick, 0.18))
    band!(ax2, tau, zeros(length(tau)), y_pristine; color = (:seagreen, 0.25))
    lines!(ax2, tau, y_pristine; color = :seagreen, linewidth = 3,
        label = "pristine: fast recovery, small burden")
    lines!(ax2, tau, y_eroded; color = :firebrick, linewidth = 3,
        label = "eroded: slow recovery, large burden  (area ∝ 1/λ)")
    axislegend(ax2; position = :rt)
    text!(ax2, 30, 0.55; text = "F = area(eroded) / area(pristine) = λ(A₀)/λ(A_t)",
        align = (:center, :center))

    save(joinpath(FIG_DIR, "two_timescale_concept.png"), fig; px_per_unit = 2)
    println("wrote two_timescale_concept.png")
end

# --------------------------------------------------------------------------
# Figure 3: after re-anchoring lambda_min to k_M, amplification tracks the
# energy investment ratio g, not the allocation fraction kappa.
# (Pre-fix, F lay on a single kappa curve; see docs/notes/lambda_min_maintenance_rate.tex
#  and examples/amp_kappa_collapse_diagnostic.jl for that original collapse.)
# --------------------------------------------------------------------------
function figure_kappa_to_g()
    lib = load_amp_species_library()
    kappa = Float64[]; gval = Float64[]; Fhalf = Float64[]
    for (_k, rec) in lib
        haskey(rec, "auxiliary_metrics") && haskey(rec["auxiliary_metrics"], "g") || continue
        local p
        try
            p = amp_record_to_deb_params(rec)
        catch
            continue
        end
        a = p.alpha_axes
        (all(isfinite, a) && p.A0 > 0 && 0 < a[3] < 1 && 0 < p.lambda_min <= p.lambda_max) || continue
        push!(kappa, a[3])
        push!(gval, Float64(rec["auxiliary_metrics"]["g"]))
        push!(Fhalf, amplification_from_margin(0.5 * p.A0, p))   # F at 50% erosion
    end

    fig = Figure(size = (940, 420))

    ax1 = Axis(fig[1, 1];
        title = "F vs κ — no longer a single curve",
        xlabel = "DEB allocation fraction κ", ylabel = "F  (A_t = 0.5 · A₀)")
    scatter!(ax1, kappa, Fhalf; color = log10.(gval), colormap = :viridis, markersize = 4)
    text!(ax1, 0.5, maximum(Fhalf) * 0.92;
        text = "Spearman(F, κ) = −0.11\n(was −1.000 before the λ_min fix)",
        align = (:center, :center))

    ax2 = Axis(fig[1, 2];
        title = "F vs g — amplification now tracks the energy investment ratio",
        xlabel = "energy investment ratio g  (log scale)", ylabel = "F",
        xscale = log10)
    sc = scatter!(ax2, gval, Fhalf; color = log10.(gval), colormap = :viridis, markersize = 4)
    vlines!(ax2, [1.0]; color = :gray, linestyle = :dash)
    text!(ax2, 1.0, maximum(Fhalf) * 0.92;
        text = " g ≤ 1: reserve-rich,\n clamped to F = 1", align = (:left, :center), color = :gray)
    Colorbar(fig[1, 3], sc; label = "log₁₀ g")

    save(joinpath(FIG_DIR, "amplification_vs_kappa_and_g.png"), fig; px_per_unit = 2)
    println("wrote amplification_vs_kappa_and_g.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    figure_restoring_force()
    figure_two_timescale()
    figure_kappa_to_g()
    println("\nFigures written to ", FIG_DIR)
end
