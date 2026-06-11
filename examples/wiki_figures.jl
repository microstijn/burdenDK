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
# Figure 3: amplification is a function of kappa alone (real library)
# --------------------------------------------------------------------------
function figure_kappa_collapse()
    lib = load_amp_species_library()
    kappa = Float64[]; Fhalf = Float64[]; logA0 = Float64[]
    for (_k, rec) in lib
        local p
        try
            p = amp_record_to_deb_params(rec)
        catch
            continue
        end
        a = p.alpha_axes
        (all(isfinite, a) && p.A0 > 0 && 0 < a[3] < 1 && 0 < p.lambda_min < p.lambda_max) || continue
        push!(kappa, a[3])
        push!(Fhalf, amplification_from_margin(0.5 * p.A0, p))   # F at 50% erosion
        push!(logA0, log10(p.A0))
    end

    fig = Figure(size = (760, 480))
    ax = Axis(fig[1, 1];
        title = "Amplification at 50% erosion is a function of κ alone",
        xlabel = "DEB allocation fraction κ", ylabel = "F  (A_t = 0.5 · A₀)")
    sc = scatter!(ax, kappa, Fhalf; color = logA0, colormap = :viridis, markersize = 5)
    Colorbar(fig[1, 2], sc; label = "log₁₀ A₀  (6 orders of magnitude — yet F doesn't move)")
    text!(ax, 0.5, maximum(Fhalf) * 0.9;
        text = "every species lies on one κ curve\n(Spearman(F, κ) = −1.000)",
        align = (:center, :center))

    save(joinpath(FIG_DIR, "kappa_collapse.png"), fig; px_per_unit = 2)
    println("wrote kappa_collapse.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    figure_restoring_force()
    figure_two_timescale()
    figure_kappa_collapse()
    println("\nFigures written to ", FIG_DIR)
end
