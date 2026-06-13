# ===========================================================================
# Landing / hero figure: the framework's core claim, margin-first and spatial.
#
# A synthetic landscape of background pressure → the model's ADAPTIVE MARGIN Aₜ (the output),
# beside a threshold-exceedance map. The point: most of the plume is BELOW the regulatory limit
# (compliant) yet its adaptive margin is already eroded — exceedance ≠ vulnerability. The margin,
# not the amplification factor, is the quantity the framework is built around and that external
# validation supports (see docs/wiki/External-Validation.md).
#
#   julia +release --project=. examples/framework_margin_map.jl
# Writes: docs/wiki/figures/framework_margin_map.png  and  docs/tex/framework_margin_map.pdf
# ===========================================================================

using TwoTimescaleResilience
using CairoMakie

const OUT_PNG = normpath(joinpath(@__DIR__, "..", "docs", "wiki", "figures", "framework_margin_map.png"))
const OUT_PDF = normpath(joinpath(@__DIR__, "..", "docs", "tex", "framework_margin_map.pdf"))
const REF    = 0.6    # threshold-free pressure reference: x = C / REF
const THRESH = 2.0    # illustrative regulatory concentration limit

function pressure_field(W, H)
    C = Array{Float64}(undef, W, H)
    for i in 1:W, j in 1:H
        u = (i - 0.5) / W; v = (j - 0.5) / H
        plume   = 2.0 * exp(-(((u - 0.18) / 0.26)^2 + ((v - 0.50) / 0.34)^2))   # broad sub-threshold plume
        hotspot = 1.7 * exp(-(((u - 0.66) / 0.07)^2 + ((v - 0.36) / 0.07)^2))   # sharp local exceedance
        bg      = 0.30 + 0.35 * u                                                # mild gradient
        C[i, j] = bg + plume + hotspot
    end
    return C
end

function build_figure()
    lib = load_amp_species_library()
    params = amp_species_deb_params(lib, "Daphnia magna")
    W, H = 96, 64
    C = pressure_field(W, H)

    relmargin = similar(C)
    for i in 1:W, j in 1:H
        x = C[i, j] / REF
        r = compute_adaptive_margin_response(
            (assimilation = x, maintenance = 0.0, growth = 0.0, reproduction = 0.0), params)
        relmargin[i, j] = r.A_t / params.A0          # fraction of pristine margin remaining (1 − Q)
    end
    exceed = Float64.(C .> THRESH)

    fig = Figure(size = (1180, 460))
    noax = (xticksvisible = false, yticksvisible = false,
            xticklabelsvisible = false, yticklabelsvisible = false)

    axa = Axis(fig[1, 1]; title = "(a)  background pressure", titlealign = :left,
        aspect = DataAspect(), noax...)
    hm1 = heatmap!(axa, C; colormap = :amp)
    Colorbar(fig[2, 1], hm1; vertical = false, label = "concentration", flipaxis = false)

    axb = Axis(fig[1, 2]; title = "(b)  adaptive margin  Aₜ / A₀   ← the model output",
        titlealign = :left, aspect = DataAspect(), noax...)
    hm2 = heatmap!(axb, relmargin; colormap = :RdBu, colorrange = (minimum(relmargin), 1.0))
    Colorbar(fig[2, 2], hm2; vertical = false, label = "margin remaining (red = eroded)", flipaxis = false)

    axc = Axis(fig[1, 3]; title = "(c)  threshold exceedance", titlealign = :left,
        aspect = DataAspect(), noax...)
    heatmap!(axc, exceed; colormap = [:gray85, :firebrick], colorrange = (0, 1))
    text!(axc, 0.5W, 0.93H; text = "grey = compliant · red = over limit",
        align = (:center, :center), fontsize = 10, color = :gray20)
    Box(fig[2, 3]; color = :transparent, strokecolor = :transparent)   # spacer to align rows

    Label(fig[0, :], "Chronic pressure erodes adaptive margin before any threshold is crossed";
        fontsize = 17, font = :bold)
    Label(fig[3, :],
        "Most of the plume is below the regulatory limit (c, grey = compliant), yet its adaptive margin is already eroded (b, red): exceedance ≠ vulnerability.";
        fontsize = 11.5, color = :gray25)

    rowsize!(fig.layout, 1, Relative(0.74))
    rowgap!(fig.layout, 6)
    return fig
end

fig = build_figure()
mkpath(dirname(OUT_PNG)); mkpath(dirname(OUT_PDF))
save(OUT_PNG, fig; px_per_unit = 2)
save(OUT_PDF, fig)
println("wrote ", OUT_PNG)
println("wrote ", OUT_PDF)
