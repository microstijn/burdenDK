# ===========================================================================
# Publication figure: Atlantic salmon migration -- movement + life stage + surface:volume.
#
# A path schematic (colour = adaptive margin Aₜ on a diverging scale, marker size ∝ body length,
# ecological phases annotated) linked to along-path profiles (exposure→burden memory with an
# illustrative limit + compliance strip, margin, recovery time, retention). The model + scenario
# live in salmon_migration_scenario.jl (shared with the table demo).
#
#   julia +release --project=. examples/salmon_migration_figure.jl
# Writes: docs/wiki/figures/salmon_migration.png  and  docs/tex/salmon_migration.pdf
# ===========================================================================

using TwoTimescaleResilience
using CairoMakie
using Statistics

include(joinpath(@__DIR__, "salmon_migration_scenario.jl"))

const OUT_PNG = normpath(joinpath(@__DIR__, "..", "docs", "wiki", "figures", "salmon_migration.png"))
const OUT_PDF = normpath(joinpath(@__DIR__, "..", "docs", "tex", "salmon_migration.pdf"))
const THRESH  = 1.0   # illustrative regulatory exposure limit (relative-pressure units)

function build_figure()
    res = run_salmon_migration()
    steps, m = res.steps, res.meta
    n   = length(steps);  mo = collect(1:n)
    labels = [s.label for s in steps]
    pos = [s.pos for s in steps];  A_t = [s.A_t for s in steps];  Lv = [s.L for s in steps]
    Cmet = [s.C[1] for s in steps]; Cpah = [s.C[2] for s in steps]
    Bmet = [s.B[1] for s in steps]; Bpah = [s.B[2] for s in steps]
    rho = [s.rho for s in steps]
    invst = [1 / s.lambda_st for s in steps]; invbl = [1 / s.lambda_bl for s in steps]
    isjuv = [s.stage == :juvenile for s in steps]
    compliant = [(Cmet[i] < THRESH) && (Cpah[i] < THRESH) for i in 1:n]

    # diverging margin scale centred on the trajectory median (red = more eroded)
    med = median(A_t); dev = max(med - minimum(A_t), maximum(A_t) - med)
    crange = (med - dev, med + dev)

    fig = Figure(size = (1180, 860))

    # ---- (a) journey schematic: habitat (x) vs month (y), colour = Aₜ, size ∝ length ------
    axA = Axis(fig[1:4, 1];
        title = "(a)  migration path\ncolour = margin Aₜ (red=eroded) · size ∝ length",
        titlealign = :left, xlabel = "habitat",
        xticks = (0:3, m.region_labels), yticks = (mo, labels), yreversed = true)
    xlims!(axA, -0.6, 3.6)
    vlines!(axA, [0.5, 1.5, 2.5]; color = (:gray, 0.25), linewidth = 1)
    lines!(axA, pos, mo; color = (:gray, 0.6), linewidth = 2)
    Lmin, Lmax = extrema(Lv)
    ms = 11 .+ 22 .* (Lv .- Lmin) ./ (Lmax - Lmin)
    scatter!(axA, pos, mo; color = A_t, colormap = :RdBu, colorrange = crange,
        markersize = ms, strokecolor = :black, strokewidth = 0.6)
    # ecological phase annotations (placed in open space)
    text!(axA, 2.5, 1.6;  text = "freshwater\nparr · smolt", align = (:center, :center),
        fontsize = 10, color = :gray30)
    text!(axA, 0.75, 6.4; text = "marine\nfeeding", align = (:center, :center),
        fontsize = 10, color = :gray30)
    text!(axA, 2.7, 10.0; text = "spawning\nreturn", align = (:center, :center),
        fontsize = 10, color = :gray30)
    Colorbar(fig[1:4, 2]; colormap = :RdBu, limits = crange,
        label = "adaptive margin  Aₜ  (red = more eroded)", height = Relative(0.5))

    function stage_band!(ax)   # shade the DEB-juvenile months
        for i in 1:n
            isjuv[i] && vspan!(ax, i - 0.5, i + 0.5; color = (:darkorange, 0.10))
        end
    end

    # ---- (b) exposure -> burden (memory), with limit + compliance ---------------------------
    axB = Axis(fig[1, 3]; title = "(b)  exposure (dashed) → internal burden (solid): memory",
        titlealign = :left, ylabel = "rel. units")
    stage_band!(axB)
    hlines!(axB, [THRESH]; color = (:black, 0.6), linestyle = :dot, linewidth = 1)
    text!(axB, n - 0.2, THRESH; text = "illustrative limit", align = (:right, :bottom),
        fontsize = 9, color = :gray30)
    lines!(axB, mo, Cmet; color = (:steelblue, 0.5), linestyle = :dash)
    lines!(axB, mo, Bmet; color = :steelblue, linewidth = 2.5, label = "metal")
    lines!(axB, mo, Cpah; color = (:firebrick, 0.5), linestyle = :dash)
    lines!(axB, mo, Bpah; color = :firebrick, linewidth = 2.5, label = "PAH")
    axislegend(axB; position = :rt, labelsize = 10, framevisible = false, orientation = :horizontal)
    hidexdecorations!(axB; grid = false); xlims!(axB, 0.5, n + 0.5)

    # ---- (c) margin, with exceedance-vs-vulnerability compliance strip ----------------------
    axC = Axis(fig[2, 3]; title = "(c)  adaptive margin Aₜ  +  compliance (▮ below / ▮ above limit)",
        titlealign = :left, ylabel = "Aₜ")
    stage_band!(axC)
    rngA = maximum(A_t) - minimum(A_t)
    ystrip = minimum(A_t) - 0.10 * rngA
    lines!(axC, mo, A_t; color = :black, linewidth = 2.5)
    scatter!(axC, mo, fill(ystrip, n);
        color = [compliant[i] ? :seagreen : :firebrick for i in 1:n],
        marker = :rect, markersize = 13)
    hidexdecorations!(axC; grid = false)
    xlims!(axC, 0.5, n + 0.5); ylims!(axC, ystrip - 0.06 * rngA, maximum(A_t) + 0.06 * rngA)

    # ---- (d) recovery time, stage-resolved vs stage-blind ----------------------------------
    axD = Axis(fig[3, 3]; title = "(d)  recovery time 1/λ : stage-resolved vs stage-blind",
        titlealign = :left, ylabel = "1/λ")
    stage_band!(axD)
    lines!(axD, mo, invst; color = :seagreen, linewidth = 2.5, label = "stage-resolved")
    lines!(axD, mo, invbl; color = :gray, linewidth = 2, linestyle = :dash, label = "stage-blind")
    axislegend(axD; position = :lt, labelsize = 10, framevisible = false)
    hidexdecorations!(axD; grid = false); xlims!(axD, 0.5, n + 0.5)

    # ---- (e) retention rho(L) : surface:volume ---------------------------------------------
    axE = Axis(fig[4, 3]; title = "(e)  monthly retention ρ(L) : surface:volume aquatic TK",
        titlealign = :left, ylabel = "ρ", xlabel = "month",
        xticks = (mo, labels), xticklabelrotation = π / 4)
    stage_band!(axE)
    lines!(axE, mo, rho; color = :purple, linewidth = 2.5)
    xlims!(axE, 0.5, n + 0.5); ylims!(axE, 0, 0.8)

    Label(fig[0, :],
        "Atlantic salmon: movement + life stage + surface:volume toxicokinetics  (orange band = juvenile)";
        fontsize = 16, font = :bold)

    colsize!(fig.layout, 1, Relative(0.30))
    colsize!(fig.layout, 2, Relative(0.05))
    rowgap!(fig.layout, 6); colgap!(fig.layout, 10)
    return fig
end

fig = build_figure()
mkpath(dirname(OUT_PNG)); mkpath(dirname(OUT_PDF))
save(OUT_PNG, fig; px_per_unit = 2)
save(OUT_PDF, fig)
println("wrote ", OUT_PNG)
println("wrote ", OUT_PDF)
