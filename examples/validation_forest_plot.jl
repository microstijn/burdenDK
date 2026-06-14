# validation_forest_plot.jl — the external-validation forest plot (manuscript figure).
# Per-anchor Spearman rho with Fisher-z 95% CIs, grouped by organisational level.
# Values from docs/notes/external_validation_synthesis.md (2026-06-13). CIs are approximate
# (Fisher-z on n; some anchors are figure-digitised / small-n — illustrative, not inferential).
#
#   julia +release --project=. examples/validation_forest_plot.jl
# Writes docs/tex/validation_forest.pdf (vector, for the manuscript) + docs/wiki/figures/validation_forest.png

using CairoMakie

# (label, rho, n, group, significance)
entries = [
    ("COMADRE: R_i → compensation",                 0.77, 197, :pop,  "**"),
    ("COMADRE: k_M → recovery (dated PGLS, rank | pace+phylo)", 0.221, 178, :pop, "*"),
    ("GlobTherm: capacity coherence",                    0.45, 664, :phys, "**"),
    ("SFG Widdows 1995 (estuary)",                       0.41, 36,  :ener, "*"),
    ("SFG Widdows 2002 (basin)",                         0.12, 23,  :ener, ""),
    ("SFG Albentosa 2012 (condition-confounded)",       -0.11, 39,  :ener, ""),
    ("SoS DOME survival (| size + condition)",           0.45, 16,  :ener, ""),
    ("Transplant: Cd-alone margin→survival",         0.90, 5,   :dyn,  ""),
    ("Phenanthrene mussel: margin→survival",         0.99, 12,  :dyn,  "**"),
    ("Phenanthrene clam: margin→survival",           0.97, 12,  :dyn,  "**"),
    ("k_M → acute sensitivity (raw)",               -0.27, 310, :neg,  "**"),
    ("k_M → acute sensitivity (| body size)",       -0.03, 310, :neg,  ""),
    ("amplification g / F (COMADRE)",                    -0.11, 197, :null, ""),
]

groupcolor = Dict(:pop=>RGBf(0.10,0.45,0.70), :phys=>RGBf(0.30,0.60,0.30),
                  :ener=>RGBf(0.85,0.55,0.10), :dyn=>RGBf(0.55,0.30,0.65),
                  :neg=>RGBf(0.75,0.25,0.25), :null=>RGBf(0.45,0.45,0.45))
grouplabel = [(:pop,"Recovery rate endpoints (population)"), (:phys,"Capacity (physiology)"),
              (:ener,"Margin state & function (individual energetics)"), (:dyn,"Margin dynamics"),
              (:neg,"Bounding / negative controls"), (:null,"Amplification scalar (predicted null)")]

function ci(r, n)
    n <= 3 && return (r, r)
    z = atanh(clamp(r, -0.999, 0.999)); s = 1.96/sqrt(n-3)
    return (tanh(z - s), tanh(z + s))
end

N = length(entries); ys = collect(N:-1:1)
fig = Figure(size = (900, 600))
ax = Axis(fig[1,1];
    xlabel = "Spearman ρ  (model quantity ↔ independent outcome)",
    xticks = -1:0.25:1, xminorgridvisible = false,
    yticks = (ys, [e[1] for e in entries]), ygridvisible = false,
    title = "External validation of the adaptive-margin model")
vlines!(ax, [0.0], color = (:gray, 0.7), linestyle = :dash)
for (i, e) in enumerate(entries)
    r = e[2]; lo, hi = ci(r, e[3]); y = ys[i]; c = groupcolor[e[4]]
    lines!(ax, [lo, hi], [y, y], color = c, linewidth = 2.5)
    scatter!(ax, [r], [y], color = c, markersize = 13, strokecolor = :white, strokewidth = 1)
    e[5] != "" && text!(ax, max(hi, r) + 0.04, y; text = e[5], align = (:left, :center), fontsize = 13, color = c)
    text!(ax, -1.02, y; text = "n=$(e[3])", align = (:left, :center), fontsize = 9, color = (:gray, 0.8))
end
xlims!(ax, -1.08, 1.12); ylims!(ax, 0.3, N + 0.7)
elems = [MarkerElement(color = groupcolor[g], marker = :circle, markersize = 12) for (g, _) in grouplabel]
Legend(fig[2,1], elems, [l for (_, l) in grouplabel]; orientation = :horizontal,
       nbanks = 2, framevisible = false, labelsize = 10, patchsize = (12,12))
rowsize!(fig.layout, 2, Auto(0.18))

mkpath("docs/wiki/figures")
save("docs/tex/validation_forest.pdf", fig)
save("docs/wiki/figures/validation_forest.png", fig; px_per_unit = 2)
println("wrote docs/tex/validation_forest.pdf + docs/wiki/figures/validation_forest.png")
