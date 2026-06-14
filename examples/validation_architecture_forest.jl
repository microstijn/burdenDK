# validation_architecture_forest.jl — the external-validation forest plot AS A VALIDATION MAP.
# Left column = the model's processing chain (capacity -> margin state -> recovery -> function ->
# dynamics -> amplification), each stage with its formula and validation STATUS; right = the forest,
# row-aligned and colour-linked, so the reader sees *which part of the model each anchor validates*.
# Values from docs/notes/external_validation_synthesis.md. CIs Fisher-z (illustrative; small-n/digitised).
#
#   julia +release --project=. examples/validation_architecture_forest.jl
# Writes docs/tex/validation_architecture.pdf + docs/wiki/figures/validation_architecture.png

using CairoMakie

# (label, rho, n, stage, significance)  -- listed top->bottom in model-processing order
entries = [
    (L"\text{GlobTherm: capacity axis coherent}",                         0.45, 664, :cap,  "**"),
    (L"\text{GlobTherm: thermal breadth}\leftarrow k_M\ \text{(gen.\ resilience refuted)}", -0.45, 664, :cap, "**"),
    (L"\text{SFG Widdows 1995 (estuary)}",                                0.41, 36,  :state, "*"),
    (L"\text{SFG Widdows 2002 (basin)}",                                  0.12, 23,  :state, ""),
    (L"\text{SFG Albentosa 2012 (condition-confounded)}",               -0.11, 39,  :state, ""),
    (L"\text{COMADRE: }R_i\to\text{ compensation}",                       0.77, 197, :rec,  "**"),
    (L"\text{COMADRE: }k_M\to\text{ recovery (dated PGLS, rank)}",        0.221,178, :rec,  "*"),
    (L"\text{SoS DOME survival }(\,|\ \text{size}+\text{condition})",     0.45, 16,  :func, ""),
    (L"\text{Transplant: Cd-alone margin}\to\text{survival}",             0.90, 5,   :dyn,  ""),
    (L"\text{Phenanthrene mussel: margin}\to\text{survival}",             0.99, 12,  :dyn,  "**"),
    (L"\text{Phenanthrene clam: margin}\to\text{survival}",               0.97, 12,  :dyn,  "**"),
    (L"\text{Amplification }g/\mathcal{F}",                              -0.11, 197, :ampl, ""),
    (L"k_M\to\text{ acute sensitivity (raw)}",                           -0.27, 310, :bound, "**"),
    (L"k_M\to\text{ acute sensitivity }(\,|\ \text{body size})",         -0.03, 310, :bound, ""),
]

stage_color = Dict(:cap=>RGBf(.30,.60,.30), :state=>RGBf(.85,.55,.10), :rec=>RGBf(.10,.45,.70),
                   :func=>RGBf(.15,.55,.55), :dyn=>RGBf(.55,.30,.65), :ampl=>RGBf(.45,.45,.45),
                   :bound=>RGBf(.75,.25,.25))

# stage display: key, title, formula, status; chain = part of the forward model flow
stage_meta = [
    (:cap,   "Capacity  (AmP/DEB)",   L"A_0,\ \ \{w_i\},\ \ \lambda_{\min}/\lambda_{\max}", true),
    (:state, "Margin state",          L"A_t = A_0(1-Q),\ \ E=\frac{x}{1+x}",                true),
    (:rec,   "Recovery rate",         L"\lambda_{\min}=k_M,\ \ \lambda_{\max}=k_M g\,;\ \ R_i", true),
    (:func,  "Margin function",       L"A_t \to \text{acute resilience}",                   true),
    (:dyn,   "Margin dynamics",       L"\text{erosion ODE},\ B_t,\ 1/\lambda",              true),
    (:ampl,  "Amplification (scalar)",L"\mathcal{F} = \lambda(A_0)/\lambda(A_t)",            true),
    (:bound, "Bounding control",      L"\text{single-trait }k_M \to \text{tox.}",           false),
]

ci(r, n) = (n <= 3 ? (r, r) : (z = atanh(clamp(r, -0.999, 0.999)); s = 1.96/sqrt(n-3); (tanh(z-s), tanh(z+s))))
N = length(entries); ys = collect(N:-1:1)
ne_ys = [ys[i] for (i, e) in enumerate(entries) if e[4] in (:ampl, :bound)]   # the not-established zone
ne_top = maximum(ne_ys) + 0.48; ne_bot = minimum(ne_ys) - 0.48; ytop = N + 1.15

fig = Figure(size = (1180, 660))
ax_a = Axis(fig[1, 1]); ax_f = Axis(fig[1, 2])
colsize!(fig.layout, 1, Relative(0.42))

# ---------------- right: forest ----------------
poly!(ax_f, Rect(-1.1, ne_bot, 2.28, ne_top - ne_bot); color = (:gray, 0.075), strokewidth = 0)
lines!(ax_f, [-1.1, 1.18], [ne_top, ne_top]; color = (:gray, 0.45), linewidth = 0.8, linestyle = :dot)
vlines!(ax_f, [0.0], color = (:gray, 0.7), linestyle = :dash)
for (i, e) in enumerate(entries)
    c = stage_color[e[4]]; r = e[2]; lo, hi = ci(r, e[3]); y = ys[i]
    lines!(ax_f, [lo, hi], [y, y], color = c, linewidth = 2.5)
    scatter!(ax_f, [r], [y], color = c, markersize = 13, strokecolor = :white, strokewidth = 1)
    e[5] != "" && text!(ax_f, max(hi, r) + 0.04, y; text = e[5], align = (:left, :center), fontsize = 13, color = c)
    text!(ax_f, -1.04, y; text = "n=$(e[3])", align = (:left, :center), fontsize = 8.5, color = (:gray, 0.85))
end
ax_f.xlabel = L"\text{Spearman }\rho\ \ (\text{model quantity}\leftrightarrow\text{independent outcome})"
ax_f.xticks = -1:0.5:1; ax_f.yticks = (ys, [e[1] for e in entries])
ax_f.yaxisposition = :right; ax_f.ygridvisible = false; ax_f.xgridvisible = false
xlims!(ax_f, -1.1, 1.18); ylims!(ax_f, 0.3, ytop)
text!(ax_f, -1.0, N + 0.45; text = "← bounds / null", align = (:left, :center), fontsize = 9, color = (:gray, 0.7))
text!(ax_f, 1.0, N + 0.45; text = "corroboration →", align = (:right, :center), fontsize = 9, color = (:gray, 0.7))

# ---------------- left: architecture ----------------
hidespines!(ax_a); hidedecorations!(ax_a)
xlims!(ax_a, 0, 1); ylims!(ax_a, 0.3, ytop)
stage_rows(s) = [ys[i] for (i, e) in enumerate(entries) if e[4] == s]
poly!(ax_a, Rect(0.0, ne_bot, 1.0, ne_top - ne_bot); color = (:gray, 0.075), strokewidth = 0)
lines!(ax_a, [0.0, 1.0], [ne_top, ne_top]; color = (:gray, 0.45), linewidth = 0.8, linestyle = :dot)
text!(ax_a, 0.5, N + 0.78; text = "scale-free claims:  relative state · rate · dynamics",
      align = (:center, :bottom), fontsize = 10.5, font = :bold, color = RGBf(.20, .35, .25))
text!(ax_a, 0.022, (ne_bot + ne_top) / 2; text = "not established", rotation = π/2,
      align = (:center, :center), fontsize = 9, font = :bold, color = (:gray, 0.5))
prev_bot = Ref{Union{Nothing,Float64}}(nothing)
for (s, title, formula, chain) in stage_meta
    rows = stage_rows(s); isempty(rows) && continue
    top = maximum(rows) + 0.46; bot = minimum(rows) - 0.46
    poly!(ax_a, Rect(0.05, bot, 0.90, top - bot); color = (stage_color[s], 0.12),
          strokecolor = stage_color[s], strokewidth = chain ? 1.6 : 1.2)
    text!(ax_a, 0.09, top - 0.10; text = title, fontsize = 12.5, font = :bold,
          align = (:left, :top), color = stage_color[s])
    text!(ax_a, 0.09, top - 0.52; text = formula, fontsize = 10, align = (:left, :top), color = (:black, 0.72))
    if chain && prev_bot[] !== nothing      # downward flow arrowhead between chain stages
        scatter!(ax_a, [0.5], [(prev_bot[] + top) / 2]; marker = :dtriangle, markersize = 11, color = (:gray, 0.55))
    end
    prev_bot[] = chain ? bot : prev_bot[]
end
mkpath("docs/wiki/figures")
save("docs/tex/validation_architecture.pdf", fig)
save("docs/wiki/figures/validation_architecture.png", fig; px_per_unit = 2)
println("wrote docs/tex/validation_architecture.pdf + docs/wiki/figures/validation_architecture.png")
