# ===========================================================================
# DOME margin-ranking figure — the "bridge" figure (manuscript).
#
# Shows the framework doing its LICENSED job: producing a RELATIVE vulnerability
# ranking of a real monitoring network. It uses only VALIDATED machinery — the
# stress-on-stress (SoS) anchor (modelled margin ~ survival-in-air, rho = +0.39
# raw / +0.45 controlling body size + condition; examples/sos_margin_validation_dome.jl)
# — and claims NOTHING about the unvalidated spatial/water-quality coupling.
#
# It is the bridge between the two fields: a water-quality monitoring network
# (ICES DOME OSPAR CEMP, 17 UK Mytilus edulis stations) ranked by the model's
# adaptive margin, with the INDEPENDENT survival-in-air outcome (colour) shown to
# co-vary with the ranking. This is the relative use the Discussion licenses — a
# mechanistically structured ranking of sites under sustained pressure — on
# already-defended ground.
#
#   julia +release --project=. examples/dome_margin_ranking_figure.jl
# Writes docs/tex/dome_margin_ranking.pdf (vector, manuscript) +
#        docs/wiki/figures/dome_margin_ranking.png
# ===========================================================================

using TwoTimescaleResilience
using CairoMakie
using Statistics, Printf

const FILE = joinpath(@__DIR__, "..", "data", "external", "sos_dome_ukcemp.csv")

parse_num(x) = (x = strip(x, ['"', ' ']); (x == "" || x == "nd") ? NaN :
                something(tryparse(Float64, x), NaN))
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
function spear(x, y)
    m = isfinite.(x) .& isfinite.(y)
    count(m) < 4 && return (NaN, 0)
    return (cor(ordinalrank(x[m]), ordinalrank(y[m])), count(m))
end

# --- read derived station table (same parse as sos_margin_validation_dome.jl) ---
function read_table(file)
    hdr = String[]; rows = Vector{Vector{String}}()
    for line in eachline(file)
        startswith(line, "#") && continue
        f = String.(split(line, ","))
        if startswith(line, "station_code,"); hdr = f; continue; end
        length(f) == length(hdr) && push!(rows, f)
    end
    return hdr, rows
end
hdr, rows = read_table(FILE)
col(name) = (i = findfirst(==(name), hdr); [parse_num(r[i]) for r in rows])

survival = col("survt_median_d")
metals = ("CD", "CU", "HG", "PB", "ZN")
pah23  = ("NAP", "ACNLE", "ACNE", "FLE", "PA", "ANT")
C = Dict(c => col(c) for c in (metals..., pah23..., "SCB7"))
totPAH = [sum(isfinite(C[p][i]) ? C[p][i] : 0.0 for p in pah23) for i in eachindex(rows)]
keep = isfinite.(survival) .& isfinite.(totPAH) .& isfinite.(C["CD"])
names = [rows[i][findfirst(==("station_name"), hdr)] for i in eachindex(rows) if keep[i]]

sub(v) = v[keep]
survival = sub(survival); totPAH = sub(totPAH)
for k in keys(C); C[k] = sub(C[k]); end

# --- routed margin per station (identical routing to the SoS anchor) ---
# assimilation <- 2-3 ring PAH; maintenance <- metals; reproduction <- PCB; growth <- baseline.
# Pressure = median-normalised relative burden (threshold-free).
params = amp_species_deb_params(load_amp_species_library(), "Mytilus_edulis")
tu(v) = (med = median(filter(isfinite, v)); [isfinite(x) ? x / med : 1.0 for x in v])
meanrows(vs...) = [mean(t) for t in zip(vs...)]

p_assim  = tu(totPAH)
p_maint  = meanrows((tu(C[m]) for m in metals)...)
p_repro  = tu(C["SCB7"])
p_growth = fill(1.0, length(survival))

A0 = compute_adaptive_margin_response((0.0, 0.0, 0.0, 0.0), params).A_t   # pristine margin
margin = [compute_adaptive_margin_response((p_assim[i], p_maint[i], p_growth[i], p_repro[i]), params).A_t
          for i in eachindex(survival)]
rel = margin ./ A0                                                        # fraction of pristine margin retained

rM, nM = spear(margin, survival)
@printf("stations = %d ; rho(margin, survival-in-air) = %+.3f (raw)\n", nM, rM)

# --- tidy station labels (region_site) ---
pretty(s) = replace(replace(s, r"_sh\d+$" => ""), "_" => " ")
labels = pretty.(names)

# --- sort by retained margin (most resilient at top) ---
ord = sortperm(rel)                       # ascending; plotted bottom->top so largest ends on top
rel_s = rel[ord]; surv_s = survival[ord]; lab_s = labels[ord]
N = length(rel_s); ys = collect(1:N)

# ===========================================================================
fig = Figure(size = (920, 560))
ax = Axis(fig[1, 1];
    xlabel = "Modelled adaptive margin retained,  A_t / A₀",
    yticks = (ys, lab_s), ygridvisible = false, xgridvisible = false,
    title = "Adaptive-margin ranking of a real monitoring network (ICES DOME, Mytilus edulis)")

cmap = :viridis
crange = (minimum(surv_s), maximum(surv_s))
colfor(s) = get(cgrad(cmap), (s - crange[1]) / (crange[2] - crange[1] + eps()))

for i in 1:N
    c = colfor(surv_s[i])
    lines!(ax, [0.0, rel_s[i]], [ys[i], ys[i]], color = (:gray, 0.45), linewidth = 1.5)
    scatter!(ax, [rel_s[i]], [ys[i]]; color = c, markersize = 15, strokecolor = :white, strokewidth = 1)
end
xlims!(ax, 0, max(1.0, maximum(rel_s) * 1.12)); ylims!(ax, 0.4, N + 0.6)

Colorbar(fig[1, 2]; colormap = cmap, limits = crange,
         label = "Independent outcome:\nsurvival-in-air (days)")

# annotate the validated correspondence (raw + controlled from the SoS anchor)
text!(ax, 0.02, N + 0.3;
      text = @sprintf("ρ(margin, survival) = %+.2f raw,  +0.45 | size + condition  (SoS anchor, n=%d)", rM, nM),
      align = (:left, :center), fontsize = 11, color = (:black, 0.75))

Label(fig[2, 1:2],
      "Each station is ranked by the model's adaptive margin (validated machinery); colour is the independent survival-in-air outcome. " *
      "Warmer (more resilient) clusters toward larger retained margin — the relative, mechanistically structured ranking the framework licenses. " *
      "No spatial/water-quality coupling is claimed here.";
      fontsize = 9.5, color = (:gray, 0.85), justification = :left, lineheight = 1.1)
rowsize!(fig.layout, 2, Auto(0.12))

mkpath("docs/wiki/figures")
save("docs/tex/dome_margin_ranking.pdf", fig)
save("docs/wiki/figures/dome_margin_ranking.png", fig; px_per_unit = 2)
println("wrote docs/tex/dome_margin_ranking.pdf + docs/wiki/figures/dome_margin_ranking.png")
