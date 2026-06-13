# ===========================================================================
# Graphical-abstract crop: a wide, map-only variant of the salmon encounter map.
# Just the GeoMakie panel (path coloured by encountered adaptive margin), landscape orientation,
# suitable as a journal graphical abstract. Model + scenario shared with the other salmon figures.
#
#   julia +release --project=. examples/salmon_migration_graphical_abstract.jl
# Writes: docs/wiki/figures/salmon_migration_abstract.png
# ===========================================================================

using TwoTimescaleResilience
using CairoMakie
using GeoMakie
using Statistics

include(joinpath(@__DIR__, "salmon_migration_scenario.jl"))

const OUT_PNG = normpath(joinpath(@__DIR__, "..", "docs", "wiki", "figures", "salmon_migration_abstract.png"))

const REGION_LONLAT = Dict(
    "natal_river" => (-4.5, 56.8), "estuary" => (-1.8, 57.6),
    "coastal_shelf" => (2.5, 59.5), "open_ocean" => (-6.0, 65.5))

function build_figure()
    res = run_salmon_migration(); steps, m = res.steps, res.meta
    reg = m.regions
    wpos(s, sel) = sum(s.occ[g] * sel(reg[g]) for g in 1:length(reg))
    lons = [wpos(s, r -> REGION_LONLAT[r][1]) for s in steps]
    lats = [wpos(s, r -> REGION_LONLAT[r][2]) for s in steps]
    A_t = [s.A_t for s in steps]; Lv = [s.L for s in steps]
    med = median(A_t); dev = max(med - minimum(A_t), maximum(A_t) - med); crange = (med - dev, med + dev)
    Lmin, Lmax = extrema(Lv); ms = 13 .+ 27 .* (Lv .- Lmin) ./ (Lmax - Lmin)

    fig = Figure(size = (1180, 560))
    ax = GeoAxis(fig[1, 1];
        source = "+proj=longlat +datum=WGS84", dest = "+proj=longlat +datum=WGS84",
        limits = (-24, 11, 52, 69),
        title = "Adaptive margin along an Atlantic salmon migration")
    glon = collect(-24.0:1.0:11.0); glat = collect(52.0:0.5:69.0)
    fieldM = [1.4 * exp(-(((lo + 3) / 7)^2 + ((la - 57.5) / 5)^2)) +
              0.5 * exp(-(((lo - 2) / 6)^2 + ((la - 59) / 4)^2)) for lo in glon, la in glat]
    heatmap!(ax, glon, glat, fieldM; colormap = :amp, alpha = 0.4, colorrange = (0, maximum(fieldM)))
    lines!(ax, GeoMakie.coastlines(); color = :gray30, linewidth = 0.7, overdraw = true)
    lines!(ax, lons, lats; color = (:gray20, 0.7), linewidth = 2.5)
    scatter!(ax, lons, lats; color = A_t, colormap = :RdBu, colorrange = crange,
        markersize = ms, strokecolor = :black, strokewidth = 0.7)
    for (lab, dx, dy) in (("natal river", 0.0, -1.0), ("estuary", 2.2, -0.5),
                          ("ocean feeding", 0.0, 1.1))
        r = lab == "ocean feeding" ? "open_ocean" : replace(lab, " " => "_")
        text!(ax, REGION_LONLAT[r][1] + dx, REGION_LONLAT[r][2] + dy; text = lab,
            fontsize = 10, color = :gray20, align = (:center, :center))
    end
    Colorbar(fig[1, 2]; colormap = :RdBu, limits = crange,
        label = "adaptive margin  Aₜ  (red = eroded)", height = Relative(0.7))
    Label(fig[2, :],
        "Chronic exposure erodes capacity along the path — before any threshold is crossed; capacity is stage-resolved, exposure occupancy-weighted.";
        fontsize = 12.5, color = :gray25)
    colsize!(fig.layout, 2, Relative(0.04))
    rowgap!(fig.layout, 4)
    return fig
end

fig = build_figure()
mkpath(dirname(OUT_PNG))
save(OUT_PNG, fig; px_per_unit = 2)
println("wrote ", OUT_PNG)
