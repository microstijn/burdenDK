# ===========================================================================
# Publication figure: Atlantic salmon migration on a MAP -- solving the time dimension.
#
# The field is C(lon, lat, t) and the path is also a function of t, so a static map cannot show a
# 4-D field plus a moving point. Two complementary views resolve this:
#
#   (a) ENCOUNTER-COLLAPSED MAP (GeoMakie, real coastline, stylized North-Atlantic geography):
#       the migration path drawn over a faint time-mean synthetic exposure field, with the path
#       COLOURED BY WHAT THE FISH ENCOUNTERS each month (its adaptive margin Aₜ). Time collapses
#       onto the path because the animal is in one place per month -- no "which month's raster?".
#
#   (b) HOVMÖLLER COMPANION: the exposure field over corridor distance (y) × time (x), with the
#       fish's distance-along-corridor trajectory overlaid. This UN-collapses time: it shows the
#       field's full space-time structure and the fish threading through it (high encounter where
#       the white path crosses the bright estuary band).
#
# Geography and the field are STYLIZED/SYNTHETIC (DynQual is freshwater-only and the salmon's ocean
# phase is uncovered; the real NetCDFs are not on disk). Coastline via GeoMakie.coastlines() (no
# external files). The model + scenario are shared with the table demo and trajectory figure.
#
#   julia +release --project=. examples/salmon_migration_map.jl
# Writes: docs/wiki/figures/salmon_migration_map.png  and  docs/tex/salmon_migration_map.pdf
# ===========================================================================

using TwoTimescaleResilience
using CairoMakie
using GeoMakie
using Statistics

include(joinpath(@__DIR__, "salmon_migration_scenario.jl"))

const OUT_PNG = normpath(joinpath(@__DIR__, "..", "docs", "wiki", "figures", "salmon_migration_map.png"))
const OUT_PDF = normpath(joinpath(@__DIR__, "..", "docs", "tex", "salmon_migration_map.pdf"))

# stylized North-Atlantic geography for the four habitat zones: (lon, lat) and corridor distance d
const REGION_LONLAT = Dict(
    "natal_river"   => (-4.5, 56.8),
    "estuary"       => (-1.8, 57.6),
    "coastal_shelf" => ( 2.5, 59.5),
    "open_ocean"    => (-6.0, 65.5),
)
const REGION_D = Dict("natal_river" => 0.0, "estuary" => 0.33, "coastal_shelf" => 0.66, "open_ocean" => 1.0)

# synthetic "total exposure" along the corridor d∈[0,1] (peaks at the estuary), from region totals
function spatial_exposure(d)
    ds = [0.0, 0.33, 0.66, 1.0]
    vs = [1.4, 3.3, 0.7, 0.1]            # natal, estuary, shelf, ocean (metal+PAH)
    d <= ds[1] && return vs[1]
    d >= ds[end] && return vs[end]
    k = findlast(x -> x <= d, ds)
    f = (d - ds[k]) / (ds[k+1] - ds[k])
    return vs[k] + f * (vs[k+1] - vs[k])
end
season(t, n) = 1.0 + 0.3 * sin(2π * (t - 1) / n)        # mild seasonal modulation
field_dt(d, t, n) = spatial_exposure(d) * season(t, n)

function build_figure()
    res = run_salmon_migration()
    steps, m = res.steps, res.meta
    n  = length(steps); mo = collect(1:n)
    labels = [s.label for s in steps]
    A_t = [s.A_t for s in steps]; Lv = [s.L for s in steps]
    reg = m.regions
    wpos(s, sel) = sum(s.occ[g] * sel(reg[g]) for g in 1:length(reg))   # occupancy-weighted
    lons = [wpos(s, r -> REGION_LONLAT[r][1]) for s in steps]
    lats = [wpos(s, r -> REGION_LONLAT[r][2]) for s in steps]
    dist = [wpos(s, r -> REGION_D[r]) for s in steps]

    med = median(A_t); dev = max(med - minimum(A_t), maximum(A_t) - med); crange = (med - dev, med + dev)
    Lmin, Lmax = extrema(Lv); ms = 12 .+ 22 .* (Lv .- Lmin) ./ (Lmax - Lmin)

    fig = Figure(size = (940, 1040))

    # ---- (a) encounter-collapsed GeoMakie map ---------------------------------------------
    axm = GeoAxis(fig[1, 1];
        source = "+proj=longlat +datum=WGS84", dest = "+proj=longlat +datum=WGS84",
        limits = (-25, 12, 52, 70),
        title = "(a)  migration path over the North Atlantic — colour = encountered margin Aₜ",
        xlabel = "longitude", ylabel = "latitude")
    glon = collect(-25.0:1.0:12.0); glat = collect(52.0:0.5:70.0)
    fieldM = [1.4 * exp(-(((lo + 3) / 7)^2 + ((la - 57.5) / 5)^2)) +
              0.5 * exp(-(((lo - 2) / 6)^2 + ((la - 59) / 4)^2)) for lo in glon, la in glat]
    heatmap!(axm, glon, glat, fieldM; colormap = :amp, alpha = 0.5, colorrange = (0, maximum(fieldM)))
    lines!(axm, GeoMakie.coastlines(); color = :gray30, linewidth = 0.7, overdraw = true)
    lines!(axm, lons, lats; color = (:gray20, 0.7), linewidth = 2)
    scatter!(axm, lons, lats; color = A_t, colormap = :RdBu, colorrange = crange,
        markersize = ms, strokecolor = :black, strokewidth = 0.6)
    for (lab, dx, dy) in (("natal river", 0.0, -1.0), ("estuary", 2.0, -0.6),
                          ("coastal shelf", 0.0, 0.9), ("open ocean", 0.0, 1.1))
        r = replace(lab, " " => "_")
        text!(axm, REGION_LONLAT[r][1] + dx, REGION_LONLAT[r][2] + dy; text = lab,
            fontsize = 9, color = :gray20, align = (:center, :center))
    end
    Colorbar(fig[1, 2]; colormap = :RdBu, limits = crange,
        label = "margin Aₜ (red = eroded)", height = Relative(0.6))

    # ---- (b) Hovmöller: exposure field over corridor distance × time + fish path -----------
    axh = Axis(fig[2, 1];
        title = "(b)  exposure field along the corridor × time (Hovmöller), with the fish's path",
        xlabel = "month", ylabel = "corridor distance",
        xticks = (mo, labels), xticklabelrotation = π / 4,
        yticks = ([0.0, 0.33, 0.66, 1.0], ["natal\nriver", "estuary", "coastal\nshelf", "open\nocean"]))
    dgrid = collect(0.0:0.02:1.0)
    Fmat = [field_dt(d, t, n) for t in mo, d in dgrid]      # [month × distance]
    hmf = heatmap!(axh, mo, dgrid, Fmat; colormap = :matter)
    lines!(axh, mo, dist; color = :white, linewidth = 3)
    scatter!(axh, mo, dist; color = :white, strokecolor = :black, strokewidth = 1, markersize = 9)
    text!(axh, 1.2, 0.05; text = "white line = where the fish is", align = (:left, :bottom),
        fontsize = 9, color = :white)
    xlims!(axh, 0.5, n + 0.5); ylims!(axh, 0, 1)
    Colorbar(fig[2, 2], hmf; label = "exposure field", height = Relative(0.85))

    Label(fig[0, :],
        "Atlantic salmon migration: encounter-collapsed map (a) + space–time field (b)";
        fontsize = 15, font = :bold)
    rowsize!(fig.layout, 1, Relative(0.56))
    colsize!(fig.layout, 2, Relative(0.04))
    rowgap!(fig.layout, 8)
    return fig
end

fig = build_figure()
mkpath(dirname(OUT_PNG)); mkpath(dirname(OUT_PDF))
save(OUT_PNG, fig; px_per_unit = 2)
save(OUT_PDF, fig)
println("wrote ", OUT_PNG)
println("wrote ", OUT_PDF)
