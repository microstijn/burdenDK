# ===========================================================================
# Shared Atlantic salmon migration scenario (single source of truth).
#
# Defines the stylized anadromous path, the synthetic per-region contaminant regime, and
# run_salmon_migration(), which returns the month-by-month trajectory of the combined
# movement + life-stage + surface:volume model. Consumed by:
#   - examples/salmon_migration_demo.jl   (prints a table)
#   - examples/salmon_migration_figure.jl (publication figure)
#
# Model pieces reused unchanged: occupancy_weighted_exposure (movement),
# waterborne_stage_retention (surface:volume rho(L)), the compound-memory recurrence,
# moa_to_deb_axes (documented MoA routing), deb_params_at_length (stage-resolved capacity).
# Contaminant levels are in threshold-free relative-pressure units. ILLUSTRATIVE, not validated.
# ===========================================================================

using TwoTimescaleResilience

# Salmon is a water-breathing fish -> the surface:volume gill-uptake mechanism applies.
const WATERBORNE = true

# Habitat regions (ordered freshwater -> open ocean) and a synthetic contaminant regime.
# Two compounds: a metal (-> maintenance) and a PAH (-> assimilation). Dirty fresh/estuary legs.
const REGIONS = ["natal_river", "estuary", "coastal_shelf", "open_ocean"]
const REGION_CONC = Dict(            # [metal, PAH] in relative-pressure units
    "natal_river"   => [0.6, 0.8],
    "estuary"       => [1.5, 1.8],
    "coastal_shelf" => [0.4, 0.3],
    "open_ocean"    => [0.05, 0.05],
)
region_concs() = [REGION_CONC[r] for r in REGIONS]

# occupancy (residence-time shares over REGIONS); shares passed here sum to 1
function occ(pairs...)
    w = zeros(Float64, length(REGIONS))
    for (name, share) in pairs
        w[findfirst(==(name), REGIONS)] = share
    end
    return w
end

# Metal -> maintenance (row 2), PAH -> assimilation (row 1); cols [metal, PAH].
const MOA_W = [0.0 1.0;
              1.0 0.0;
              0.0 0.0;
              0.0 0.0]

const RHO_REF = 0.70          # monthly retention at reference length L_ref = L_i (adult)
const KMAG = [1.0, 1.0]       # per-compound magnification

function build_schedule(L_b, L_p, L_i)
    L_parr  = 0.5 * (L_b + L_p)
    L_smolt = L_p
    ramp(frac) = L_p + frac * (L_i - L_p)
    return [
        (label = "Apr parr",       occ = occ("natal_river"=>1.0),                   L = L_parr),
        (label = "May parr",       occ = occ("natal_river"=>1.0),                   L = L_parr),
        (label = "Jun smolt",      occ = occ("natal_river"=>0.5, "estuary"=>0.5),   L = L_smolt),
        (label = "Jul post-smolt", occ = occ("estuary"=>0.5, "coastal_shelf"=>0.5), L = ramp(0.10)),
        (label = "Aug coastal",    occ = occ("coastal_shelf"=>1.0),                 L = ramp(0.25)),
        (label = "Sep ocean",      occ = occ("open_ocean"=>1.0),                    L = ramp(0.55)),
        (label = "Dec ocean",      occ = occ("open_ocean"=>1.0),                    L = ramp(0.80)),
        (label = "Mar ocean",      occ = occ("open_ocean"=>1.0),                    L = ramp(0.95)),
        (label = "Apr return",     occ = occ("coastal_shelf"=>0.5, "estuary"=>0.5), L = L_i),
        (label = "May estuary",    occ = occ("estuary"=>1.0),                       L = L_i),
        (label = "Jun spawn",      occ = occ("natal_river"=>1.0),                   L = L_i),
    ]
end

axes_from_burden(B) = moa_to_deb_axes(B, MoAToDEBMapping(W = MOA_W))
step_memory(B, C, rho) = [rho * B[j] + (1.0 - rho) * KMAG[j] * C[j] for j in eachindex(B)]

"""
    run_salmon_migration(; species="Salmo_salar")

Run the combined movement + life-stage + surface:volume model along the stylized path and return
`(steps, meta)`. `steps` is a vector of per-month NamedTuples
`(label, stage, L, pos, C, B, rho, A_t, lambda_st, lambda_bl, F_st, F_bl)`, where `pos` is the
occupancy-weighted habitat index (0=natal river .. 3=open ocean). `meta` carries the ontogeny and
the parr/adult recovery extremes for annotation.
"""
function run_salmon_migration(; species::AbstractString = "Salmo_salar")
    lib = load_amp_species_library()
    haskey(lib, species) || error("$species not in AmP library")
    rec = lib[species]
    has_ontogeny(rec) || error("$species has no ontogeny block; regenerate the library")
    ont = rec["ontogeny"]
    L_b, L_p, L_i = Float64(ont["L_b"]), Float64(ont["L_p"]), Float64(ont["L_i"])

    sched = build_schedule(L_b, L_p, L_i)
    rc = region_concs()
    whole = amp_record_to_deb_params(rec)
    rho_at(L) = waterborne_stage_retention(RHO_REF, L, L_i; waterborne = WATERBORNE)

    rho_parr = rho_at(sched[1].L)
    B = [analytical_initial_burden(rho_parr, KMAG[j], REGION_CONC["natal_river"][j], 24) for j in 1:2]

    steps = NamedTuple[]
    for s in sched
        rho   = rho_at(s.L)
        C     = occupancy_weighted_exposure(rc, s.occ)
        B     = step_memory(B, C, rho)
        axes  = axes_from_burden(B)
        stage = deb_params_at_length(rec, s.L)
        r_st  = compute_adaptive_margin_response(axes, stage.params)
        r_bl  = compute_adaptive_margin_response(axes, whole)
        pos   = sum(s.occ[g] * (g - 1) for g in 1:length(REGIONS))
        push!(steps, (label = s.label, stage = stage.stage, L = s.L, pos = pos,
                      occ = copy(s.occ),
                      C = C, B = copy(B), rho = rho,
                      A_t = r_st.A_t, lambda_st = r_st.lambda_t, lambda_bl = r_bl.lambda_t,
                      F_st = r_st.F_t, F_bl = r_bl.F_t))
    end

    parr  = deb_params_at_length(rec, sched[1].L)
    adult = deb_params_for_stage(rec, :adult)
    meta = (L_b = L_b, L_p = L_p, L_i = L_i,
            model = String(ont["model"]), s_M = Float64(ont["s_M"]),
            regions = REGIONS, region_labels = ["natal\nriver", "estuary", "coastal\nshelf", "open\nocean"],
            rho_parr = rho_at(sched[1].L), rho_adult = rho_at(L_i),
            parr_lambda_max = parr.params.lambda_max, adult_lambda_max = adult.params.lambda_max)
    return (steps = steps, meta = meta)
end
