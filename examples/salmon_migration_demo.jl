# ===========================================================================
# Migration / mobile-target proof-of-concept: Atlantic salmon (Salmo salar).
#
# The resident framework exposes a target to the field at a single cell. A mobile animal instead
# experiences the residence-time-weighted field along its trajectory. This demo runs a stylized
# anadromous salmon migration through a synthetic per-region contaminant regime, and shows that:
#
#   (1) MEMORY carries burden BETWEEN regions: burden picked up in the dirty estuary/coastal legs
#       persists, so the returning adult hits the spawning grounds already eroded -- something an
#       instantaneous single-region map cannot represent.
#   (2) STAGE-RESOLVED recovery matters along the path: salmon is an abj model (strong metabolic
#       acceleration), so the small freshwater parr recovers far faster (high lambda_max = v_eff/L)
#       than the large returning adult. A stage-blind model would badly misstate the parr's
#       recovery timescale.
#
# Reuses, unchanged: occupancy_weighted_exposure (movement), the compound-memory recurrence
# B_t = rho*B + (1-rho)*K*C, the documented MoA->axis routing (moa_to_deb_axes), and the
# stage-resolved capacity deb_params_at_length. Contaminant levels are in threshold-free relative
# pressure units (concentration / reference), as in the SFG margin examples.
#
# This is an ILLUSTRATIVE proof-of-concept (like the transplant anchor), NOT a validated movement
# model: the path and the regime are stylized, and the structure is built to accept real occupancy
# data and a real stressor field later.
#
#   julia +release --project=. examples/salmon_migration_demo.jl
# ===========================================================================

using TwoTimescaleResilience
using Printf

# --- Regions and a synthetic contaminant regime (relative-pressure units) -------------------
# Two compounds with documented modes of action: a metal (-> maintenance) and a PAH/hydrocarbon
# (-> assimilation). The freshwater and estuarine legs are dirty; the open ocean is clean.
const REGIONS = ["natal_river", "estuary", "coastal_shelf", "open_ocean"]
#                              [metal, PAH]   (relative pressure, dimensionless)
const REGION_CONC = Dict(
    "natal_river"   => [0.6, 0.8],
    "estuary"       => [1.5, 1.8],
    "coastal_shelf" => [0.4, 0.3],
    "open_ocean"    => [0.05, 0.05],
)
region_concs() = [REGION_CONC[r] for r in REGIONS]

# one-hot / split occupancy helper over REGIONS
function occ(pairs...)
    w = zeros(Float64, length(REGIONS))
    for (name, share) in pairs
        w[findfirst(==(name), REGIONS)] = share
    end
    return w
end

# Metal -> maintenance (axis 2), PAH -> assimilation (axis 1); growth/reproduction unaffected.
# Columns are [metal, PAH]; rows are [assimilation, maintenance, growth, reproduction].
const MOA_W = [0.0 1.0;
               1.0 0.0;
               0.0 0.0;
               0.0 0.0]

# Compound memory: monthly retention rho, magnification K (per compound). Illustrative.
const RHO = 0.70
const KMAG = [1.0, 1.0]

# --- Stylized migration schedule (month -> phase, occupancy, structural length) -------------
# Lengths come from the species' own ontogeny (L_b, L_p, L_i); see build_schedule below.
function build_schedule(L_b, L_p, L_i)
    L_parr  = 0.5 * (L_b + L_p)        # freshwater juvenile
    L_smolt = L_p                       # downstream migrant at puberty length
    # marine growth ramp L_p -> L_i across the ocean phase
    ramp(frac) = L_p + frac * (L_i - L_p)
    return [
        (label = "Apr parr",        phase = :juvenile, occ = occ("natal_river"=>1.0),                  L = L_parr),
        (label = "May parr",        phase = :juvenile, occ = occ("natal_river"=>1.0),                  L = L_parr),
        (label = "Jun smolt",       phase = :juvenile, occ = occ("natal_river"=>0.5, "estuary"=>0.5),  L = L_smolt),
        (label = "Jul post-smolt",  phase = :juvenile, occ = occ("estuary"=>0.5, "coastal_shelf"=>0.5),L = ramp(0.10)),
        (label = "Aug coastal",     phase = :adult,    occ = occ("coastal_shelf"=>1.0),               L = ramp(0.25)),
        (label = "Sep ocean",       phase = :adult,    occ = occ("open_ocean"=>1.0),                  L = ramp(0.55)),
        (label = "Dec ocean",       phase = :adult,    occ = occ("open_ocean"=>1.0),                  L = ramp(0.80)),
        (label = "Mar ocean",       phase = :adult,    occ = occ("open_ocean"=>1.0),                  L = ramp(0.95)),
        (label = "Apr return",      phase = :adult,    occ = occ("coastal_shelf"=>0.5, "estuary"=>0.5),L = L_i),
        (label = "May estuary",     phase = :adult,    occ = occ("estuary"=>1.0),                     L = L_i),
        (label = "Jun spawn",       phase = :adult,    occ = occ("natal_river"=>1.0),                 L = L_i),
    ]
end

# axis pressures from internal burden via the documented MoA routing
axes_from_burden(B) = moa_to_deb_axes(B, MoAToDEBMapping(W = MOA_W))

# step the compound-memory recurrence one month
step_memory(B, C) = [RHO * B[j] + (1.0 - RHO) * KMAG[j] * C[j] for j in eachindex(B)]

function main()
    lib = load_amp_species_library()
    name = "Salmo_salar"
    haskey(lib, name) || error("$name not in AmP library")
    rec = lib[name]
    has_ontogeny(rec) || error("$name has no ontogeny block; regenerate the library")
    ont = rec["ontogeny"]
    L_b, L_p, L_i = Float64(ont["L_b"]), Float64(ont["L_p"]), Float64(ont["L_i"])
    sched = build_schedule(L_b, L_p, L_i)
    rc = region_concs()

    whole = amp_record_to_deb_params(rec)   # stage-blind (whole-organism) capacity

    # Warm-start internal burden at natal-river equilibrium (the fish grew up there).
    B = [analytical_initial_burden(RHO, KMAG[j], REGION_CONC["natal_river"][j], 24) for j in 1:2]

    @printf("Atlantic salmon (Salmo salar): model=%s  s_M=%.2f  L_b=%.3f L_p=%.3f L_i=%.3f cm\n",
            ont["model"], Float64(ont["s_M"]), L_b, L_p, L_i)
    println("Contaminants in relative-pressure units; metal->maintenance, PAH->assimilation.\n")
    @printf("%-13s %-9s %6s | %5s %5s | %5s %5s | %6s | %7s %7s | %6s %6s\n",
            "month", "stage", "L(cm)", "Cmet", "Cpah", "Bmet", "Bpah",
            "A_t", "1/λ_st", "1/λ_bl", "F_st", "F_bl")
    println(repeat("-", 104))

    for s in sched
        C = occupancy_weighted_exposure(rc, s.occ)          # experienced (path-weighted) exposure
        B = step_memory(B, C)                                # carry + update internal burden
        axes = axes_from_burden(B)

        stage = deb_params_at_length(rec, s.L)               # stage-resolved capacity at length L
        r_st  = compute_adaptive_margin_response(axes, stage.params)
        r_bl  = compute_adaptive_margin_response(axes, whole)   # stage-blind

        @printf("%-13s %-9s %6.2f | %5.2f %5.2f | %5.2f %5.2f | %6.1f | %7.1f %7.1f | %6.3f %6.3f\n",
                s.label, String(stage.stage), s.L, C[1], C[2], B[1], B[2],
                r_st.A_t, 1/r_st.lambda_t, 1/r_bl.lambda_t, r_st.F_t, r_bl.F_t)
    end

    # ---- headline contrasts -------------------------------------------------------------------
    println("\n--- what movement + stage add ---")

    # (1) Memory/path: returning spawner burden vs a resident of the (clean-ish) natal river.
    Bmig = [analytical_initial_burden(RHO, KMAG[j], REGION_CONC["natal_river"][j], 24) for j in 1:2]
    for s in sched
        Bmig = step_memory(Bmig, occupancy_weighted_exposure(rc, s.occ))
    end
    Bres = [analytical_initial_burden(RHO, KMAG[j], REGION_CONC["natal_river"][j], 24) for j in 1:2]
    for _ in sched
        Bres = step_memory(Bres, REGION_CONC["natal_river"])   # never leaves the natal river
    end
    @printf("(1) MEMORY: at spawning, migrant burden [metal,PAH]=[%.2f,%.2f] vs natal-river resident [%.2f,%.2f]\n",
            Bmig[1], Bmig[2], Bres[1], Bres[2])
    println("    -> the migrant arrives to spawn carrying burden from the estuarine leg; a single-region view misses it.")

    # (2) Stage: parr vs adult recovery timescale at the SAME exposure (stage-blind = adult-like).
    parr  = deb_params_at_length(rec, 0.5 * (L_b + L_p))
    adult = deb_params_for_stage(rec, :adult)
    @printf("(2) STAGE: parr recovers ~%.1fx faster than the returning adult (1/λmax: parr=%.1f vs adult=%.1f).\n",
            parr.params.lambda_max / adult.params.lambda_max,
            1 / parr.params.lambda_max, 1 / adult.params.lambda_max)
    println("    -> a stage-blind model uses the adult (whole-organism) rate everywhere, badly overstating the parr's recovery time.")
end

main()
