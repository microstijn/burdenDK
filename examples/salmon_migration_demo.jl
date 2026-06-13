# ===========================================================================
# Migration / mobile-target proof-of-concept: Atlantic salmon (Salmo salar).
#
# The resident framework exposes a target to the field at a single cell. A mobile animal instead
# experiences the residence-time-weighted field along its trajectory. This demo runs a stylized
# anadromous salmon migration through a synthetic per-region contaminant regime, and shows:
#
#   (1) MOVEMENT + MEMORY: burden is picked up in the dirty estuary/coastal legs and carried
#       BETWEEN regions, so the returning adult reaches the spawning grounds already eroded --
#       an instantaneous single-region map cannot represent this.
#   (2) LIFE STAGE (recovery): salmon is an abj model (strong metabolic acceleration), so the
#       small freshwater parr recovers far faster (high lambda_max = v_eff/L) than the large
#       returning adult.
#   (3) SURFACE:VOLUME aquatic toxicokinetics: a waterborne contaminant exchanges across the
#       gill/skin SURFACE and dilutes into body VOLUME, so the exchange rate scales ∝ 1/L. Both
#       uptake and elimination scale, so steady-state burden is size-independent but the RATE is
#       not: the small parr equilibrates fast (low monthly retention rho, tracks the ambient
#       water), the large adult LAGS (high rho, carries burden longer). This is a GILL/SKIN
#       mechanism -- it applies to water-breathing aquatic targets ONLY (set WATERBORNE=false for
#       a terrestrial / air-breathing / dietary-exposed species and rho stays size-independent).
#
# Reuses, unchanged: occupancy_weighted_exposure (movement), waterborne_stage_retention
# (surface:volume rho), the compound-memory recurrence B_t = rho*B + (1-rho)*K*C, the documented
# MoA->axis routing (moa_to_deb_axes), and the stage-resolved capacity deb_params_at_length.
# Contaminant levels are in threshold-free relative pressure units, as in the SFG margin examples.
#
# ILLUSTRATIVE proof-of-concept (like the transplant anchor), NOT a validated movement model.
#
#   julia +release --project=. examples/salmon_migration_demo.jl
# ===========================================================================

using TwoTimescaleResilience
using Printf

# Salmon is a water-breathing fish -> the surface:volume gill-uptake mechanism applies.
const WATERBORNE = true

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

# Compound memory: reference monthly retention rho_ref (at the reference length L_i = adult) and
# per-compound magnification K. rho becomes length-dependent via surface:volume (see header).
const RHO_REF = 0.70
const KMAG = [1.0, 1.0]

# --- Stylized migration schedule (month -> phase, occupancy, structural length) -------------
function build_schedule(L_b, L_p, L_i)
    L_parr  = 0.5 * (L_b + L_p)        # freshwater juvenile
    L_smolt = L_p                       # downstream migrant at puberty length
    ramp(frac) = L_p + frac * (L_i - L_p)   # marine growth L_p -> L_i
    return [
        (label = "Apr parr",        occ = occ("natal_river"=>1.0),                   L = L_parr),
        (label = "May parr",        occ = occ("natal_river"=>1.0),                   L = L_parr),
        (label = "Jun smolt",       occ = occ("natal_river"=>0.5, "estuary"=>0.5),   L = L_smolt),
        (label = "Jul post-smolt",  occ = occ("estuary"=>0.5, "coastal_shelf"=>0.5), L = ramp(0.10)),
        (label = "Aug coastal",     occ = occ("coastal_shelf"=>1.0),                 L = ramp(0.25)),
        (label = "Sep ocean",       occ = occ("open_ocean"=>1.0),                    L = ramp(0.55)),
        (label = "Dec ocean",       occ = occ("open_ocean"=>1.0),                    L = ramp(0.80)),
        (label = "Mar ocean",       occ = occ("open_ocean"=>1.0),                    L = ramp(0.95)),
        (label = "Apr return",      occ = occ("coastal_shelf"=>0.5, "estuary"=>0.5), L = L_i),
        (label = "May estuary",     occ = occ("estuary"=>1.0),                       L = L_i),
        (label = "Jun spawn",       occ = occ("natal_river"=>1.0),                   L = L_i),
    ]
end

# axis pressures from internal burden via the documented MoA routing
axes_from_burden(B) = moa_to_deb_axes(B, MoAToDEBMapping(W = MOA_W))

# one month of the compound-memory recurrence at length-dependent retention rho
step_memory(B, C, rho) = [rho * B[j] + (1.0 - rho) * KMAG[j] * C[j] for j in eachindex(B)]

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
    whole = amp_record_to_deb_params(rec)              # stage-blind (whole-organism) capacity
    rho_at(L) = waterborne_stage_retention(RHO_REF, L, L_i; waterborne = WATERBORNE)

    # Warm-start internal burden at natal-river equilibrium, using the PARR's fast-equilibrating
    # retention (the fish grew up there as a small juvenile).
    rho_parr = rho_at(sched[1].L)
    B = [analytical_initial_burden(rho_parr, KMAG[j], REGION_CONC["natal_river"][j], 24) for j in 1:2]
    Bres = copy(B)                                      # resident: same fish, never leaves natal river

    @printf("Atlantic salmon (Salmo salar): model=%s  s_M=%.2f  L_b=%.3f L_p=%.3f L_i=%.3f cm  waterborne=%s\n",
            ont["model"], Float64(ont["s_M"]), L_b, L_p, L_i, WATERBORNE)
    println("Relative-pressure units; metal->maintenance, PAH->assimilation; rho(L)=rho_ref^(L_i/L).\n")
    @printf("%-13s %6s %5s | %5s %5s | %5s %5s | %6s | %7s %7s | %6s\n",
            "month", "L(cm)", "rho", "Cmet", "Cpah", "Bmet", "Bpah", "A_t", "1/λ_st", "1/λ_bl", "F_st")
    println(repeat("-", 92))

    for s in sched
        rho   = rho_at(s.L)
        C     = occupancy_weighted_exposure(rc, s.occ)         # experienced (path-weighted) exposure
        B     = step_memory(B, C, rho)                          # carry + update internal burden
        Bres  = step_memory(Bres, REGION_CONC["natal_river"], rho)
        axes  = axes_from_burden(B)
        stage = deb_params_at_length(rec, s.L)                  # stage-resolved capacity at length L
        r_st  = compute_adaptive_margin_response(axes, stage.params)
        r_bl  = compute_adaptive_margin_response(axes, whole)   # stage-blind

        @printf("%-13s %6.2f %5.2f | %5.2f %5.2f | %5.2f %5.2f | %6.1f | %7.1f %7.1f | %6.3f\n",
                s.label, s.L, rho, C[1], C[2], B[1], B[2],
                r_st.A_t, 1 / r_st.lambda_t, 1 / r_bl.lambda_t, r_st.F_t)
    end

    # ---- headline contrasts -------------------------------------------------------------------
    println("\n--- what movement + stage + surface:volume add ---")
    @printf("(1) MOVEMENT/MEMORY: at spawning, migrant burden [metal,PAH]=[%.2f,%.2f] vs natal-river resident [%.2f,%.2f]\n",
            B[1], B[2], Bres[1], Bres[2])
    println("    -> the migrant arrives to spawn carrying burden from the estuarine leg; a single-region view misses it.")

    parr  = deb_params_at_length(rec, sched[1].L)
    adult = deb_params_for_stage(rec, :adult)
    @printf("(2) STAGE (recovery): parr recovers ~%.1fx faster than the returning adult (1/λmax: parr=%.1f vs adult=%.1f).\n",
            parr.params.lambda_max / adult.params.lambda_max,
            1 / parr.params.lambda_max, 1 / adult.params.lambda_max)

    @printf("(3) SURFACE:VOLUME (aquatic TK): monthly retention rises with size -- parr rho=%.2f (equilibrates fast, tracks the water) vs adult rho=%.2f (lags, carries burden between regions).\n",
            rho_at(sched[1].L), rho_at(L_i))
    println("    -> gill/skin uptake only: for a terrestrial/air-breathing/dietary species (WATERBORNE=false) rho stays size-independent.")
end

main()
