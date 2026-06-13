# ===========================================================================
# Migration / mobile-target proof-of-concept: Atlantic salmon (Salmo salar) -- TABLE.
#
# Prints the month-by-month trajectory of the combined movement + life-stage + surface:volume
# model. The model + scenario live in salmon_migration_scenario.jl (shared with the figure script
# salmon_migration_figure.jl). See that file for the full mechanism description.
#
# Shows: (1) MOVEMENT/MEMORY -- burden is carried between regions, so the returning adult reaches
# the spawning grounds already eroded; (2) LIFE STAGE -- the small parr recovers far faster than
# the adult; (3) SURFACE:VOLUME -- monthly retention rho rises with size (parr equilibrates fast,
# tracks the water; adult lags). ILLUSTRATIVE proof-of-concept, not a validated movement model.
#
#   julia +release --project=. examples/salmon_migration_demo.jl
# ===========================================================================

using TwoTimescaleResilience
using Printf

include(joinpath(@__DIR__, "salmon_migration_scenario.jl"))

function main()
    res = run_salmon_migration()
    steps, m = res.steps, res.meta

    @printf("Atlantic salmon (Salmo salar): model=%s  s_M=%.2f  L_b=%.3f L_p=%.3f L_i=%.3f cm  waterborne=%s\n",
            m.model, m.s_M, m.L_b, m.L_p, m.L_i, WATERBORNE)
    println("Relative-pressure units; metal->maintenance, PAH->assimilation; rho(L)=rho_ref^(L_i/L).\n")
    @printf("%-13s %-9s %6s %5s | %5s %5s | %5s %5s | %6s | %7s %7s | %6s\n",
            "month", "stage", "L(cm)", "rho", "Cmet", "Cpah", "Bmet", "Bpah", "A_t", "1/λ_st", "1/λ_bl", "F_st")
    println(repeat("-", 100))
    for s in steps
        @printf("%-13s %-9s %6.2f %5.2f | %5.2f %5.2f | %5.2f %5.2f | %6.1f | %7.1f %7.1f | %6.3f\n",
                s.label, String(s.stage), s.L, s.rho, s.C[1], s.C[2], s.B[1], s.B[2],
                s.A_t, 1 / s.lambda_st, 1 / s.lambda_bl, s.F_st)
    end

    last = steps[end]
    println("\n--- what movement + stage + surface:volume add ---")
    @printf("(1) MOVEMENT/MEMORY: at spawning the migrant carries burden [metal,PAH]=[%.2f,%.2f] from the estuarine leg;\n",
            last.B[1], last.B[2])
    println("    a single-region (resident) view of the spawning grounds misses it.")
    @printf("(2) STAGE (recovery): parr recovers ~%.1fx faster than the returning adult (1/λmax: parr=%.1f vs adult=%.1f).\n",
            m.parr_lambda_max / m.adult_lambda_max, 1 / m.parr_lambda_max, 1 / m.adult_lambda_max)
    @printf("(3) SURFACE:VOLUME (aquatic TK): retention rises with size -- parr rho=%.2f (fast, tracks water) vs adult rho=%.2f (lags).\n",
            m.rho_parr, m.rho_adult)
    println("    gill/skin uptake only: a terrestrial/air-breathing/dietary species (WATERBORNE=false) keeps rho size-independent.")
end

main()
