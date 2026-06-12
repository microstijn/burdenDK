# ===========================================================================
# External-ish check (#2): does the model's vulnerability ordering line up with
# any recognized life-history axis?
#
# READ-ONLY. The amplification factor F now tracks the energy investment ratio g
# (see docs/notes/lambda_min_maintenance_rate.tex). g is an internal DEB number;
# this script asks whether it -- and the resulting F -- rank-correlate with
# *interpretable* AmP life-history traits (body size, pace of life), read directly
# from allStat.mat. Strong alignment is weak corroboration; orthogonality or
# inversion is informative before investing in real chronic-then-acute data.
#
#   julia +release --project=. examples/amp_g_lifehistory_check.jl
# ===========================================================================

using TwoTimescaleResilience
using MAT
using Statistics
using Printf

function ordinalrank(v::AbstractVector)
    p = sortperm(v); r = similar(p); r[p] = 1:length(v); return Float64.(r)
end
spearman(x, y) = cor(ordinalrank(x), ordinalrank(y))

# --- model side: species -> (kappa, g, F at 50% erosion) ---
function model_indices(lib)
    kappa = Dict{String,Float64}(); gval = Dict{String,Float64}(); Fhalf = Dict{String,Float64}()
    for (sp, rec) in lib
        haskey(rec, "auxiliary_metrics") && haskey(rec["auxiliary_metrics"], "g") || continue
        local p
        try
            p = amp_record_to_deb_params(rec)
        catch
            continue
        end
        a = p.alpha_axes
        (all(isfinite, a) && p.A0 > 0 && 0 < a[3] < 1 && 0 < p.lambda_min <= p.lambda_max) || continue
        kappa[sp] = a[3]
        gval[sp]  = Float64(rec["auxiliary_metrics"]["g"])
        Fhalf[sp] = amplification_from_margin(0.5 * p.A0, p)
    end
    return kappa, gval, Fhalf
end

# Life-history traits to test (allStat keys) and what axis they capture.
const TRAITS = [
    ("Ww_i", "ultimate wet weight  (body size)"),
    ("L_i",  "ultimate length      (body size)"),
    ("a_m",  "life span            (slow life)"),
    ("a_p",  "age at puberty       (slow life)"),
    ("r_B",  "von Bertalanffy rate (fast life)"),
    ("k_M",  "maintenance rate k_M (fast metab)"),
    ("p_Am", "assimilation {p_Am}  (energy intake)"),
]

function main()
    lib = load_amp_species_library()
    kappa, gval, Fhalf = model_indices(lib)

    f = matopen(joinpath(@__DIR__, "..", "data", "allStat.mat"))
    allStat = read(f, "allStat"); close(f)

    @printf("Model species with g/F: %d\n\n", length(gval))
    @printf("%-22s %6s %14s %14s\n", "life-history trait", "n", "rho(g, trait)", "rho(F@0.5, trait)")
    println("-"^60)
    for (key, label) in TRAITS
        gx = Float64[]; gy = Float64[]; fx = Float64[]; fy = Float64[]
        for sp in keys(gval)
            haskey(allStat, sp) || continue
            d = allStat[sp]
            (d isa Dict && haskey(d, key)) || continue
            val = d[key]
            (val isa Real && isfinite(val) && val > 0) || continue
            push!(gx, gval[sp]); push!(gy, Float64(val))
            push!(fx, Fhalf[sp]); push!(fy, Float64(val))
        end
        n = length(gx)
        if n > 10
            @printf("%-22s %6d %14.3f %14.3f\n", label, n, spearman(gx, gy), spearman(fx, fy))
        else
            @printf("%-22s %6d %14s %14s\n", label, n, "--", "--")
        end
    end

    # reference: g vs kappa (should be weak, confirming g != kappa)
    ks = Float64[]; gs = Float64[]
    for sp in keys(gval); push!(ks, kappa[sp]); push!(gs, gval[sp]); end
    @printf("\nreference  rho(g, kappa) = %+.3f   (weak => g carries info beyond allocation)\n",
            spearman(gs, ks))
    println("\nReading: a strong +/- rho means the model's 'vulnerable' species (high g)")
    println("coincide with one end of a known life-history axis; ~0 means the index is")
    println("orthogonal to that axis. Direction matters for whether it is biologically sane.")
end

main()
