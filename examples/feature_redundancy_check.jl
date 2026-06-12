# ===========================================================================
# M1: how much independent signal do the margin / capacity / axis features carry,
# vs the (derived) F-features?
#
# READ-ONLY. Builds a synthetic-but-realistic spatial scenario -- real AmP species
# spanning capacity A0, a spatial species-composition gradient, and a spatially +
# seasonally + axis-varying stress field -- then runs the actual feature builder
# and the redundancy-dropping standardizer, and reports:
#   (1) how correlated the F-features are with the Q / relative-margin features,
#   (2) whether the capacity-aware absolute-margin features carry signal that the
#       relative margin does not (they only can when composition varies in space),
#   (3) what drop_redundant removes.
#
#   julia +release --project=. examples/feature_redundancy_check.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics
using Printf

# ---- pick a panel of real species spanning capacity A0 ----
function pick_species(lib; n = 8)
    params = Tuple{Float64, Any}[]   # (A0, DEBAxisParams)
    for (_k, rec) in lib
        local p
        try
            p = amp_record_to_deb_params(rec)
        catch
            continue
        end
        a = p.alpha_axes
        (all(isfinite, a) && p.A0 > 0 && 0 < a[3] < 1 && 0 < p.lambda_min <= p.lambda_max) || continue
        push!(params, (p.A0, p))
    end
    sort!(params; by = first)
    idx = round.(Int, range(1, length(params); length = n))
    return [params[i][2] for i in idx]
end

# ---- build the 5-D response arrays directly from the model functions ----
function build_scenario(species; nx = 12, ny = 12, nm = 6)
    nspec = length(species)
    nmod = 1
    Q  = fill(NaN, nx, ny, nm, nspec, nmod)
    F  = fill(NaN, nx, ny, nm, nspec, nmod)
    At = fill(NaN, nx, ny, nm, nspec, nmod)
    EA = zeros(nx, ny, nm, nmod); EM = zeros(nx, ny, nm, nmod)
    EG = zeros(nx, ny, nm, nmod); ER = zeros(nx, ny, nm, nmod)
    A0vec = [p.A0 for p in species]
    weights = [axis_weights_for_species(p) for p in species]

    Eshape(x) = x / (1 + x)   # bounded per-axis impairment

    for x in 1:nx, y in 1:ny, m in 1:nm
        season = 0.5 * (1 + sin(2pi * (m - 1) / nm))
        # axis-varying stress field: maintenance grows with y, assimilation in the west,
        # a little growth/repro pressure that varies with x -> varied axis composition
        xM = 1.5 * (y / ny) * (0.5 + season)
        xA = (x <= nx ÷ 2 ? 1.2 : 0.1) * (0.5 + 0.5 * season)
        xG = 0.4 * (x / nx)
        xR = 0.3 * (1 - x / nx)
        eA = Eshape(xA); eM = Eshape(xM); eG = Eshape(xG); eR = Eshape(xR)
        EA[x, y, m, 1] = eA; EM[x, y, m, 1] = eM; EG[x, y, m, 1] = eG; ER[x, y, m, 1] = eR

        # species-composition gradient: a moving window of present species along x
        # (low-capacity species in the west, high-capacity in the east)
        center = (x - 1) / (nx - 1) * (nspec - 1) + 1
        present = [s for s in 1:nspec if abs(s - center) <= 1.5]
        for s in present
            w = weights[s]
            Q_t = w.w_assimilation * eA + w.w_maintenance * eM + w.w_growth * eG + w.w_reproduction * eR
            Q[x, y, m, s, 1]  = Q_t
            At[x, y, m, s, 1] = species[s].A0 * (1 - Q_t)
            F[x, y, m, s, 1]  = amplification_from_margin(At[x, y, m, s, 1], species[s])
        end
    end

    return (Q_t = Q, F_t = F, A_t = At, A0 = A0vec,
            E_assimilation = EA, E_maintenance = EM, E_growth = EG, E_reproduction = ER)
end

function pearson(a, b)
    mask = isfinite.(a) .& isfinite.(b)
    sum(mask) > 2 ? cor(a[mask], b[mask]) : NaN
end

function main()
    lib = load_amp_species_library()
    species = pick_species(lib; n = 8)
    @printf("Panel: %d species, A0 from %.1f to %.3g (capacity spread)\n",
            length(species), minimum(p -> p.A0, species), maximum(p -> p.A0, species))

    resp = build_scenario(species)
    feat = build_threshold_free_vulnerability_features(
        resp; mixture_model_names = ["grouped_ca_then_ia_axis_effects"], validate = false)
    fm = feat.feature_matrix; names = feat.feature_names
    col(nm) = (i = findfirst(==(nm), names); i === nothing ? fill(NaN, size(fm, 1)) : fm[:, i])

    println("\n(1) Are the F-features just the Q / margin features? (Pearson over cells)")
    @printf("   cor(mean_F,  mean_Q)                       = %+.3f\n", pearson(col("mean_F_grouped"), col("mean_Q_grouped")))
    @printf("   cor(p95_F,   p95_Q)                        = %+.3f\n", pearson(col("p95_F_grouped"), col("p95_Q_grouped")))
    @printf("   cor(mean_F,  mean_relative_margin)         = %+.3f\n", pearson(col("mean_F_grouped"), col("mean_relative_margin_remaining")))

    println("\n(2) Do the capacity-aware features carry signal the relative margin does NOT?")
    @printf("   cor(mean_log_absolute_margin, mean_relative_margin) = %+.3f\n",
            pearson(col("mean_log_absolute_margin_remaining"), col("mean_relative_margin_remaining")))
    @printf("   cor(mean_log_baseline_capacity, mean_relative_margin) = %+.3f   (~0 => independent)\n",
            pearson(col("mean_log_baseline_capacity"), col("mean_relative_margin_remaining")))
    @printf("   cor(mean_log_baseline_capacity, mean_F)               = %+.3f\n",
            pearson(col("mean_log_baseline_capacity"), col("mean_F_grouped")))

    println("\n(3) What does drop_redundant (|cor| >= 0.95) remove?")
    std = standardize_threshold_free_vulnerability_features(
        fm, names; drop_redundant = true, correlation_threshold = 0.95)
    fcount(p) = count(n -> occursin(p, n), names)
    fkept(p)  = count(n -> occursin(p, n), std.kept_feature_names)
    @printf("   features: %d input -> %d kept, %d dropped\n",
            length(names), length(std.kept_feature_names), length(std.dropped_features))
    @printf("   F-derived features:        %d -> %d kept\n", fcount("_F_"), fkept("_F_"))
    @printf("   capacity/abs-margin feats: %d -> %d kept\n", fcount("absolute_margin") + fcount("baseline_capacity"),
            fkept("absolute_margin") + fkept("baseline_capacity"))
    println("\n   dropped features and what they duplicated:")
    for d in std.dropped_features
        if d.reason == "high_correlation"
            @printf("     - %-34s (%s)\n", d.feature_name, d.detail)
        end
    end
end

main()
