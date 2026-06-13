using Test
using TwoTimescaleResilience

# Stage-resolved capacity (life-stage integration).
#
# These tests use a self-consistent SYNTHETIC AmP record (a Dict shaped like JSON.parsefile
# output) so they are deterministic and independent of regenerating data/AmP_Species_Library.json.
# Consistency baked in: v = p_Am/A0; whole-organism lambda_max = v/L_m; for std L_i = L_m and for
# abj L_i = s_M*L_m (the DEB identity), so adult-at-L_i reproduces the whole-organism lambda_max
# in both cases. Juveniles (smaller L) recover faster: lambda_max(L) = v_eff(L)/L.

function _synthetic_record(; model, s_M, L_b, L_j, L_p, L_i)
    A0    = 1000.0            # E_m
    p_Am  = 20.0             # => v = p_Am/A0 = 0.02
    p_M   = 8.0
    E_G   = 4000.0
    k_M   = 0.002            # = p_M/E_G
    kappa = 0.8
    L_m   = 2.0             # = kappa*p_Am/p_M
    v     = p_Am / A0       # 0.02
    lambda_max = v / L_m    # 0.01 (translator uses base L_m, even for abj)
    lambda_min = min(k_M, lambda_max)
    return Dict{String, Any}(
        "A0" => A0,
        "alpha_axes" => [1.0 / A0, p_M / (kappa * p_Am), kappa, 1.0 - kappa],  # [0.001,0.5,0.8,0.2]
        "lambda_bounds" => Dict{String, Any}("lambda_min" => lambda_min, "lambda_max" => lambda_max),
        "auxiliary_metrics" => Dict{String, Any}(
            "L_m" => L_m, "p_Am" => p_Am, "p_M" => p_M, "k_M" => k_M, "E_G" => E_G,
            "g" => (v / L_m) / k_M),
        "ontogeny" => Dict{String, Any}(
            "model" => model, "s_M" => s_M, "L_b" => L_b, "L_j" => L_j, "L_p" => L_p, "L_i" => L_i)
    )
end

std_rec = _synthetic_record(model = "std", s_M = 1.0, L_b = 0.05, L_j = 0.05, L_p = 1.0, L_i = 2.0)
abj_rec = _synthetic_record(model = "abj", s_M = 3.0, L_b = 0.05, L_j = 0.15, L_p = 1.0, L_i = 6.0)

@testset "AmP life-stage (stage-resolved capacity)" begin
    whole = amp_record_to_deb_params(std_rec)

    @testset "Additive: stage API needs an ontogeny block" begin
        no_ont = deepcopy(std_rec)
        delete!(no_ont, "ontogeny")
        @test amp_record_to_deb_params(no_ont) isa DEBAxisParams   # whole-organism path intact
        @test has_ontogeny(std_rec)
        @test !has_ontogeny(no_ont)
        @test_throws ArgumentError validate_ontogeny_record(no_ont)
        @test_throws ArgumentError deb_params_for_stage(no_ont, :adult)
    end

    @testset "Regression: adult-at-L_i reproduces whole-organism recovery" begin
        adult = deb_params_for_stage(std_rec, :adult)
        @test adult isa DEBStageProfile
        @test adult.stage == :adult
        @test adult.params.lambda_max ≈ whole.lambda_max
        @test adult.params.lambda_min ≈ whole.lambda_min
        @test adult.params.A0 == whole.A0
        @test adult.params.alpha_axes == whole.alpha_axes
        @test adult.s_M_effective ≈ 1.0

        # Holds for abj too, because L_i = s_M*L_m: v*s_M / (s_M*L_m) = v/L_m.
        whole_abj = amp_record_to_deb_params(abj_rec)
        adult_abj = deb_params_for_stage(abj_rec, :adult)
        @test adult_abj.params.lambda_max ≈ whole_abj.lambda_max
        @test adult_abj.s_M_effective ≈ 3.0
    end

    @testset "Monotonicity: juvenile recovers faster than adult" begin
        juv   = deb_params_for_stage(std_rec, :juvenile)
        adult = deb_params_for_stage(std_rec, :adult)
        @test juv.stage == :juvenile
        @test juv.params.lambda_max > adult.params.lambda_max
        # A0, alpha-axes, lambda_min are stage-invariant
        @test juv.params.A0 == adult.params.A0 == whole.A0
        @test juv.params.alpha_axes == whole.alpha_axes
        @test juv.params.lambda_min ≈ whole.lambda_min
    end

    @testset "Acceleration: lambda_max constant through the s_M window" begin
        # within (L_b, L_j) = (0.05, 0.15): v_eff(L)/L = v/L_b, independent of L
        p1 = deb_params_at_length(abj_rec, 0.08)
        p2 = deb_params_at_length(abj_rec, 0.12)
        @test p1.params.lambda_max ≈ p2.params.lambda_max
        @test p1.params.lambda_max ≈ 0.02 / 0.05      # v / L_b = 0.4
        @test p1.s_M_effective > 1.0
        # std species: no acceleration anywhere
        @test deb_params_at_length(std_rec, 0.08).s_M_effective ≈ 1.0
        # v_eff primitive directly
        @test v_eff_at_length(0.10; v = 0.02, L_b = 0.05, L_j = 0.15, s_M = 3.0) ≈ 0.02 * (0.10 / 0.05)
        @test v_eff_at_length(6.0;  v = 0.02, L_b = 0.05, L_j = 0.15, s_M = 3.0) ≈ 0.02 * 3.0
        @test v_eff_at_length(1.0;  v = 0.02, L_b = 0.05, L_j = 0.05, s_M = 1.0) ≈ 0.02
    end

    @testset "Repro axis role switches at puberty; weight unchanged" begin
        pre  = deb_params_at_length(std_rec, 0.5)   # < L_p = 1.0
        post = deb_params_at_length(std_rec, 1.5)   # >= L_p
        @test pre.stage == :juvenile
        @test pre.reproductive_axis_role == :maturation
        @test post.reproductive_axis_role == :reproduction
        # kappa-rule weights identical across stages (alpha-axes unchanged ⇒ same weight vector)
        w_whole = axis_weights_for_species(whole)
        w_pre   = axis_weights_for_species(pre.params)
        w_post  = axis_weights_for_species(post.params)
        @test w_pre.w_reproduction ≈ w_whole.w_reproduction
        @test w_post.w_reproduction ≈ w_whole.w_reproduction
        @test w_pre.w_assimilation ≈ w_whole.w_assimilation
    end

    @testset "Continuous length API + guards" begin
        @test deb_params_at_length(std_rec, 1.0).params.lambda_max ≈ 0.02 / 1.0
        @test deb_params_at_length(std_rec, 0.02).stage == :embryo   # L < L_b
        @test deb_params_at_length(std_rec, 0.05).stage == :juvenile # L == L_b is birth -> juvenile
        @test_throws ArgumentError deb_params_at_length(std_rec, -1.0)
        @test_throws ArgumentError deb_params_at_length(std_rec, 0.0)
        @test_throws ArgumentError deb_params_for_stage(std_rec, :embryo)   # out of scope (no feeding)
        @test_throws ArgumentError deb_params_for_stage(std_rec, :nonsense)
    end
end
