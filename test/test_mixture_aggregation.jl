using Test
using TwoTimescaleResilience

@testset "Mixture Aggregation: Tranche 1" begin
    # Toy input with known compound burdens
    burdens = [
        (burden_assimilation = 1.0, burden_maintenance = 2.0, burden_growth = 3.0, burden_reproduction = 4.0),
        (burden_assimilation = 0.5, burden_maintenance = 0.0, burden_growth = 1.5, burden_reproduction = 2.0)
    ]

    # additive_axis_burden equals manual summation
    res = aggregate_deb_axis_burdens(burdens; mixture_method="additive_axis_burden")

    @test res.total_burden_assimilation ≈ 1.5
    @test res.total_burden_maintenance ≈ 2.0
    @test res.total_burden_growth ≈ 4.5
    @test res.total_burden_reproduction ≈ 6.0

    # all returned numeric values are finite
    @test isfinite(res.total_burden_assimilation)
    @test isfinite(res.total_burden_maintenance)
    @test isfinite(res.total_burden_growth)
    @test isfinite(res.total_burden_reproduction)

    # axis_toxic_unit_sum equals manual summation
    res_toxic = aggregate_deb_axis_burdens(burdens; mixture_method="axis_toxic_unit_sum")
    @test res_toxic.total_burden_assimilation ≈ 1.5
    @test res_toxic.total_burden_maintenance ≈ 2.0
    @test res_toxic.total_burden_growth ≈ 4.5
    @test res_toxic.total_burden_reproduction ≈ 6.0

    # axis_toxic_unit_sum equals additive_axis_burden for the same input
    @test res_toxic == res

    # zero-burden input returns zero totals
    zero_burdens = [
        (burden_assimilation = 0.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0),
        (burden_assimilation = 0.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)
    ]
    res_zero = aggregate_deb_axis_burdens(zero_burdens; mixture_method="additive_axis_burden")
    @test res_zero.total_burden_assimilation == 0.0
    @test res_zero.total_burden_maintenance == 0.0
    @test res_zero.total_burden_growth == 0.0
    @test res_zero.total_burden_reproduction == 0.0

    # unknown mixture_method throws ArgumentError
    @test_throws ArgumentError aggregate_deb_axis_burdens(burdens; mixture_method="unknown_method")
end

@testset "Mixture Aggregation: Tranche 3 Diagnostics" begin
    # Toy case
    diag_burdens = [
        (chemical_name = "A", burden_assimilation = 2.0, burden_maintenance = 2.0, burden_growth = 0.0, burden_reproduction = 1.0),
        (chemical_name = "B", burden_assimilation = 1.0, burden_maintenance = 2.0, burden_growth = 0.0, burden_reproduction = 3.0),
        (chemical_name = "C", burden_assimilation = 0.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)
    ]

    res_diag = mixture_contribution_diagnostics(diag_burdens)

    # Assimilation: A(2.0), B(1.0), C(0.0) => 2 contributing, dominant A, max frac 2/3
    @test res_diag.n_compounds_contributing_assimilation == 2
    @test res_diag.dominant_compound_assimilation == "A"
    @test res_diag.max_single_compound_fraction_assimilation ≈ 2/3

    # Growth: zero burden
    @test res_diag.n_compounds_contributing_growth == 0
    @test res_diag.dominant_compound_growth == "none"
    @test res_diag.max_single_compound_fraction_growth == 0.0

    # Maintenance: Tie between A and B (both 2.0). First encountered ("A") wins.
    @test res_diag.n_compounds_contributing_maintenance == 2
    @test res_diag.dominant_compound_maintenance == "A"
    @test res_diag.max_single_compound_fraction_maintenance ≈ 2.0/4.0

    # Reproduction: B(3.0) is dominant, A(1.0).
    @test res_diag.n_compounds_contributing_reproduction == 2
    @test res_diag.dominant_compound_reproduction == "B"
    @test res_diag.max_single_compound_fraction_reproduction ≈ 3.0/4.0

    # Finite values
    @test isfinite(res_diag.max_single_compound_fraction_assimilation)
    @test isfinite(res_diag.max_single_compound_fraction_maintenance)
    @test isfinite(res_diag.max_single_compound_fraction_growth)
    @test isfinite(res_diag.max_single_compound_fraction_reproduction)
end

@testset "Mixture Effects: Tranche 1" begin
    # Test 1: Zero burdens
    zeros_b = [
        (burden_assimilation = 0.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0),
        (burden_assimilation = 0.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)
    ]
    res1_tu = aggregate_axis_mixture_effects(zeros_b; mixture_effect_model="axis_toxic_unit_sum")
    @test res1_tu.X_assimilation == 0.0
    @test res1_tu.E_assimilation == 0.0
    @test res1_tu.E_maintenance == 0.0

    res1_ia = aggregate_axis_mixture_effects(zeros_b; mixture_effect_model="independent_action_axis_effects")
    @test res1_ia.X_assimilation == 0.0
    @test res1_ia.E_assimilation == 0.0

    # Test 2: Single compound x = 1 gives E = 0.5
    ones_b = [
        (burden_assimilation = 1.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)
    ]
    res2_tu = aggregate_axis_mixture_effects(ones_b; mixture_effect_model="axis_toxic_unit_sum")
    @test res2_tu.X_assimilation == 1.0
    @test res2_tu.E_assimilation == 0.5

    res2_ia = aggregate_axis_mixture_effects(ones_b; mixture_effect_model="independent_action_axis_effects")
    @test res2_ia.X_assimilation == 1.0
    @test res2_ia.E_assimilation == 0.5

    # Test 3: Two compounds x1 = 1, x2 = 1
    two_ones_b = [
        (burden_assimilation = 1.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0),
        (burden_assimilation = 1.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)
    ]
    res3_tu = aggregate_axis_mixture_effects(two_ones_b; mixture_effect_model="axis_toxic_unit_sum")
    @test res3_tu.X_assimilation == 2.0
    @test res3_tu.E_assimilation ≈ (2.0 / 3.0)

    res3_ia = aggregate_axis_mixture_effects(two_ones_b; mixture_effect_model="independent_action_axis_effects")
    @test res3_ia.E_assimilation ≈ 1.0 - (1.0 - 0.5) * (1.0 - 0.5) # 0.75

    # Test 4: Two compounds x1 = 1, x2 = 4
    mixed_b = [
        (burden_assimilation = 1.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0),
        (burden_assimilation = 4.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)
    ]
    res4_tu = aggregate_axis_mixture_effects(mixed_b; mixture_effect_model="axis_toxic_unit_sum")
    @test res4_tu.X_assimilation == 5.0
    @test res4_tu.E_assimilation ≈ (5.0 / 6.0)

    res4_ia = aggregate_axis_mixture_effects(mixed_b; mixture_effect_model="independent_action_axis_effects")
    @test res4_ia.E_assimilation ≈ 1.0 - (0.5) * (0.2) # 0.9

    # Test 5: Multiple axes are handled independently
    multi_b = [
        (burden_assimilation = 1.0, burden_maintenance = 2.0, burden_growth = 3.0, burden_reproduction = 4.0)
    ]
    res5_tu = aggregate_axis_mixture_effects(multi_b; mixture_effect_model="axis_toxic_unit_sum")
    @test res5_tu.E_assimilation ≈ 1/2
    @test res5_tu.E_maintenance ≈ 2/3
    @test res5_tu.E_growth ≈ 3/4
    @test res5_tu.E_reproduction ≈ 4/5

    # Negative burden throws
    neg_b = [(burden_assimilation = -1.0, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)]
    @test_throws ArgumentError aggregate_axis_mixture_effects(neg_b; mixture_effect_model="axis_toxic_unit_sum")

    # NaN burden throws
    nan_b = [(burden_assimilation = NaN, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)]
    @test_throws ArgumentError aggregate_axis_mixture_effects(nan_b; mixture_effect_model="axis_toxic_unit_sum")

    # Inf burden throws
    inf_b = [(burden_assimilation = Inf, burden_maintenance = 0.0, burden_growth = 0.0, burden_reproduction = 0.0)]
    @test_throws ArgumentError aggregate_axis_mixture_effects(inf_b; mixture_effect_model="axis_toxic_unit_sum")

    # Unknown model throws
    @test_throws ArgumentError aggregate_axis_mixture_effects(ones_b; mixture_effect_model="unknown_model")
end
