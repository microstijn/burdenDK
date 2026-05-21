using Test
using TwoTimescaleResilience

@testset "ModeOfAction" begin
    @test mode_names() == (:thermal, :oxygen, :osmotic, :immune, :eutrophication, :toxic, :feeding, :physical)
    @test isimip_variable_names() == (:WT, :BOD, :TDS, :FC, :Nutrient, :Chemical, :Plastic)

    # 1. Simple linear mapping produces expected modes.
    U = zeros(8, 7)
    U[1, 1] = 1.0 # WT -> thermal
    U[2, 2] = 2.0 # BOD -> oxygen
    mapping = ModeOfActionMapping(U=U, clamp_nonnegative=false, clamp_unit=false)
    vals = [1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0]
    m = mode_of_action(vals, mapping)
    @test m.thermal ≈ 1.0
    @test m.oxygen ≈ 1.0
    @test m.osmotic ≈ 0.0

    # 2. Strict upper-triangular interaction works.
    # 3. Lower triangle is ignored.
    int_mats = [zeros(7,7) for _ in 1:8]
    int_mats[1][1, 2] = 0.5 # WT x BOD -> thermal (upper triangular)
    int_mats[1][2, 1] = 99.0 # BOD x WT (lower triangular, should be ignored)
    int_mats[1][3, 3] = 99.0 # TDS x TDS (diagonal, should be ignored)

    mapping_int = ModeOfActionMapping(U=U, interactions=int_mats, clamp_nonnegative=false, clamp_unit=false)
    vals2 = [1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    m_int = mode_of_action(vals2, mapping_int)
    # thermal = U[1,1]*1.0 + int[1,2]*1.0*2.0 = 1.0 + 0.5*2.0 = 2.0
    @test m_int.thermal ≈ 2.0

    # 4. clamp_nonnegative works.
    U_neg = fill(-1.0, 8, 7)
    mapping_neg_clamp = ModeOfActionMapping(U=U_neg, clamp_nonnegative=true)
    mapping_neg_noclamp = ModeOfActionMapping(U=U_neg, clamp_nonnegative=false)
    vals_pos = fill(1.0, 7)
    @test mode_of_action(vals_pos, mapping_neg_clamp).thermal ≈ 0.0
    @test mode_of_action(vals_pos, mapping_neg_noclamp).thermal ≈ -7.0

    # 5. clamp_unit works.
    U_large = fill(2.0, 8, 7)
    mapping_unit = ModeOfActionMapping(U=U_large, clamp_unit=true)
    @test mode_of_action(vals_pos, mapping_unit).thermal ≈ 1.0

    # 6. Wrong matrix sizes throw errors.
    @test_throws ArgumentError mode_of_action([1.0, 2.0], mapping) # wrong vals length
    mapping_bad_int_len = ModeOfActionMapping(U=U, interactions=[zeros(7,7)])
    @test_throws ArgumentError mode_of_action(vals, mapping_bad_int_len) # wrong number of interaction matrices
    mapping_bad_int_size = ModeOfActionMapping(U=U, interactions=[zeros(6,6) for _ in 1:8])
    @test_throws ArgumentError mode_of_action(vals, mapping_bad_int_size) # wrong size of interaction matrices

    # 7. mode_of_action_grid preserves shape.
    layers = [fill(0.1, 2, 2) for _ in 1:7]
    grids = mode_of_action_grid(layers, mapping)
    @test size(grids.thermal) == (2, 2)
    @test grids.thermal[1,1] ≈ 0.1
    @test grids.oxygen[1,1] ≈ 0.2

    # 8. mode_of_action_grid propagates NaN.
    layers_nan = [fill(0.1, 2, 2) for _ in 1:7]
    layers_nan[1][1, 2] = NaN
    grids_nan = mode_of_action_grid(layers_nan, mapping)
    @test isnan(grids_nan.thermal[1, 2])
    @test !isnan(grids_nan.thermal[1, 1])
    @test isnan(grids_nan.oxygen[1, 2])

    # mismatched layers
    @test_throws ArgumentError mode_of_action_grid([layers[1]], mapping)
    layers_diff_size = copy(layers)
    layers_diff_size[2] = zeros(3,3)
    @test_throws ArgumentError mode_of_action_grid(layers_diff_size, mapping)
end
