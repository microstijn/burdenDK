using Test
using TwoTimescaleResilience

@testset "MoAToDEBMapping" begin
    # 1. Simple matrix multiplication gives expected axes.
    W = fill(0.1, 4, 8)
    W[1, 1] = 1.0 # thermal -> assimilation
    mapping = MoAToDEBMapping(W=W, clamp_nonnegative=false, clamp_unit=false)
    modes = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    axes = moa_to_deb_axes(modes, mapping)
    @test axes.assimilation ≈ 1.0
    @test axes.maintenance ≈ 0.1
    @test axes.growth ≈ 0.1
    @test axes.reproduction ≈ 0.1

    # 2. Strict upper-triangular interactions work.
    int_mats = [zeros(8,8) for _ in 1:4]
    int_mats[1][1, 2] = 0.5 # thermal x oxygen -> assimilation
    mapping_int = MoAToDEBMapping(W=W, interactions=int_mats, clamp_nonnegative=false, clamp_unit=false)
    modes_int = [1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    axes_int = moa_to_deb_axes(modes_int, mapping_int)
    # assimilation = W[1,1]*1.0 + W[1,2]*2.0 + 0.5*1.0*2.0 = 1.0 + 0.2 + 1.0 = 2.2
    @test axes_int.assimilation ≈ 2.2

    # 3. Default mapping has size 4 x 8.
    def_mapping = default_moa_to_deb_mapping()
    @test size(def_mapping.W) == (4, 8)

    # 4. Default mapping returns nonnegative axes for nonnegative inputs.
    modes_pos = fill(1.0, 8)
    axes_def = moa_to_deb_axes(modes_pos, def_mapping)
    @test axes_def.assimilation >= 0.0
    @test axes_def.maintenance >= 0.0

    # 5. Unit clamp works.
    mapping_unit = MoAToDEBMapping(W=fill(2.0, 4, 8), clamp_unit=true)
    axes_unit = moa_to_deb_axes(modes_pos, mapping_unit)
    @test axes_unit.assimilation ≈ 1.0

    # 6. Grid mapping matches scalar mapping for constant 2x2 grids.
    modes_grid = (
        thermal = fill(1.0, 2, 2),
        oxygen = fill(0.0, 2, 2),
        osmotic = fill(0.0, 2, 2),
        immune = fill(0.0, 2, 2),
        eutrophication = fill(0.0, 2, 2),
        toxic = fill(0.0, 2, 2),
        feeding = fill(0.0, 2, 2),
        physical = fill(0.0, 2, 2)
    )
    axes_grid = moa_to_deb_axes_grid(modes_grid, mapping)
    @test axes_grid.assimilation[1,1] ≈ 1.0
    @test axes_grid.maintenance[1,1] ≈ 0.1

    # 7. NaN propagation works.
    modes_grid_nan = (
        thermal = [1.0 NaN; 1.0 1.0],
        oxygen = fill(0.0, 2, 2),
        osmotic = fill(0.0, 2, 2),
        immune = fill(0.0, 2, 2),
        eutrophication = fill(0.0, 2, 2),
        toxic = fill(0.0, 2, 2),
        feeding = fill(0.0, 2, 2),
        physical = fill(0.0, 2, 2)
    )
    axes_grid_nan = moa_to_deb_axes_grid(modes_grid_nan, mapping)
    @test isnan(axes_grid_nan.assimilation[1, 2])
    @test !isnan(axes_grid_nan.assimilation[1, 1])

    # 8. Timeseries mapping returns NamedTuple of vectors compatible with event_cost_from_axes.
    modes_ts = (
        thermal = [1.0, 2.0],
        oxygen = [0.0, 0.0],
        osmotic = [0.0, 0.0],
        immune = [0.0, 0.0],
        eutrophication = [0.0, 0.0],
        toxic = [0.0, 0.0],
        feeding = [0.0, 0.0],
        physical = [0.0, 0.0]
    )
    axes_ts = moa_to_deb_axes_timeseries(modes_ts, mapping)
    @test length(axes_ts.assimilation) == 2
    @test axes_ts.assimilation[1] ≈ 1.0
    @test axes_ts.assimilation[2] ≈ 2.0
    @test axes_ts.maintenance[1] ≈ 0.1

    # Check compatibility with event_cost_from_axes
    params = DEBAxisParams()
    cost = event_cost_from_axes(axes_ts, params)
    @test length(cost) == 2
end
