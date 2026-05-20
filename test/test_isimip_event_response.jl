using Test
using TwoTimescaleResilience

@testset "ISIMIP Event Response" begin
    deb_params = DEBAxisParams()
    moa_map = default_isimip_moa_mapping()
    deb_map = default_moa_to_deb_mapping()

    t = collect(0.0:1.0:10.0)
    n_time = length(t)
    n_vars = 7

    A_bg = 1.0
    axes_bg = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)

    # 1. Zero D gives zero event cost and zero response if y0=0.
    D_zero = zeros(n_time, n_vars)
    res_zero = simulate_isimip_deb_event_response(t, D_zero, moa_map, deb_map, A_bg, axes_bg, deb_params; y0=0.0)
    @test all(res_zero.C_event .== 0.0)
    @test all(res_zero.response.y .== 0.0)

    # 2. Positive D gives nonnegative response.
    D_pos = zeros(n_time, n_vars)
    D_pos[2:5, 1] .= 1.0 # WT event
    res_pos = simulate_isimip_deb_event_response(t, D_pos, moa_map, deb_map, A_bg, axes_bg, deb_params; y0=0.0)
    @test all(res_pos.response.y .>= 0.0)
    @test maximum(res_pos.response.y) > 0.0

    # 3. Same event under lower lambda gives larger AUC.
    A_bg_low = 0.5 # lower background margin -> lower restoring force lambda -> slower recovery -> larger area
    res_pos_low = simulate_isimip_deb_event_response(t, D_pos, moa_map, deb_map, A_bg_low, axes_bg, deb_params; y0=0.0)
    auc_normal = sum(res_pos.response.y)
    auc_low = sum(res_pos_low.response.y)
    @test auc_low > auc_normal

    # 4. Event burden to modes works with known simple matrix.
    moa_map_simple = ModeOfActionMapping(U=fill(0.1, 8, 7))
    modes = event_burdens_to_modes(D_pos, moa_map_simple)
    @test length(modes.thermal) == n_time
    @test modes.thermal[3] ≈ 0.1
    @test modes.thermal[1] ≈ 0.0

    # 5. axes_event returned by event_burdens_to_deb_axes works with event_cost_from_axes.
    axes_event = event_burdens_to_deb_axes(modes, deb_map)
    C_ev = event_cost_from_axes(axes_event, deb_params)
    @test length(C_ev) == n_time
    @test C_ev[3] > 0.0
    @test C_ev[1] == 0.0
end
