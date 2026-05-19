module TwoTimescaleResilience

export load_nc_layer, normalise_layer, write_ascii_grid, read_ascii_grid, plot_scenario_comparison, plot_grid, plot_background_layers, plot_amplification_grid,
       trapezoid_auc, recovery_time_after_last_pulse, compute_metrics, TwoTimescaleSimulationResult, simulate_two_timescale, PulseStressor, rectangular_pulse, pulse_exposure_matrix, burden_matrix_exact_update, acute_mixture_burden,
       synthetic_background_layers, run_synthetic_raster_demo, compute_background_index_grid, adaptive_margin_grid, restoring_force_grid, amplification_factor_grid,
       positive_part, phi_background, adaptive_margin, restoring_force, amplification_factor, BackgroundParams, BackgroundStressor, background_index

include("background.jl")
include("grids.jl")
include("pulses.jl")
include("simulation.jl")
include("metrics.jl")
include("plotting.jl")
include("ascii.jl")
include("netcdf.jl")

end # module
