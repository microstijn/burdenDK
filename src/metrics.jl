function trapezoid_auc(t::Vector{Float64}, y::Vector{Float64})
    n = length(t)
    auc = 0.0
    for i in 1:(n-1)
        dt = t[i+1] - t[i]
        auc += (y[i] + y[i+1]) / 2.0 * dt
    end
    return auc
end

function recovery_time_after_last_pulse(sim::TwoTimescaleSimulationResult; fraction::Float64=0.05)
    t = sim.t
    y = sim.y
    pulses = sim.pulse_stressors

    if isempty(pulses)
        return 0.0
    end

    # Find the end of the last pulse
    last_pulse_end = maximum([p.t_end for p in pulses])

    # Find max response
    ymax = maximum(y)
    threshold = ymax * fraction

    # Start looking after the last pulse ends
    start_idx = findfirst(t .>= last_pulse_end)
    if start_idx === nothing
        return 0.0 # simulation ends before pulse ends
    end

    # Find when it drops below threshold
    for i in start_idx:length(t)
        if y[i] < threshold
            return t[i] - last_pulse_end
        end
    end

    return Inf # Did not recover within simulation window
end

function compute_metrics(sim::TwoTimescaleSimulationResult)
    Mmax, idx_Mmax = findmax(sim.M)
    t_Mmax = sim.t[idx_Mmax]

    ymax, idx_ymax = findmax(sim.y)
    t_ymax = sim.t[idx_ymax]

    auc_M = trapezoid_auc(sim.t, sim.M)
    auc_y = trapezoid_auc(sim.t, sim.y)

    q = 1.0 # default used in simulate_two_timescale
    analytical_auc_y = (q / sim.lambda) * auc_M

    auc_ratio_numeric_to_analytical = analytical_auc_y > 0 ? auc_y / analytical_auc_y : 1.0

    recovery_time = recovery_time_after_last_pulse(sim)
    residual_y = sim.y[end]

    lambda_below_warning = sim.lambda <= 0.05 ? true : false

    return (
        B = sim.B,
        A = sim.A,
        lambda = sim.lambda,
        amplification_factor = sim.amplification_factor,
        Mmax = Mmax,
        t_Mmax = t_Mmax,
        ymax = ymax,
        t_ymax = t_ymax,
        auc_M = auc_M,
        auc_y = auc_y,
        analytical_auc_y = analytical_auc_y,
        auc_ratio_numeric_to_analytical = auc_ratio_numeric_to_analytical,
        recovery_time = recovery_time,
        residual_y = residual_y,
        lambda_below_warning = lambda_below_warning
    )
end
