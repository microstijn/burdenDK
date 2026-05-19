Base.@kwdef struct TwoTimescaleSimulationResult
    name::String
    background_stressors::Vector{BackgroundStressor}
    pulse_stressors::Vector{PulseStressor}
    B::Float64
    A::Float64
    lambda::Float64
    amplification_factor::Float64
    t::Vector{Float64}
    P::Matrix{Float64}
    D::Matrix{Float64}
    M::Vector{Float64}
    y::Vector{Float64}
end

function simulate_two_timescale(
    name::String,
    background_stressors::Vector{BackgroundStressor},
    pulse_stressors::Vector{PulseStressor},
    background_params::BackgroundParams;
    background_interaction=nothing,
    pulse_interaction=nothing,
    tmax::Float64=100.0,
    dt::Float64=0.1,
    q::Float64=1.0
)
    # Background layer
    B = background_index(background_stressors; interaction=background_interaction)
    A = adaptive_margin(B, background_params)
    lam = restoring_force(B, background_params)
    F = amplification_factor(B, background_params)

    # Pulse layer setup
    t_vec = collect(0.0:dt:tmax)
    n_times = length(t_vec)

    # Pulse exposures
    P = pulse_exposure_matrix(t_vec, pulse_stressors)

    # TKTD Burdens
    D = burden_matrix_exact_update(t_vec, P, pulse_stressors, dt)

    # Acute mixture burden
    M = zeros(Float64, n_times)
    for i in 1:n_times
        M[i] = acute_mixture_burden(D[i, :], pulse_stressors; interaction=pulse_interaction)
    end

    # Response dynamics (exact update assuming M is piecewise constant over [t_i, t_{i+1}))
    # dy/dt = -lambda * y + q * M(t)
    # Exact solution for y_{i+1} given constant M_i over dt:
    # y_{i+1} = y_i * exp(-lambda * dt) + (q * M_i / lambda) * (1 - exp(-lambda * dt))

    y = zeros(Float64, n_times)
    y[1] = 0.0 # Initial condition

    exp_term = exp(-lam * dt)
    for i in 1:(n_times - 1)
        y[i+1] = y[i] * exp_term + (q * M[i] / lam) * (1.0 - exp_term)
    end

    return TwoTimescaleSimulationResult(
        name=name,
        background_stressors=background_stressors,
        pulse_stressors=pulse_stressors,
        B=B,
        A=A,
        lambda=lam,
        amplification_factor=F,
        t=t_vec,
        P=P,
        D=D,
        M=M,
        y=y
    )
end
