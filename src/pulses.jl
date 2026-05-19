Base.@kwdef struct PulseStressor
    name::String
    amplitude::Float64
    t_start::Float64
    t_end::Float64
    k::Float64
    weight::Float64
end

function rectangular_pulse(t::Float64, amplitude::Float64, t_start::Float64, t_end::Float64)
    if t >= t_start && t <= t_end
        return amplitude
    else
        return 0.0
    end
end

function pulse_exposure_matrix(t_vec::Vector{Float64}, pulses::Vector{PulseStressor})
    n_times = length(t_vec)
    n_pulses = length(pulses)
    P = zeros(Float64, n_times, n_pulses)

    for j in 1:n_pulses
        p = pulses[j]
        for i in 1:n_times
            P[i, j] = rectangular_pulse(t_vec[i], p.amplitude, p.t_start, p.t_end)
        end
    end

    return P
end

function burden_matrix_exact_update(t_vec::Vector{Float64}, P::Matrix{Float64}, pulses::Vector{PulseStressor}, dt::Float64)
    n_times = length(t_vec)
    n_pulses = length(pulses)
    D = zeros(Float64, n_times, n_pulses)

    # D_{i+1,j} = D_{i,j} exp(-k_j dt) + P_{i,j}(1-exp(-k_j dt))
    # where P_{i,j} is assumed constant over [t_i, t_{i+1})

    for i in 1:(n_times - 1)
        for j in 1:n_pulses
            k_j = pulses[j].k
            exp_term = exp(-k_j * dt)
            D[i+1, j] = D[i, j] * exp_term + P[i, j] * (1.0 - exp_term)
        end
    end

    return D
end

function acute_mixture_burden(Drow::Vector{Float64}, pulses::Vector{PulseStressor}; interaction=nothing)
    n_pulses = length(pulses)
    w = [p.weight for p in pulses]

    M_additive = sum(w .* Drow)

    M_interaction = 0.0
    if interaction !== nothing
        for j in 1:n_pulses
            for k in (j+1):n_pulses
                M_interaction += interaction[j, k] * Drow[j] * Drow[k]
            end
        end
    end

    return M_additive + M_interaction
end
