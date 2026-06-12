# Background parameter types and scalar functions

Base.@kwdef struct BackgroundParams
    A0::Float64 = 1.0
    alpha::Float64 = 1.15
    KB::Float64 = 0.55
    hill::Float64 = 2.0
    use_saturating_phi::Bool = true

    lambda_min::Float64 = 0.04
    lambda_max::Float64 = 1.0
end

positive_part(x::Real) = max(float(x), 0.0)

function phi_background(B::Real, params::BackgroundParams)
    if params.use_saturating_phi
        num = params.alpha * (B ^ params.hill)
        den = (params.KB ^ params.hill) + (B ^ params.hill)
        return num / den
    else
        return params.alpha * B
    end
end

function adaptive_margin(B::Real, params::BackgroundParams)
    return params.A0 - phi_background(B, params)
end

function restoring_force(B::Real, params::BackgroundParams)
    A_val = adaptive_margin(B, params)
    A_pos = positive_part(A_val)

    # Linear recovery curve (see deb_axes.jl): the restoring force scales with the
    # fraction of pristine margin retained, between the two rate bounds.
    #   lambda(B) = lambda_min + (lambda_max - lambda_min) * clamp(A(B)/A0, 0, 1)
    diff = params.lambda_max - params.lambda_min
    return params.lambda_min + diff * clamp(A_pos / params.A0, 0.0, 1.0)
end

function amplification_factor(B::Real, params::BackgroundParams)
    l_0 = restoring_force(0.0, params)
    l_B = restoring_force(B, params)
    return l_0 / l_B
end

Base.@kwdef struct BackgroundStressor
    name::String
    value::Float64
    weight::Float64
end

function background_values(stressors::Vector{BackgroundStressor})
    return [s.value for s in stressors]
end

function background_weights(stressors::Vector{BackgroundStressor})
    return [s.weight for s in stressors]
end

function pairwise_sum(M::AbstractMatrix{Float64}, x::Vector{Float64})
    n = length(x)
    total = 0.0
    for j in 1:n
        for k in (j+1):n
            total += M[j, k] * x[j] * x[k]
        end
    end
    return total
end

function background_index(stressors::Vector{BackgroundStressor}; interaction=nothing)
    vals = background_values(stressors)
    w = background_weights(stressors)

    B_additive = sum(w .* vals)

    B_interaction = 0.0
    if interaction !== nothing
        B_interaction = pairwise_sum(interaction, vals)
    end

    B = B_additive + B_interaction
    return max(B, 0.0)
end
