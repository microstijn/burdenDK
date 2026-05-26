export DEBAxisParams, DEBAxisMapping, deb_axes, deb_adaptive_margin, ec50_anchored_fractional_impairment, axis_weights_for_species, compute_adaptive_margin_response, restoring_force_from_margin, amplification_from_margin, restoring_force_from_margin_and_axes
export deb_axes_grid, deb_adaptive_margin_grid, restoring_force_from_margin_grid, amplification_from_margin_grid, restoring_force_from_margin_and_axes_grid
export default_pathogen_organic_deb_mapping
export deb_amplification_pipeline
export pulse_deb_axes_timeseries, total_deb_margin_timeseries
export simulate_deb_axis_response

Base.@kwdef struct DEBAxisParams
    A0::Float64 = 1.0
    alpha_axes::NTuple{4, Float64} = (0.30, 0.35, 0.20, 0.15)
    lambda_min::Float64 = 0.04
    lambda_max::Float64 = 1.0
    KA::Float64 = 0.30

    recovery_axes::NTuple{4, Float64} = (0.10, 0.80, 0.10, 0.05)
    use_axis_recovery_penalty::Bool = false
    use_buffer_recovery_factor::Bool = false
    beta_Z::Float64 = 0.0
end

Base.@kwdef struct DEBAxisMapping
    # rows correspond to axes: assimilation, maintenance, growth, reproduction
    W::Matrix{Float64}
    # optional vector of four interaction matrices; each matrix is J x J
    interactions::Union{Nothing, Vector{Matrix{Float64}}} = nothing
    clamp_nonnegative::Bool = true
    clamp_unit::Bool = false
end

function pairwise_axis_interaction(Gamma::Matrix{Float64}, values::AbstractVector{<:Real})
    s = 0.0
    J = length(values)
    for j in 1:J
        for k in (j+1):J
            s += Gamma[j, k] * values[j] * values[k]
        end
    end
    return s
end

function deb_axes_vector(values::AbstractVector{<:Real}, mapping::DEBAxisMapping)::Vector{Float64}
    J = length(values)
    W = mapping.W
    s = W * values
    
    if !isnothing(mapping.interactions)
        for a in 1:4
            s[a] += pairwise_axis_interaction(mapping.interactions[a], values)
        end
    end
    
    if mapping.clamp_nonnegative
        s .= max.(s, 0.0)
    end
    
    if mapping.clamp_unit
        s .= min.(s, 1.0)
    end
    
    return s
end

function deb_axes(values::AbstractVector{<:Real}, mapping::DEBAxisMapping)
    s = deb_axes_vector(values, mapping)
    return (
        assimilation = s[1],
        maintenance = s[2],
        growth = s[3],
        reproduction = s[4]
    )
end

function _deb_axes_to_vector(axes)
    if axes isa NamedTuple
        return [
            Float64(axes.assimilation),
            Float64(axes.maintenance),
            Float64(axes.growth),
            Float64(axes.reproduction),
        ]
    elseif axes isa AbstractVector || axes isa Tuple
        if length(axes) != 4
            throw(ArgumentError("axes must have length 4"))
        end
        return Float64.(collect(axes))
    else
        throw(ArgumentError("axes must be a NamedTuple or length-4 vector/tuple"))
    end
end

function deb_adaptive_margin(axes, params::DEBAxisParams)
    s = _deb_axes_to_vector(axes)
    alpha = collect(params.alpha_axes)
    return params.A0 - sum(alpha .* s)
end

function restoring_force_from_margin(A::Real, params::DEBAxisParams)
    Ap = max(A, 0.0)
    return params.lambda_min + (params.lambda_max - params.lambda_min) * Ap / (params.KA + Ap)
end

function restoring_force_from_margin_and_axes(A::Real, axes, params::DEBAxisParams; Z=nothing)
    Ap = max(Float64(A), 0.0)

    base_lambda =
        params.lambda_min +
        (params.lambda_max - params.lambda_min) * Ap / (params.KA + Ap)

    penalty = 1.0

    if params.use_axis_recovery_penalty
        s = _deb_axes_to_vector(axes)
        beta = collect(params.recovery_axes)
        penalty *= exp(-sum(beta .* s))
    end

    if params.use_buffer_recovery_factor && Z !== nothing
        penalty *= exp(params.beta_Z * Float64(Z))
    end

    lambda = base_lambda * penalty

    return clamp(lambda, params.lambda_min, params.lambda_max)
end

function amplification_from_margin(A::Real, params::DEBAxisParams; A_control=params.A0)
    lambda_control = restoring_force_from_margin(A_control, params)
    lambda_A = restoring_force_from_margin(A, params)
    return lambda_control / lambda_A
end

function deb_axes_grid(layers::Vector{Matrix{Float64}}, mapping::DEBAxisMapping)
    nrows, ncols = size(layers[1])
    J = length(layers)
    
    sAgrid = zeros(Float64, nrows, ncols)
    sMgrid = zeros(Float64, nrows, ncols)
    sGgrid = zeros(Float64, nrows, ncols)
    sRgrid = zeros(Float64, nrows, ncols)
    
    for r in 1:nrows
        for c in 1:ncols
            has_nan = false
            for j in 1:J
                if isnan(layers[j][r, c])
                    has_nan = true
                    break
                end
            end
            
            if has_nan
                sAgrid[r, c] = NaN
                sMgrid[r, c] = NaN
                sGgrid[r, c] = NaN
                sRgrid[r, c] = NaN
            else
                b = zeros(Float64, J)
                for j in 1:J
                    b[j] = layers[j][r, c]
                end
                
                s = deb_axes_vector(b, mapping)
                sAgrid[r, c] = s[1]
                sMgrid[r, c] = s[2]
                sGgrid[r, c] = s[3]
                sRgrid[r, c] = s[4]
            end
        end
    end
    
    return (assimilation=sAgrid, maintenance=sMgrid, growth=sGgrid, reproduction=sRgrid)
end

function deb_adaptive_margin_grid(axes, params::DEBAxisParams)
    sA = axes.assimilation
    sM = axes.maintenance
    sG = axes.growth
    sR = axes.reproduction
    
    nrows, ncols = size(sA)
    Agrid = zeros(Float64, nrows, ncols)
    
    alphaA, alphaM, alphaG, alphaR = params.alpha_axes
    
    for r in 1:nrows
        for c in 1:ncols
            if isnan(sA[r, c]) || isnan(sM[r, c]) || isnan(sG[r, c]) || isnan(sR[r, c])
                Agrid[r, c] = NaN
            else
                Agrid[r, c] = params.A0 - (alphaA * sA[r, c] + alphaM * sM[r, c] + alphaG * sG[r, c] + alphaR * sR[r, c])
            end
        end
    end
    return Agrid
end

function restoring_force_from_margin_grid(Agrid::Matrix{Float64}, params::DEBAxisParams)
    nrows, ncols = size(Agrid)
    lambdagrid = zeros(Float64, nrows, ncols)
    
    for r in 1:nrows
        for c in 1:ncols
            if isnan(Agrid[r, c])
                lambdagrid[r, c] = NaN
            else
                lambdagrid[r, c] = restoring_force_from_margin(Agrid[r, c], params)
            end
        end
    end
    return lambdagrid
end

function restoring_force_from_margin_and_axes_grid(Agrid::Matrix{Float64}, axes, params::DEBAxisParams; Zgrid=nothing)
    sA = axes.assimilation
    sM = axes.maintenance
    sG = axes.growth
    sR = axes.reproduction

    nrows, ncols = size(Agrid)
    lambdagrid = zeros(Float64, nrows, ncols)

    for r in 1:nrows
        for c in 1:ncols
            if isnan(Agrid[r, c]) || isnan(sA[r, c]) || isnan(sM[r, c]) || isnan(sG[r, c]) || isnan(sR[r, c])
                lambdagrid[r, c] = NaN
            elseif Zgrid !== nothing && isnan(Zgrid[r, c])
                lambdagrid[r, c] = NaN
            else
                cell_axes = (
                    assimilation = sA[r, c],
                    maintenance = sM[r, c],
                    growth = sG[r, c],
                    reproduction = sR[r, c],
                )

                Z = Zgrid === nothing ? nothing : Zgrid[r, c]

                lambdagrid[r, c] =
                    restoring_force_from_margin_and_axes(
                        Agrid[r, c],
                        cell_axes,
                        params;
                        Z = Z
                    )
            end
        end
    end

    return lambdagrid
end

function amplification_from_margin_grid(Agrid::Matrix{Float64}, params::DEBAxisParams; A_control=params.A0)
    nrows, ncols = size(Agrid)
    Fgrid = zeros(Float64, nrows, ncols)
    lambda_control = restoring_force_from_margin(A_control, params)
    
    for r in 1:nrows
        for c in 1:ncols
            if isnan(Agrid[r, c])
                Fgrid[r, c] = NaN
            else
                lambda_A = restoring_force_from_margin(Agrid[r, c], params)
                Fgrid[r, c] = lambda_control / lambda_A
            end
        end
    end
    return Fgrid
end

"""
    default_pathogen_organic_deb_mapping(; interaction_strength=0.25, clamp_unit=false)

Returns a hypothetical, reduced-order DEB-inspired mapping for the `pathogen` and `organic` 
background stressors. 

- `pathogen` is interpreted as a faecal/pathogen proxy.
- `organic` is interpreted as an organic pollution/BOD proxy.

Note: This mapping is meant for hypothesis generation and sensitivity testing. No calibration 
claim is made.
"""
function default_pathogen_organic_deb_mapping(; interaction_strength=0.25, clamp_unit=false)
    W = [
        0.00 0.50;
        0.35 0.40;
        0.00 0.25;
        0.25 0.15
    ]
    
    Gamma_A = zeros(2, 2)
    Gamma_M = zeros(2, 2)
    Gamma_G = zeros(2, 2)
    Gamma_R = zeros(2, 2)
    
    Gamma_M[1, 2] = interaction_strength
    
    interactions = [Gamma_A, Gamma_M, Gamma_G, Gamma_R]
    
    return DEBAxisMapping(
        W=W,
        interactions=interactions,
        clamp_nonnegative=true,
        clamp_unit=clamp_unit
    )
end

function deb_amplification_pipeline(
    layers::Vector{Matrix{Float64}},
    mapping::DEBAxisMapping,
    params::DEBAxisParams;
    buffer_grid=nothing,
    buffer_params=nothing
)
    axes = deb_axes_grid(layers, mapping)
    
    nrows, ncols = size(layers[1])
    Agrid = zeros(Float64, nrows, ncols)

    if buffer_grid !== nothing && buffer_params !== nothing
        sA = axes.assimilation
        sM = axes.maintenance
        sG = axes.growth
        sR = axes.reproduction
        for r in 1:nrows
            for c in 1:ncols
                if isnan(buffer_grid[r, c]) || isnan(sA[r, c]) || isnan(sM[r, c]) || isnan(sG[r, c]) || isnan(sR[r, c])
                    Agrid[r, c] = NaN
                else
                    cell_axes = (
                        assimilation = sA[r, c],
                        maintenance = sM[r, c],
                        growth = sG[r, c],
                        reproduction = sR[r, c],
                    )
                    Agrid[r, c] = adaptive_margin_with_buffer(cell_axes, buffer_grid[r, c], params, buffer_params)
                end
            end
        end
    else
        Agrid = deb_adaptive_margin_grid(axes, params)
    end

    if params.use_axis_recovery_penalty || params.use_buffer_recovery_factor
        lambdagrid = restoring_force_from_margin_and_axes_grid(Agrid, axes, params; Zgrid=buffer_grid)

        lambda_control = restoring_force_from_margin_and_axes(
            params.A0,
            (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0),
            params;
            Z = 0.0
        )
    else
        lambdagrid = restoring_force_from_margin_grid(Agrid, params)
        lambda_control = restoring_force_from_margin(params.A0, params)
    end

    Fgrid = zeros(Float64, nrows, ncols)
    for r in 1:nrows
        for c in 1:ncols
            if isnan(lambdagrid[r, c]) || (buffer_grid !== nothing && isnan(buffer_grid[r, c]))
                Fgrid[r, c] = NaN
            else
                Fgrid[r, c] = lambda_control / lambdagrid[r, c]
            end
        end
    end

    if buffer_grid !== nothing
        return (
            axes = axes,
            Z = buffer_grid,
            A = Agrid,
            lambda = lambdagrid,
            amplification = Fgrid
        )
    else
        return (
            axes = axes,
            A = Agrid,
            lambda = lambdagrid,
            amplification = Fgrid
        )
    end
end

function pulse_deb_axes_timeseries(D::Matrix{Float64}, mapping::DEBAxisMapping)
    ntime, nstressors = size(D)
    
    sA_pulse = zeros(Float64, ntime)
    sM_pulse = zeros(Float64, ntime)
    sG_pulse = zeros(Float64, ntime)
    sR_pulse = zeros(Float64, ntime)
    
    for t in 1:ntime
        b = D[t, :]
        s = deb_axes_vector(b, mapping)
        
        sA_pulse[t] = s[1]
        sM_pulse[t] = s[2]
        sG_pulse[t] = s[3]
        sR_pulse[t] = s[4]
    end
    
    return (
        assimilation = sA_pulse,
        maintenance = sM_pulse,
        growth = sG_pulse,
        reproduction = sR_pulse
    )
end

function total_deb_margin_timeseries(A_background::Real, pulse_axes, params::DEBAxisParams)
    sA_pulse = pulse_axes.assimilation
    sM_pulse = pulse_axes.maintenance
    sG_pulse = pulse_axes.growth
    sR_pulse = pulse_axes.reproduction
    
    ntime = length(sA_pulse)
    A_total = zeros(Float64, ntime)
    
    alphaA, alphaM, alphaG, alphaR = params.alpha_axes
    
    for t in 1:ntime
        cost = alphaA * sA_pulse[t] + alphaM * sM_pulse[t] + alphaG * sG_pulse[t] + alphaR * sR_pulse[t]
        A_total[t] = A_background - cost
    end
    
    return A_total
end

function simulate_deb_axis_response(t::Vector{Float64}, pulse_axes, A_background::Float64, params::DEBAxisParams, q::Float64; y0::Float64 = 0.0, dt::Float64)
    ntime = length(t)
    y = zeros(Float64, ntime)
    
    sA_pulse = pulse_axes.assimilation
    sM_pulse = pulse_axes.maintenance
    sG_pulse = pulse_axes.growth
    sR_pulse = pulse_axes.reproduction
    
    alphaA, alphaM, alphaG, alphaR = params.alpha_axes
    
    y_current = y0
    if ntime > 0
        y[1] = y_current
    end
    
    for i in 1:(ntime - 1)
        cost_i = alphaA * sA_pulse[i] + alphaM * sM_pulse[i] + alphaG * sG_pulse[i] + alphaR * sR_pulse[i]
        A_total_i = A_background - cost_i
        lambda_i = restoring_force_from_margin(A_total_i, params)
        
        if lambda_i > 0
            y_next = y_current * exp(-lambda_i * dt) + (q * cost_i / lambda_i) * (1.0 - exp(-lambda_i * dt))
        else
            y_next = y_current + q * cost_i * dt
        end
        
        y[i+1] = y_next
        y_current = y_next
    end
    
    return y
end

function ec50_anchored_fractional_impairment(axis_pressures)
    s = _deb_axes_to_vector(axis_pressures)
    
    for x in s
        if !isfinite(x) || x < 0.0
            throw(ArgumentError("Axis pressure must be finite and >= 0"))
        end
    end
    
    return (
        assimilation = s[1] / (1.0 + s[1]),
        maintenance = s[2] / (1.0 + s[2]),
        growth = s[3] / (1.0 + s[3]),
        reproduction = s[4] / (1.0 + s[4])
    )
end

function axis_weights_for_species(params::DEBAxisParams; method::String = "auto", axis_weight_scope::String = "all_axes")
    if method == "auto" || method == "normalized_alpha_axes"
        alphas = params.alpha_axes
        sum_alphas = sum(alphas)
        
        # Validation
        for a in alphas
            if !isfinite(a) || a < 0.0
                if method == "normalized_alpha_axes"
                    throw(ArgumentError("Alpha axes values must be finite and non-negative"))
                else
                    return (
                        w_assimilation = 0.25,
                        w_maintenance = 0.25,
                        w_growth = 0.25,
                        w_reproduction = 0.25,
                        axis_weight_method = "equal_weight_diagnostic_fallback",
                        axis_weight_scope = axis_weight_scope
                    )
                end
            end
        end
        
        if sum_alphas <= 0.0
            if method == "normalized_alpha_axes"
                throw(ArgumentError("Sum of alpha axes must be > 0"))
            else
                return (
                    w_assimilation = 0.25,
                    w_maintenance = 0.25,
                    w_growth = 0.25,
                    w_reproduction = 0.25,
                    axis_weight_method = "equal_weight_diagnostic_fallback",
                    axis_weight_scope = axis_weight_scope
                )
            end
        end
        
        return (
            w_assimilation = alphas[1] / sum_alphas,
            w_maintenance = alphas[2] / sum_alphas,
            w_growth = alphas[3] / sum_alphas,
            w_reproduction = alphas[4] / sum_alphas,
            axis_weight_method = "normalized_alpha_axes",
            axis_weight_scope = axis_weight_scope
        )
    elseif method == "equal_weight_diagnostic_fallback"
        return (
            w_assimilation = 0.25,
            w_maintenance = 0.25,
            w_growth = 0.25,
            w_reproduction = 0.25,
            axis_weight_method = "equal_weight_diagnostic_fallback",
            axis_weight_scope = axis_weight_scope
        )
    else
        throw(ArgumentError("Unsupported axis weight method: $method"))
    end
end

function compute_adaptive_margin_response(axis_pressures, params::DEBAxisParams; 
    response_mode::String = "raw_margin_subtraction", 
    A_floor_fraction::Float64 = 1e-6)
    
    if !isfinite(A_floor_fraction) || A_floor_fraction <= 0.0 || A_floor_fraction > 1.0
        throw(ArgumentError("A_floor_fraction must be finite and satisfy 0 < A_floor_fraction <= 1"))
    end
    
    s = _deb_axes_to_vector(axis_pressures)
    
    if response_mode == "raw_margin_subtraction"
        A_t = deb_adaptive_margin(axis_pressures, params)
        lambda_t = restoring_force_from_margin(A_t, params)
        lambda0 = restoring_force_from_margin(params.A0, params)
        F_t = lambda0 / lambda_t
        
        # Calculate derived metrics for audit consistency
        E_res = ec50_anchored_fractional_impairment(axis_pressures)
        w_res = axis_weights_for_species(params)
        
        Q_t = w_res.w_assimilation * E_res.assimilation +
              w_res.w_maintenance * E_res.maintenance +
              w_res.w_growth * E_res.growth +
              w_res.w_reproduction * E_res.reproduction
              
        return (
            response_mode = response_mode,
            axis_weight_method = w_res.axis_weight_method,
            axis_weight_scope = w_res.axis_weight_scope,
            
            X_assimilation = s[1],
            X_maintenance = s[2],
            X_growth = s[3],
            X_reproduction = s[4],
            
            E_assimilation = E_res.assimilation,
            E_maintenance = E_res.maintenance,
            E_growth = E_res.growth,
            E_reproduction = E_res.reproduction,
            
            w_assimilation = w_res.w_assimilation,
            w_maintenance = w_res.w_maintenance,
            w_growth = w_res.w_growth,
            w_reproduction = w_res.w_reproduction,
            
            Q_t = Q_t,
            A0 = params.A0,
            A_t = A_t,
            lambda0 = lambda0,
            lambda_t = lambda_t,
            F_t = F_t,
            
            # Aliases for backwards compatibility
            A = A_t,
            lambda = lambda_t,
            amplification = F_t
        )
    elseif response_mode == "ec50_anchored_fractional_impairment"
        E_res = ec50_anchored_fractional_impairment(axis_pressures)
        w_res = axis_weights_for_species(params)
        
        Q_t = w_res.w_assimilation * E_res.assimilation +
              w_res.w_maintenance * E_res.maintenance +
              w_res.w_growth * E_res.growth +
              w_res.w_reproduction * E_res.reproduction
              
        A_t = params.A0 * max(A_floor_fraction, 1.0 - Q_t)
        lambda_t = restoring_force_from_margin(A_t, params)
        lambda0 = restoring_force_from_margin(params.A0, params)
        F_t = lambda0 / lambda_t
        
        return (
            response_mode = response_mode,
            axis_weight_method = w_res.axis_weight_method,
            axis_weight_scope = w_res.axis_weight_scope,
            
            X_assimilation = s[1],
            X_maintenance = s[2],
            X_growth = s[3],
            X_reproduction = s[4],
            
            E_assimilation = E_res.assimilation,
            E_maintenance = E_res.maintenance,
            E_growth = E_res.growth,
            E_reproduction = E_res.reproduction,
            
            w_assimilation = w_res.w_assimilation,
            w_maintenance = w_res.w_maintenance,
            w_growth = w_res.w_growth,
            w_reproduction = w_res.w_reproduction,
            
            Q_t = Q_t,
            A0 = params.A0,
            A_t = A_t,
            lambda0 = lambda0,
            lambda_t = lambda_t,
            F_t = F_t,
            
            # Aliases for backwards compatibility
            A = A_t,
            lambda = lambda_t,
            amplification = F_t
        )
    else
        throw(ArgumentError("Unknown response mode: $response_mode"))
    end
end
