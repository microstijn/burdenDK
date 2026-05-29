function analytical_initial_burden(rho::Real, K::Real, C_bg::Real, spinup_months::Integer; B0::Real = 0.0)
    rho_f = Float64(rho)
    K_f = Float64(K)
    C_bg_f = Float64(C_bg)
    B0_f = Float64(B0)

    if !(0.0 <= rho_f < 1.0) || !isfinite(rho_f)
        throw(ArgumentError("retention rho must be finite and satisfy 0.0 <= rho < 1.0"))
    end
    if K_f <= 0.0 || !isfinite(K_f)
        throw(ArgumentError("bioaccumulation factor K must be positive and finite"))
    end
    if C_bg_f < 0.0 || !isfinite(C_bg_f)
        throw(ArgumentError("background concentration C_bg must be >= 0 and finite"))
    end
    if spinup_months < 0
        throw(ArgumentError("spinup_months must be >= 0"))
    end
    if B0_f < 0.0 || !isfinite(B0_f)
        throw(ArgumentError("initial burden B0 must be >= 0 and finite"))
    end

    if spinup_months == 0
        return B0_f
    end

    if rho_f == 0.0
        return K_f * C_bg_f
    end

    rho_n = rho_f^spinup_months
    return rho_n * B0_f + K_f * C_bg_f * (1.0 - rho_n)
end

function background_for_target_burden(rho::Real, K::Real, B_target::Real, spinup_months::Integer; B0::Real = 0.0)
    rho_f = Float64(rho)
    K_f = Float64(K)
    B_target_f = Float64(B_target)
    B0_f = Float64(B0)

    if !(0.0 <= rho_f < 1.0) || !isfinite(rho_f)
        throw(ArgumentError("retention rho must be finite and satisfy 0.0 <= rho < 1.0"))
    end
    if K_f <= 0.0 || !isfinite(K_f)
        throw(ArgumentError("bioaccumulation factor K must be positive and finite"))
    end
    if B_target_f < 0.0 || !isfinite(B_target_f)
        throw(ArgumentError("target burden B_target must be >= 0 and finite"))
    end
    if spinup_months <= 0
        throw(ArgumentError("spinup_months must be > 0"))
    end
    if B0_f < 0.0 || !isfinite(B0_f)
        throw(ArgumentError("initial burden B0 must be >= 0 and finite"))
    end

    if rho_f == 0.0
        return B_target_f / K_f
    end

    rho_n = rho_f^spinup_months

    denom = K_f * (1.0 - rho_n)
    C_bg = (B_target_f - rho_n * B0_f) / denom

    if C_bg < 0.0
        throw(ArgumentError("unreachable target requiring negative C_bg"))
    end

    return C_bg
end

function analytical_periodic_initial_burden(rho::Real, K::Real, C_cycle::AbstractVector{<:Real})
    rho_f = Float64(rho)
    K_f = Float64(K)

    if !(0.0 <= rho_f < 1.0) || !isfinite(rho_f)
        throw(ArgumentError("retention rho must be finite and satisfy 0.0 <= rho < 1.0"))
    end
    if K_f <= 0.0 || !isfinite(K_f)
        throw(ArgumentError("bioaccumulation factor K must be positive and finite"))
    end
    if isempty(C_cycle)
        throw(ArgumentError("C_cycle must be non-empty"))
    end
    for c in C_cycle
        if !isfinite(c) || c < 0.0
            throw(ArgumentError("all values in C_cycle must be >= 0 and finite"))
        end
    end

    n = length(C_cycle)

    if rho_f == 0.0
        return K_f * Float64(C_cycle[end])
    end

    denom = 1.0 - rho_f^n
    weighted_sum = 0.0

    for i in 1:n
        weighted_sum += (rho_f^(n - i)) * Float64(C_cycle[i])
    end

    return ((1.0 - rho_f) * K_f * weighted_sum) / denom
end