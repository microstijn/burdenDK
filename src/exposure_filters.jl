export ExposureFilter
export apply_exposure_filter, apply_exposure_filter_grid
export aquatic_exposure_filter, human_exposure_filter, default_exposure_filter
export surface_volume_retention, waterborne_stage_retention

Base.@kwdef struct ExposureFilter
    multipliers::Vector{Float64}
    name::String = "generic"
    description::String = ""
end

function apply_exposure_filter(values::AbstractVector{<:Real}, filter::ExposureFilter)
    if length(values) != length(filter.multipliers)
        throw(ArgumentError("Length of values ($(length(values))) must equal length of multipliers ($(length(filter.multipliers)))"))
    end
    return filter.multipliers .* values
end

function apply_exposure_filter_grid(layers::Vector{Matrix{Float64}}, filter::ExposureFilter)
    n_vars = length(layers)
    if n_vars == 0
        throw(ArgumentError("At least one layer must be provided"))
    end
    if n_vars != length(filter.multipliers)
        throw(ArgumentError("Number of layers ($n_vars) must equal length of multipliers ($(length(filter.multipliers)))"))
    end

    dims = size(layers[1])
    for i in 2:n_vars
        if size(layers[i]) != dims
            throw(ArgumentError("All layers must have the same size"))
        end
    end

    out_layers = [Matrix{Float64}(undef, dims) for _ in 1:n_vars]

    for c in 1:dims[2]
        for r in 1:dims[1]
            for i in 1:n_vars
                v = layers[i][r, c]
                if isnan(v)
                    out_layers[i][r, c] = NaN
                else
                    out_layers[i][r, c] = v * filter.multipliers[i]
                end
            end
        end
    end

    return out_layers
end

function aquatic_exposure_filter(nvars::Int)
    return ExposureFilter(
        multipliers = fill(1.0, nvars),
        name = "aquatic",
        description = "Direct aquatic exposure, multipliers are 1.0"
    )
end

function default_exposure_filter(nvars::Int)
    return aquatic_exposure_filter(nvars)
end

function human_exposure_filter()
    return ExposureFilter(
        multipliers = [0.10, 0.25, 0.50, 0.60, 0.35, 0.40, 0.30],
        name = "human",
        description = "This is a placeholder contact/use filter, not calibrated human health exposure."
    )
end

# --- Surface:volume waterborne toxicokinetics (aquatic, water-breathing targets only) -------
# A waterborne contaminant exchanges across the gill/skin SURFACE and dilutes into body VOLUME,
# so the one-compartment exchange rate scales with surface:volume ∝ 1/L. Because BOTH uptake and
# elimination scale the same way, the uptake/elimination ratio (steady-state burden) is
# size-independent; what changes with length is the RATE. Over a fixed month the exchange rate
# enters through the retention: rho(L) = exp(-k(L) Δt) with k(L) ∝ 1/L. Anchoring at a reference
# length L_ref (retention rho_ref there) gives rho(L) = rho_ref^(L_ref/L): a small juvenile
# (L < L_ref) equilibrates FAST (low retention, tracks the ambient water), a large adult LAGS
# (high retention, carries burden longer). No new knob -- just the surface:volume exponent on the
# existing per-compound retention. This mechanism is GILL/SKIN uptake and does NOT apply to
# terrestrial, air-breathing, or dietary-dominated routes.

"""
    surface_volume_retention(rho_ref, L, L_ref) -> Float64

Length-dependent monthly retention for waterborne (gill/skin) toxicokinetics: `rho_ref^(L_ref/L)`.
`rho_ref` is the retention at reference structural length `L_ref`. Smaller stages (`L < L_ref`)
get lower retention (faster equilibration); larger stages lag. Steady-state burden is unchanged
(uptake/elimination both scale ∝ 1/L), so only the rate moves with length.
"""
function surface_volume_retention(rho_ref::Real, L::Real, L_ref::Real)
    rr = Float64(rho_ref); Lf = Float64(L); Lr = Float64(L_ref)
    if !isfinite(rr) || !(0.0 <= rr < 1.0)
        throw(ArgumentError("rho_ref must be finite and satisfy 0 <= rho_ref < 1"))
    end
    if !isfinite(Lf) || Lf <= 0.0 || !isfinite(Lr) || Lr <= 0.0
        throw(ArgumentError("L and L_ref must be finite and > 0"))
    end
    rr == 0.0 && return 0.0
    return rr ^ (Lr / Lf)
end

"""
    waterborne_stage_retention(rho_ref, L, L_ref; waterborne) -> Float64

Gated surface:volume retention. Returns [`surface_volume_retention`](@ref) for a waterborne /
water-breathing aquatic target, or `rho_ref` unchanged otherwise -- the surface:volume mechanism
is gill/skin uptake and must not fire on terrestrial, air-breathing, or dietary-dominated targets.
"""
function waterborne_stage_retention(rho_ref::Real, L::Real, L_ref::Real; waterborne::Bool)
    return waterborne ? surface_volume_retention(rho_ref, L, L_ref) : Float64(rho_ref)
end
