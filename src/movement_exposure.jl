export occupancy_weighted_exposure

# Movement / mobile-target exposure.
#
# The resident model exposes a target to the field at a single cell. A mobile animal instead
# experiences the residence-time-weighted field along its trajectory. This primitive computes that
# weighted exposure; downstream it feeds the SAME compound-memory recurrence and (stage-resolved)
# adaptive-margin machinery as the resident case. A single occupied region reduces exactly to that
# region's concentration vector, so the resident case is the degenerate special case.

"""
    occupancy_weighted_exposure(region_concs, occupancy) -> Vector{Float64}

Residence-time-weighted exposure for a mobile target: given per-region contaminant concentration
vectors and an occupancy (residence-time) distribution over those regions, return the experienced
concentration vector `C = Σ_g π_g C_g` -- the path analogue of resident single-cell exposure.

`region_concs` is a `Vector` of equal-length per-region concentration vectors, or a `regions × J`
matrix. `occupancy` is a length-`regions` weight vector; it is validated (finite, non-negative,
not all-zero) and normalised to sum 1. A single occupied region returns that region's vector
exactly.
"""
function occupancy_weighted_exposure(region_concs::AbstractVector{<:AbstractVector{<:Real}},
                                     occupancy::AbstractVector{<:Real})
    nregions = length(region_concs)
    if nregions == 0
        throw(ArgumentError("region_concs must contain at least one region"))
    end
    if length(occupancy) != nregions
        throw(ArgumentError("occupancy length ($(length(occupancy))) must equal number of regions ($nregions)"))
    end

    J = length(region_concs[1])
    for (g, c) in enumerate(region_concs)
        if length(c) != J
            throw(ArgumentError("all region concentration vectors must have equal length; region $g has $(length(c)), expected $J"))
        end
        for v in c
            if !isfinite(v) || v < 0
                throw(ArgumentError("region concentrations must be finite and non-negative"))
            end
        end
    end

    w = _normalised_occupancy(occupancy)
    out = zeros(Float64, J)
    for g in 1:nregions
        wg = w[g]
        wg == 0.0 && continue
        cg = region_concs[g]
        for j in 1:J
            out[j] += wg * Float64(cg[j])
        end
    end
    return out
end

function occupancy_weighted_exposure(region_concs::AbstractMatrix{<:Real},
                                     occupancy::AbstractVector{<:Real})
    rows = [vec(region_concs[g, :]) for g in 1:size(region_concs, 1)]
    return occupancy_weighted_exposure(rows, occupancy)
end

function _normalised_occupancy(occupancy::AbstractVector{<:Real})
    if isempty(occupancy)
        throw(ArgumentError("occupancy must be non-empty"))
    end
    s = 0.0
    for w in occupancy
        if !isfinite(w) || w < 0
            throw(ArgumentError("occupancy weights must be finite and non-negative"))
        end
        s += w
    end
    if s <= 0.0
        throw(ArgumentError("occupancy weights must not be all zero"))
    end
    return [Float64(w) / s for w in occupancy]
end
