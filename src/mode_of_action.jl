export ModeOfActionMapping
export mode_of_action, mode_of_action_vector, mode_of_action_grid
export mode_names, isimip_variable_names

mode_names() = (:thermal, :oxygen, :osmotic, :immune, :eutrophication, :toxic, :feeding, :physical)
isimip_variable_names() = (:WT, :BOD, :TDS, :FC, :Nutrient, :Chemical, :Plastic)

Base.@kwdef struct ModeOfActionMapping
    U::Matrix{Float64}
    interactions::Union{Nothing, Vector{Matrix{Float64}}} = nothing
    clamp_nonnegative::Bool = true
    clamp_unit::Bool = false
end

function mode_of_action_vector(values::AbstractVector{<:Real}, mapping::ModeOfActionMapping)::Vector{Float64}
    n_modes, n_vars = size(mapping.U)
    if length(values) != n_vars
        throw(ArgumentError("Length of values ($(length(values))) must equal number of columns in U ($n_vars)"))
    end

    modes = mapping.U * values

    if mapping.interactions !== nothing
        if length(mapping.interactions) != n_modes
            throw(ArgumentError("Length of interactions ($(length(mapping.interactions))) must equal number of rows in U ($n_modes)"))
        end
        for r in 1:n_modes
            int_mat = mapping.interactions[r]
            if size(int_mat, 1) != n_vars || size(int_mat, 2) != n_vars
                throw(ArgumentError("Interaction matrix $r size ($(size(int_mat))) must be ($n_vars, $n_vars)"))
            end
            for j in 1:(n_vars-1)
                for k in (j+1):n_vars
                    modes[r] += int_mat[j, k] * values[j] * values[k]
                end
            end
        end
    end

    if mapping.clamp_nonnegative
        modes .= max.(0.0, modes)
    end
    if mapping.clamp_unit
        modes .= min.(1.0, modes)
    end

    return modes
end

function mode_of_action(values::AbstractVector{<:Real}, mapping::ModeOfActionMapping)
    m = mode_of_action_vector(values, mapping)
    if length(m) == 8
        return (
            thermal = m[1],
            oxygen = m[2],
            osmotic = m[3],
            immune = m[4],
            eutrophication = m[5],
            toxic = m[6],
            feeding = m[7],
            physical = m[8]
        )
    else
        # For non-canonical mappings, just return named tuple with indexed names
        return NamedTuple{Tuple(Symbol("mode_$i") for i in 1:length(m))}(Tuple(m))
    end
end

function mode_of_action_grid(layers::Vector{Matrix{Float64}}, mapping::ModeOfActionMapping)
    n_vars = length(layers)
    if n_vars == 0
        throw(ArgumentError("At least one layer must be provided"))
    end

    dims = size(layers[1])
    for i in 2:n_vars
        if size(layers[i]) != dims
            throw(ArgumentError("All layers must have the same size"))
        end
    end

    n_modes = size(mapping.U, 1)
    if n_vars != size(mapping.U, 2)
         throw(ArgumentError("Number of layers ($n_vars) must equal number of columns in U ($(size(mapping.U, 2)))"))
    end

    out_grids = [Matrix{Float64}(undef, dims) for _ in 1:n_modes]

    val_buf = zeros(Float64, n_vars)

    for c in 1:dims[2]
        for r in 1:dims[1]
            has_nan = false
            for i in 1:n_vars
                v = layers[i][r, c]
                if isnan(v)
                    has_nan = true
                    break
                end
                val_buf[i] = v
            end

            if has_nan
                for m in 1:n_modes
                    out_grids[m][r, c] = NaN
                end
            else
                mv = mode_of_action_vector(val_buf, mapping)
                for m in 1:n_modes
                    out_grids[m][r, c] = mv[m]
                end
            end
        end
    end

    if n_modes == 8
        return (
            thermal = out_grids[1],
            oxygen = out_grids[2],
            osmotic = out_grids[3],
            immune = out_grids[4],
            eutrophication = out_grids[5],
            toxic = out_grids[6],
            feeding = out_grids[7],
            physical = out_grids[8]
        )
    else
        return NamedTuple{Tuple(Symbol("mode_$i") for i in 1:n_modes)}(Tuple(out_grids))
    end
end
