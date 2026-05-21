export MoAToDEBMapping
export moa_to_deb_axes, moa_to_deb_axes_grid
export moa_to_deb_axes_timeseries
export default_moa_to_deb_mapping

Base.@kwdef struct MoAToDEBMapping
    W::Matrix{Float64}
    interactions::Union{Nothing, Vector{Matrix{Float64}}} = nothing
    clamp_nonnegative::Bool = true
    clamp_unit::Bool = false
end

function moa_to_deb_axes(modes, mapping::MoAToDEBMapping)
    if modes isa NamedTuple
        if haskey(modes, :thermal) && length(modes) == 8
            # Canonical tuple
            mode_vals = [modes.thermal, modes.oxygen, modes.osmotic, modes.immune, modes.eutrophication, modes.toxic, modes.feeding, modes.physical]
        else
            mode_vals = [v for v in values(modes)]
        end
    else
        mode_vals = [v for v in modes]
    end

    n_axes, n_modes = size(mapping.W)
    if length(mode_vals) != n_modes
        throw(ArgumentError("Number of mode values ($(length(mode_vals))) must equal number of columns in W ($n_modes)"))
    end

    axes_vals = mapping.W * mode_vals

    if mapping.interactions !== nothing
        if length(mapping.interactions) != n_axes
            throw(ArgumentError("Length of interactions ($(length(mapping.interactions))) must equal number of rows in W ($n_axes)"))
        end
        for a in 1:n_axes
            int_mat = mapping.interactions[a]
            if size(int_mat, 1) != n_modes || size(int_mat, 2) != n_modes
                throw(ArgumentError("Interaction matrix $a size ($(size(int_mat))) must be ($n_modes, $n_modes)"))
            end
            for j in 1:(n_modes-1)
                for k in (j+1):n_modes
                    axes_vals[a] += int_mat[j, k] * mode_vals[j] * mode_vals[k]
                end
            end
        end
    end

    if mapping.clamp_nonnegative
        axes_vals .= max.(0.0, axes_vals)
    end
    if mapping.clamp_unit
        axes_vals .= min.(1.0, axes_vals)
    end

    if n_axes == 4
        return (
            assimilation = axes_vals[1],
            maintenance = axes_vals[2],
            growth = axes_vals[3],
            reproduction = axes_vals[4]
        )
    else
        return NamedTuple{Tuple(Symbol("axis_$i") for i in 1:n_axes)}(Tuple(axes_vals))
    end
end

function moa_to_deb_axes_grid(modes, mapping::MoAToDEBMapping)
    if modes isa NamedTuple
        if haskey(modes, :thermal) && length(modes) == 8
            layers = [modes.thermal, modes.oxygen, modes.osmotic, modes.immune, modes.eutrophication, modes.toxic, modes.feeding, modes.physical]
        else
            layers = [v for v in values(modes)]
        end
    else
        layers = [v for v in modes]
    end

    n_modes = length(layers)
    if n_modes == 0
        throw(ArgumentError("At least one mode layer must be provided"))
    end

    dims = size(layers[1])
    for i in 2:n_modes
        if size(layers[i]) != dims
            throw(ArgumentError("All mode layers must have the same size"))
        end
    end

    n_axes = size(mapping.W, 1)
    if n_modes != size(mapping.W, 2)
         throw(ArgumentError("Number of mode layers ($n_modes) must equal number of columns in W ($(size(mapping.W, 2)))"))
    end

    out_grids = [Matrix{Float64}(undef, dims) for _ in 1:n_axes]

    val_buf = zeros(Float64, n_modes)

    for c in 1:dims[2]
        for r in 1:dims[1]
            has_nan = false
            for i in 1:n_modes
                v = layers[i][r, c]
                if isnan(v)
                    has_nan = true
                    break
                end
                val_buf[i] = v
            end

            if has_nan
                for a in 1:n_axes
                    out_grids[a][r, c] = NaN
                end
            else
                ax_v = moa_to_deb_axes(val_buf, mapping)
                if ax_v isa NamedTuple
                    ax_vals = [v for v in values(ax_v)]
                else
                    ax_vals = ax_v
                end
                for a in 1:n_axes
                    out_grids[a][r, c] = ax_vals[a]
                end
            end
        end
    end

    if n_axes == 4
        return (
            assimilation = out_grids[1],
            maintenance = out_grids[2],
            growth = out_grids[3],
            reproduction = out_grids[4]
        )
    else
        return NamedTuple{Tuple(Symbol("axis_$i") for i in 1:n_axes)}(Tuple(out_grids))
    end
end

function moa_to_deb_axes_timeseries(mode_timeseries, mapping::MoAToDEBMapping)
    if mode_timeseries isa NamedTuple
        if haskey(mode_timeseries, :thermal) && length(mode_timeseries) == 8
            series = [mode_timeseries.thermal, mode_timeseries.oxygen, mode_timeseries.osmotic, mode_timeseries.immune, mode_timeseries.eutrophication, mode_timeseries.toxic, mode_timeseries.feeding, mode_timeseries.physical]
        else
            series = [v for v in values(mode_timeseries)]
        end
    else
        series = [v for v in mode_timeseries]
    end

    n_modes = length(series)
    if n_modes == 0
        throw(ArgumentError("At least one mode timeseries must be provided"))
    end

    n_time = length(series[1])
    for i in 2:n_modes
        if length(series[i]) != n_time
            throw(ArgumentError("All mode timeseries must have the same length"))
        end
    end

    n_axes = size(mapping.W, 1)
    if n_modes != size(mapping.W, 2)
         throw(ArgumentError("Number of mode timeseries ($n_modes) must equal number of columns in W ($(size(mapping.W, 2)))"))
    end

    out_series = [Vector{Float64}(undef, n_time) for _ in 1:n_axes]

    val_buf = zeros(Float64, n_modes)

    for t in 1:n_time
        for i in 1:n_modes
            val_buf[i] = series[i][t]
        end
        ax_v = moa_to_deb_axes(val_buf, mapping)
        if ax_v isa NamedTuple
            ax_vals = [v for v in values(ax_v)]
        else
            ax_vals = ax_v
        end
        for a in 1:n_axes
            out_series[a][t] = ax_vals[a]
        end
    end

    if n_axes == 4
        return (
            assimilation = out_series[1],
            maintenance = out_series[2],
            growth = out_series[3],
            reproduction = out_series[4]
        )
    else
        return NamedTuple{Tuple(Symbol("axis_$i") for i in 1:n_axes)}(Tuple(out_series))
    end
end

function default_moa_to_deb_mapping(; clamp_unit=false)
    W = [
        0.25  0.35  0.05  0.15  0.25  0.30  0.60  0.25;
        0.55  0.65  0.60  0.65  0.45  0.65  0.25  0.40;
        0.35  0.30  0.30  0.10  0.30  0.35  0.25  0.35;
        0.25  0.10  0.15  0.35  0.15  0.40  0.10  0.10
    ]
    return MoAToDEBMapping(
        W = W,
        interactions = nothing,
        clamp_nonnegative = true,
        clamp_unit = clamp_unit
    )
end
