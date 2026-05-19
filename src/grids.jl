function compute_background_index_grid(layers::Vector{Matrix{Float64}}, weights::Vector{Float64}; interaction=nothing)
    if isempty(layers)
        error("No layers provided")
    end

    nrows, ncols = size(layers[1])
    n_layers = length(layers)

    for i in 2:n_layers
        if size(layers[i]) != (nrows, ncols)
            error("All layers must have the same dimensions")
        end
    end

    if length(weights) != n_layers
        error("Number of weights must match number of layers")
    end

    Bgrid = zeros(Float64, nrows, ncols)

    for i in 1:nrows
        for j in 1:ncols
            vals = zeros(Float64, n_layers)
            for k in 1:n_layers
                vals[k] = layers[k][i, j]
            end

            B_additive = sum(weights .* vals)

            B_interaction = 0.0
            if interaction !== nothing
                for a in 1:n_layers
                    for b in (a+1):n_layers
                        B_interaction += interaction[a, b] * vals[a] * vals[b]
                    end
                end
            end

            Bgrid[i, j] = max(B_additive + B_interaction, 0.0)
        end
    end

    return Bgrid
end

function adaptive_margin_grid(Bgrid::Matrix{Float64}, params::BackgroundParams)
    return map(B -> adaptive_margin(B, params), Bgrid)
end

function restoring_force_grid(Bgrid::Matrix{Float64}, params::BackgroundParams)
    return map(B -> restoring_force(B, params), Bgrid)
end

function amplification_factor_grid(Bgrid::Matrix{Float64}, params::BackgroundParams)
    return map(B -> amplification_factor(B, params), Bgrid)
end

function synthetic_background_layers(nrows::Int, ncols::Int)
    layers = Matrix{Float64}[]

    # Layer 1: TDS-like gradient (west to east)
    L1 = zeros(Float64, nrows, ncols)
    for i in 1:nrows
        for j in 1:ncols
            L1[i, j] = (j - 1) / (ncols - 1)
        end
    end
    push!(layers, L1)

    # Layer 2: BOD-like gradient (south to north)
    L2 = zeros(Float64, nrows, ncols)
    for i in 1:nrows
        for j in 1:ncols
            L2[i, j] = (nrows - i) / (nrows - 1)
        end
    end
    push!(layers, L2)

    # Layer 3: FC-like hotspot in center
    L3 = zeros(Float64, nrows, ncols)
    ci, cj = nrows / 2.0, ncols / 2.0
    max_dist = sqrt(ci^2 + cj^2)
    for i in 1:nrows
        for j in 1:ncols
            dist = sqrt((i - ci)^2 + (j - cj)^2)
            L3[i, j] = max(0.0, 1.0 - dist / (max_dist * 0.5))
        end
    end
    push!(layers, L3)

    # Layer 4: temperature-like gradient
    L4 = zeros(Float64, nrows, ncols)
    for i in 1:nrows
        for j in 1:ncols
            L4[i, j] = 0.5 + 0.5 * sin(pi * j / ncols) * cos(pi * i / nrows)
        end
    end
    # Ensure bounds 0-1
    clamp!(L4, 0.0, 1.0)
    push!(layers, L4)

    return layers
end

function run_synthetic_raster_demo(params::BackgroundParams; output_dir="output/raster_demo")
    mkpath(output_dir)

    nrows, ncols = 50, 50
    layers = synthetic_background_layers(nrows, ncols)
    weights = [0.25, 0.25, 0.25, 0.25]

    # Mild interaction between L1 and L2
    interaction = zeros(4, 4)
    interaction[1, 2] = 0.5

    Bgrid = compute_background_index_grid(layers, weights; interaction=interaction)
    Agrid = adaptive_margin_grid(Bgrid, params)
    lambdagrid = restoring_force_grid(Bgrid, params)
    Fgrid = amplification_factor_grid(Bgrid, params)

    # Write ASCII grids
    write_ascii_grid(joinpath(output_dir, "background_index_B.asc"), Bgrid)
    write_ascii_grid(joinpath(output_dir, "adaptive_margin_A.asc"), Agrid)
    write_ascii_grid(joinpath(output_dir, "restoring_force_lambda.asc"), lambdagrid)
    write_ascii_grid(joinpath(output_dir, "amplification_factor.asc"), Fgrid)

    # Write PNGs
    plot_grid(Bgrid; title="Background Index (B)", filename=joinpath(output_dir, "background_index_B.png"), colormap=:viridis)
    plot_amplification_grid(Fgrid; filename=joinpath(output_dir, "amplification_factor.png"))

    return Bgrid, Agrid, lambdagrid, Fgrid
end
