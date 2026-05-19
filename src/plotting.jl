using CairoMakie

function plot_scenario_comparison(sims::Vector{TwoTimescaleSimulationResult}; filename="scenario_comparison.png")
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Response (y)", title = "Scenario Comparison")

    for sim in sims
        lines!(ax, sim.t, sim.y, label = sim.name, linewidth=2)
    end

    axislegend(ax)
    save(filename, fig)
    return fig
end

function plot_grid(grid::Matrix{Float64}; title="", filename="grid.png", colormap=:viridis)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title = title)

    hm = heatmap!(ax, grid, colormap=colormap)
    Colorbar(fig[1, 2], hm)

    save(filename, fig)
    return fig
end

function plot_background_layers(layers::Vector{Matrix{Float64}}, names::Vector{String}; filename="background_layers.png")
    n = length(layers)
    ncols = ceil(Int, sqrt(n))
    nrows = ceil(Int, n / ncols)

    fig = Figure(size = (800, 600))

    for i in 1:n
        r = div(i - 1, ncols) + 1
        c = rem(i - 1, ncols) + 1

        ax = Axis(fig[r, c], title = names[i])
        hm = heatmap!(ax, layers[i], colormap=:viridis)
        Colorbar(fig[r, c+1], hm)
    end

    save(filename, fig)
    return fig
end

function plot_amplification_grid(F::Matrix{Float64}; filename="amplification_factor.png")
    return plot_grid(F; title="Amplification Factor", filename=filename, colormap=:plasma)
end
