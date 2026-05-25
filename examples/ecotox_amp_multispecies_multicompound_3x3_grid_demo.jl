using Pkg
try
    Pkg.activate(joinpath(@__DIR__, ".."))
catch err
    @warn "Could not activate project" exception=(err, catch_backtrace())
end

using TwoTimescaleResilience
using CairoMakie
using Statistics

function finite_mean(A)
    vals = filter(isfinite, vec(A))
    isempty(vals) && return NaN
    return mean(vals)
end

function finite_max(A)
    vals = filter(isfinite, vec(A))
    isempty(vals) && return NaN
    return maximum(vals)
end

function finite_min(A)
    vals = filter(isfinite, vec(A))
    isempty(vals) && return NaN
    return minimum(vals)
end

function main()
    println("============================================================")
    println("TRANCHE 2: MULTI-SPECIES × MULTI-COMPOUND 3x3 GRID DEMO")
    println("============================================================")

    output_dir = joinpath(@__DIR__, "..", "output", "ecotox_amp_multispecies_multicompound_3x3_grid_demo")
    mkpath(output_dir)

    amp_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
    ecotox_path = joinpath(@__DIR__, "..", "data", "ECOTOX_Toxicity_Library.json")

    amp = load_amp_species_library(amp_path)
    ecotox = load_ecotox_library(ecotox_path)

    # ------------------------------------------------------------
    # Species selection
    # ------------------------------------------------------------
    preferred_species = ["Abatus_cordatus", "Podarcis_muralis", "Thalia_democratica"]
    selected_species_names = String[]
    species_params = Dict{String, Any}()

    for sp_name in preferred_species
        if haskey(amp, sp_name)
            push!(selected_species_names, sp_name)
            species_params[sp_name] = amp_species_deb_params(amp, sp_name)
        else
            @warn "Preferred species $sp_name missing from AmP library."
        end
    end

    if length(selected_species_names) < 3
        for (sp_name, data) in amp
            if !(sp_name in selected_species_names)
                try
                    params = amp_species_deb_params(amp, sp_name)
                    push!(selected_species_names, sp_name)
                    species_params[sp_name] = params
                    if length(selected_species_names) == 3
                        break
                    end
                catch
                end
            end
        end
    end

    println("\n--- Selected Species ---")
    for sp in selected_species_names
        println("- $sp")
    end

    # ------------------------------------------------------------
    # Compound selection
    # ------------------------------------------------------------
    preferred_compounds = ["7647-14-5", "7440-43-9", "7440-50-8"]
    selected_records = []

    available_cas = unique([get(rec, "cas_hyphenated", get(rec, "cas_norm", get(rec, "cas", ""))) for rec in ecotox])

    function find_valid_record(target_cas)
        recs = TwoTimescaleResilience.ecotox_filter_records(ecotox; cas=target_cas)
        for eff in ["MOR", "GRO", "REP"]
            for rec in recs
                if get(rec, "effect_code", "") == eff && validate_ecotox_record(rec)
                    return rec
                end
            end
        end
        for rec in recs
            if validate_ecotox_record(rec)
                return rec
            end
        end
        return nothing
    end

    selected_cas_list = String[]
    for cas in preferred_compounds
        rec = find_valid_record(cas)
        if rec !== nothing
            push!(selected_records, rec)
            push!(selected_cas_list, cas)
        else
            @warn "Preferred compound CAS $cas has no valid records. Will select fallback."
        end
    end

    if length(selected_records) < 3
        for cas in available_cas
            if !(cas in selected_cas_list) && cas != ""
                rec = find_valid_record(cas)
                if rec !== nothing
                    push!(selected_records, rec)
                    push!(selected_cas_list, cas)
                    if length(selected_records) == 3
                        break
                    end
                end
            end
        end
    end

    println("\n--- Selected Compound Records ---")
    for rec in selected_records
        cas_hyph = get(rec, "cas_hyphenated", "")
        eff = get(rec, "effect_code", "")
        routed = ecotox_effect_to_deb_axis(eff)
        println("- CAS: $cas_hyph | Effect: $eff | Routed Axis: $routed")
    end

    # ------------------------------------------------------------
    # 3x3 concentration grids & active stress grids
    # ------------------------------------------------------------
    compound_grids = Dict{String, Matrix{Float64}}()
    active_grids = Dict{String, Matrix{Float64}}()

    for rec in selected_records
        cas = get(rec, "cas_hyphenated", "")
        noec = Float64(get(rec, "NOEC_median", 0.0))
        ec50 = Float64(get(rec, "EC50_median", 0.0))

        Cgrid = [
            0.0              noec              (noec + ec50) / 2.0;
            noec             ec50              2.0 * ec50;
            (noec + ec50)/2.0  2.0 * ec50      5.0 * ec50
        ]

        compound_grids[cas] = Cgrid

        active_grid = zeros(Float64, 3, 3)
        for i in 1:3, j in 1:3
            active_grid[i, j] = ecotox_active_stress(Cgrid[i, j], noec, ec50)
        end

        active_grids[cas] = active_grid

        println("\n--- Grid summaries for $cas ---")
        println("Concentration grid:")
        display(Cgrid)
        println("Active stress grid:")
        display(active_grid)
    end

    # ------------------------------------------------------------
    # Compute aggregated burden grids
    # ------------------------------------------------------------
    assimilation_grid = zeros(Float64, 3, 3)
    maintenance_grid = zeros(Float64, 3, 3)
    growth_grid = zeros(Float64, 3, 3)
    reproduction_grid = zeros(Float64, 3, 3)

    for i in 1:3, j in 1:3
        concentrations = Dict{String, Float64}()
        for cas in selected_cas_list
            concentrations[cas] = compound_grids[cas][i, j]
        end
        burden = ecotox_records_to_deb_burden(concentrations, selected_records)
        assimilation_grid[i, j] = burden.assimilation
        maintenance_grid[i, j] = burden.maintenance
        growth_grid[i, j] = burden.growth
        reproduction_grid[i, j] = burden.reproduction
    end

    println("\n--- Aggregated Burden Grids ---")
    println("Assimilation:")
    display(assimilation_grid)
    println("Maintenance:")
    display(maintenance_grid)
    println("Growth:")
    display(growth_grid)
    println("Reproduction:")
    display(reproduction_grid)

    # ------------------------------------------------------------
    # Compute species-specific response grids
    # ------------------------------------------------------------
    species_A_grids = Dict{String, Matrix{Float64}}()
    species_lambda_grids = Dict{String, Matrix{Float64}}()
    species_F_grids = Dict{String, Matrix{Float64}}()

    for sp_name in selected_species_names
        Agrid = zeros(Float64, 3, 3)
        lambdagrid = zeros(Float64, 3, 3)
        Fgrid = zeros(Float64, 3, 3)

        params = species_params[sp_name]

        for i in 1:3, j in 1:3
            b = (
                assimilation = assimilation_grid[i, j],
                maintenance = maintenance_grid[i, j],
                growth = growth_grid[i, j],
                reproduction = reproduction_grid[i, j]
            )
            response = ecotox_burden_to_response(b, params)
            Agrid[i, j] = response.A
            lambdagrid[i, j] = response.lambda
            Fgrid[i, j] = response.amplification
        end

        species_A_grids[sp_name] = Agrid
        species_lambda_grids[sp_name] = lambdagrid
        species_F_grids[sp_name] = Fgrid

        println("\n--- Species Response Summary for $sp_name ---")
        println("  Mean Amplification: ", finite_mean(Fgrid))
        println("  Max Amplification:  ", finite_max(Fgrid))
        println("  Min A:              ", finite_min(Agrid))
        println("  Min lambda:         ", finite_min(lambdagrid))
    end

    # ------------------------------------------------------------
    # CairoMakie figures
    # ------------------------------------------------------------
    # These are synthetic 3x3 conceptual grids, so matrices are plotted as-is.

    # Figure 1: active stress grids per compound
    fig1 = Figure(size = (1200, 400))
    for (idx, cas) in enumerate(selected_cas_list)
        ax = Axis(fig1[1, idx], title="Active Stress\n$(cas)")
        hm = heatmap!(ax, active_grids[cas])
        Colorbar(fig1[1, idx][1, 2], hm)
    end
    save(joinpath(output_dir, "ecotox_active_stress_compounds.png"), fig1)

    # Figure 2: aggregated DEB burden grids
    fig2 = Figure(size = (800, 800))
    axes_burden = [
        (1, 1, "Assimilation", assimilation_grid),
        (1, 2, "Maintenance", maintenance_grid),
        (2, 1, "Growth", growth_grid),
        (2, 2, "Reproduction", reproduction_grid)
    ]
    for (r, c, name, grid) in axes_burden
        ax = Axis(fig2[r, c], title=name)
        hm = heatmap!(ax, grid)
        Colorbar(fig2[r, c][1, 2], hm)
    end
    save(joinpath(output_dir, "ecotox_aggregated_deb_burdens.png"), fig2)

    # Figure 3: species-specific amplification grids
    fig3 = Figure(size = (1200, 400))
    # Share colorrange
    min_F = minimum([finite_min(species_F_grids[sp]) for sp in selected_species_names])
    max_F = maximum([finite_max(species_F_grids[sp]) for sp in selected_species_names])
    # ensure it's at least [1, 1+eps] to avoid singular colormap if everything is 1.0
    if max_F - min_F < 1e-4
        max_F = min_F + 1e-4
    end

    for (idx, sp) in enumerate(selected_species_names)
        ax = Axis(fig3[1, idx], title=sp)
        hm = heatmap!(ax, species_F_grids[sp], colormap=:inferno, colorrange=(min_F, max_F))
        if idx == length(selected_species_names)
            Colorbar(fig3[1, idx][1, 2], hm)
        end
    end
    save(joinpath(output_dir, "ecotox_amp_species_amplification_grids.png"), fig3)

    # Figure 4: species-specific A grids
    fig4 = Figure(size = (1200, 400))
    for (idx, sp) in enumerate(selected_species_names)
        ax = Axis(fig4[1, idx], title="A grid: $sp")
        hm = heatmap!(ax, species_A_grids[sp])
        Colorbar(fig4[1, idx][1, 2], hm)
    end
    save(joinpath(output_dir, "ecotox_amp_species_A_grids.png"), fig4)


    println("\n--- Expected scientific behavior ---")
    println("- Each compound grid is scaled to that compound's own NOEC and EC50.")
    println("- Active stress is zero below or at NOEC.")
    println("- Active stress is near one at EC50.")
    println("- Aggregated DEB burden depends on ECOTOX effect-code routing.")
    println("- The same attack burden grid is applied to each AmP species.")
    println("- Species-specific AmP parameters produce different A, lambda, and amplification grids.")
    println("- This demonstrates that the framework supports multiple compounds and multiple species before full ISIMIP raster integration.")

    # ------------------------------------------------------------
    # Sanity checks
    # ------------------------------------------------------------
    println("\n--- Sanity Checks ---")
    checks_passed = true

    # 1. Active stress grids are finite
    for cas in selected_cas_list
        if !all(isfinite.(active_grids[cas]))
            @warn "Active stress grid for $cas contains non-finite values."
            checks_passed = false
        end
    end

    # 2. Maximum active stress >= 1 for each compound
    for cas in selected_cas_list
        if finite_max(active_grids[cas]) < 0.999
            @warn "Max active stress for $cas is less than 1."
            checks_passed = false
        end
    end

    # 3. Aggregated burden grids are nonnegative
    for grid in [assimilation_grid, maintenance_grid, growth_grid, reproduction_grid]
        if any(x -> x < 0.0, grid)
            @warn "Aggregated burden grid contains negative values."
            checks_passed = false
        end
    end

    # 4. Amplification grids are >= 1 within tolerance
    for sp in selected_species_names
        if finite_min(species_F_grids[sp]) < 0.9999
            @warn "Amplification grid for $sp contains values < 1."
            checks_passed = false
        end
    end

    # 5. Species-specific amplification grids are not all identical
    # Only if species params actually differ (likely they do)
    if length(selected_species_names) > 1
        sp1, sp2 = selected_species_names[1], selected_species_names[2]
        if species_F_grids[sp1] ≈ species_F_grids[sp2]
            @warn "Species amplification grids for $sp1 and $sp2 are exactly identical."
            # Not strictly failing the check, but worth noting
        end
    end

    if checks_passed
        println("3x3 multi-species multi-compound sanity checks passed.")
    else
        println("Some sanity checks failed. Check warnings.")
    end

    println("\nDone.")
end

main()
