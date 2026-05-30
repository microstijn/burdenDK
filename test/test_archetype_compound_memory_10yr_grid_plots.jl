using Test
using NCDatasets

@testset "Archetype Compound Memory 10yr Grid Plotting" begin
    # Gate the test
    run_plots = get(ENV, "TTR_RUN_PLOTTING_TESTS", "false") == "true" ||
                get(ENV, "TTR_RUN_EXAMPLE_TESTS", "false") == "true" ||
                get(ENV, "TTR_RUN_EXTENDED_TESTS", "false") == "true"

    if run_plots
        mktempdir() do temp_dir
            # Create a tiny temporary NetCDF file with mock variables
            nc_path = joinpath(temp_dir, "vulnerability_regime_outputs.nc")

            # Dimensions
            nx, ny = 5, 5

            NCDataset(nc_path, "c") do ds
                defDim(ds, "x", nx)
                defDim(ds, "y", ny)

                # Create a cluster_id array (categorical, integer IDs)
                cluster_data = rand(1:3, nx, ny)
                v_cluster = defVar(ds, "cluster_id", Int32, ("x", "y"))
                v_cluster[:, :] = cluster_data

                # Mock optional variables
                v_f = defVar(ds, "p95_F_grouped", Float64, ("x", "y"))
                v_f[:, :] = rand(Float64, nx, ny)

                v_q = defVar(ds, "p95_Q_grouped", Float64, ("x", "y"))
                v_q[:, :] = rand(Float64, nx, ny)

                v_ent = defVar(ds, "axis_entropy", Float64, ("x", "y"))
                v_ent[:, :] = rand(Float64, nx, ny)

                v_sens = defVar(ds, "max_delta_F_grouped_minus_TU", Float64, ("x", "y"))
                v_sens[:, :] = rand(Float64, nx, ny) .- 0.5

                v_seas = defVar(ds, "month_of_max_p95_F_grouped", Float64, ("x", "y"))
                v_seas[:, :] = rand(1.0:12.0, nx, ny)
            end

            # Run the script via command line to avoid world-age issues
            script_path = normpath(joinpath(@__DIR__, "..", "examples", "plot_archetype_compound_memory_10yr_grid_demo.jl"))

            # Use julia executable
            julia_cmd = Base.julia_cmd()
            env_vars = copy(ENV)
            env_vars["TTR_10YR_DEMO_OUTPUT_DIR"] = temp_dir
            cmd = setenv(`$julia_cmd --project=. $script_path`, env_vars)

            # Call the script
            run(cmd)

            # verify required PNG files are created and non-empty
            figures_dir = joinpath(temp_dir, "figures")
            @test isdir(figures_dir)

            expected_files = [
                "cluster_map.png",
                "p95_F_grouped_map.png",
                "p95_Q_grouped_map.png",
                "axis_entropy_map.png",
                "mixture_sensitivity_map.png",
                "seasonality_map.png",
                "vulnerability_regime_summary_panel.png"
            ]

            for file in expected_files
                path = joinpath(figures_dir, file)
                @test isfile(path)
                @test filesize(path) > 0
            end

            # CSV tests are implicitly skipped and handled gracefully since we didn't create them
        end
    else
        @info "Skipping TTR_RUN_PLOTTING_TESTS because it is not enabled."
        @test true
    end
end
