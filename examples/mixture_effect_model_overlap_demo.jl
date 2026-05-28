# Mixture Effect Model Overlap Demo
#
# This example uses synthetic axis-burden records to demonstrate mixture-effect arithmetic only.
# It does not create or modify ECOTOX records.
# It does not implement synergism or antagonism.
# It does not implement fitted interaction coefficients.
# It does not implement DEBtox scaled damage.
# It does not implement raster integration.
#
# The example demonstrates that:
# - with one contributor, all models agree;
# - with distinct same-axis groups, IA and grouped can exceed TU;
# - with one shared group, grouped equals TU;
# - with mixed grouping, grouped can lie between TU and IA.

using TwoTimescaleResilience
using CairoMakie
using CSV
using DataFrames

println("Starting Mixture Effect Model Overlap Demo...")

# Ensure output directory exists
out_dir = normpath(joinpath(@__DIR__, "..", "output", "mixture_effect_model_overlap_demo"))
mkpath(out_dir)

# Define scenarios
scenarios_config = [
    (
        id = "single_contributor_control",
        desc = "One compound on maintenance",
        records = [
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "A", chemical_name = "compound_A", cas_norm = "synthetic_A", cas_hyphenated = "synthetic_A")
        ]
    ),
    (
        id = "two_distinct_effect_groups",
        desc = "Two compounds on maintenance, distinct groups",
        records = [
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "A", chemical_name = "compound_A", cas_norm = "synthetic_A", cas_hyphenated = "synthetic_A"),
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "B", chemical_name = "compound_B", cas_norm = "synthetic_B", cas_hyphenated = "synthetic_B")
        ]
    ),
    (
        id = "two_same_effect_group",
        desc = "Two compounds on maintenance, same group",
        records = [
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "A", chemical_name = "compound_A", cas_norm = "synthetic_A", cas_hyphenated = "synthetic_A"),
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "A", chemical_name = "compound_B", cas_norm = "synthetic_B", cas_hyphenated = "synthetic_B")
        ]
    ),
    (
        id = "mixed_grouping",
        desc = "Three compounds on maintenance, mixed grouping",
        records = [
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "A", chemical_name = "compound_A", cas_norm = "synthetic_A", cas_hyphenated = "synthetic_A"),
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "A", chemical_name = "compound_B", cas_norm = "synthetic_B", cas_hyphenated = "synthetic_B"),
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "B", chemical_name = "compound_C", cas_norm = "synthetic_C", cas_hyphenated = "synthetic_C")
        ]
    ),
    (
        id = "multiple_axes",
        desc = "Compounds on maintenance and growth",
        records = [
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "A", chemical_name = "compound_A", cas_norm = "synthetic_A", cas_hyphenated = "synthetic_A"),
            (burden_assimilation = 0.0, burden_maintenance = 1.0, burden_growth = 0.0, burden_reproduction = 0.0, effect_code = "B", chemical_name = "compound_B", cas_norm = "synthetic_B", cas_hyphenated = "synthetic_B"),
            (burden_assimilation = 0.0, burden_maintenance = 0.0, burden_growth = 4.0, burden_reproduction = 0.0, effect_code = "C", chemical_name = "compound_C", cas_norm = "synthetic_C", cas_hyphenated = "synthetic_C")
        ]
    )
]

models = [
    "axis_toxic_unit_sum",
    "independent_action_axis_effects",
    "grouped_ca_then_ia_axis_effects"
]

summary_rows = []
long_rows = []

for s in scenarios_config
    for m in models
        res = aggregate_axis_mixture_effects(s.records; mixture_effect_model = m)

        push!(summary_rows, (
            scenario = s.id,
            scenario_description = s.desc,
            mixture_effect_model = m,
            X_assimilation = res.X_assimilation,
            X_maintenance = res.X_maintenance,
            X_growth = res.X_growth,
            X_reproduction = res.X_reproduction,
            E_assimilation = res.E_assimilation,
            E_maintenance = res.E_maintenance,
            E_growth = res.E_growth,
            E_reproduction = res.E_reproduction
        ))

        for ax in ["assimilation", "maintenance", "growth", "reproduction"]
            push!(long_rows, (
                scenario = s.id,
                scenario_description = s.desc,
                mixture_effect_model = m,
                axis = ax,
                X_axis = getproperty(res, Symbol("X_$ax")),
                E_axis = getproperty(res, Symbol("E_$ax"))
            ))
        end
    end
end

df_summary = DataFrame(summary_rows)
df_long = DataFrame(long_rows)

summary_csv_path = joinpath(out_dir, "mixture_effect_model_overlap_summary.csv")
long_csv_path = joinpath(out_dir, "mixture_effect_model_overlap_axis_long.csv")

CSV.write(summary_csv_path, df_summary)
CSV.write(long_csv_path, df_long)

println("Saved summary CSV to: $summary_csv_path")
println("Saved long CSV to: $long_csv_path")

# Plotting
fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1],
    title = "Mixture Effect Model Overlap (Maintenance Axis)",
    xlabel = "Scenario",
    ylabel = "E_maintenance",
    xticks = (1:length(scenarios_config), [s.id for s in scenarios_config]),
    xticklabelrotation = pi/4
)

# Extract points for plotting
for (i, m) in enumerate(models)
    y_vals = Float64[]
    for s in scenarios_config
        val = filter(row -> row.scenario == s.id && row.mixture_effect_model == m, df_summary)[1, :E_maintenance]
        push!(y_vals, val)
    end

    # Give clear styles
    linestyle = if m == "axis_toxic_unit_sum"
        :solid
    elseif m == "independent_action_axis_effects"
        :dash
    else
        :dot
    end

    marker = if m == "axis_toxic_unit_sum"
        :circle
    elseif m == "independent_action_axis_effects"
        :rect
    else
        :diamond
    end

    lines!(ax, 1:length(scenarios_config), y_vals, label = m, linestyle = linestyle, linewidth = 2)
    scatter!(ax, 1:length(scenarios_config), y_vals, label = m, marker = marker, markersize = 12)
end

axislegend(ax, position = :lt)

plot_path = joinpath(out_dir, "mixture_effect_model_overlap_comparison.png")
save(plot_path, fig)
println("Saved plot to: $plot_path")

println("Done.")
