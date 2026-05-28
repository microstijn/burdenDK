# src/mixture_aggregation.jl

export aggregate_deb_axis_burdens, mixture_contribution_diagnostics, aggregate_axis_mixture_effects

"""
    aggregate_axis_mixture_effects(compound_axis_burdens; mixture_effect_model = "axis_toxic_unit_sum")

Aggregates compound-specific DEB-axis burdens into axis-level mixture effects using the specified effect model.
Outputs bounded fractional axis impairments E_a.
"""
function aggregate_axis_mixture_effects end

"""
    aggregate_deb_axis_burdens(compound_axis_burdens; mixture_method = "additive_axis_burden")

Aggregates compound-specific DEB-axis burdens according to the specified mixture method.

Supported methods:
- "additive_axis_burden": Legacy/default baseline preserving the previous exact axis-wise summation behavior.
- "axis_toxic_unit_sum": Explicit analytical same-axis summation of model-defined active stress units after compound memory and ECOTOX threshold scaling.

Note:
This is exactly same-axis summation.
It does NOT implement formal Concentration Addition over raw external concentrations.
It does NOT implement Independent Action.
It does NOT implement bounded nonlinear aggregation.
It does NOT implement synergism or antagonism.
It does NOT implement low-effect approximation flags.
It does NOT implement physiological condition memory Z_t.
It does NOT implement raster integration.

The input `compound_axis_burdens` should be an iterable of objects (e.g., NamedTuples) containing:
- `burden_assimilation`
- `burden_maintenance`
- `burden_growth`
- `burden_reproduction`

Returns a NamedTuple:
`(total_burden_assimilation, total_burden_maintenance, total_burden_growth, total_burden_reproduction)`
"""
function aggregate_deb_axis_burdens(compound_axis_burdens; mixture_method::String = "additive_axis_burden")
    if mixture_method == "additive_axis_burden" || mixture_method == "axis_toxic_unit_sum"
        total_A = 0.0
        total_M = 0.0
        total_G = 0.0
        total_R = 0.0

        for r in compound_axis_burdens
            total_A += Float64(r.burden_assimilation)
            total_M += Float64(r.burden_maintenance)
            total_G += Float64(r.burden_growth)
            total_R += Float64(r.burden_reproduction)
        end

        return (
            total_burden_assimilation = total_A,
            total_burden_maintenance = total_M,
            total_burden_growth = total_G,
            total_burden_reproduction = total_R
        )
    else
        throw(ArgumentError("Unknown mixture_method: $mixture_method"))
    end
end

"""
    mixture_contribution_diagnostics(compound_axis_burdens)

Computes exact algebraic diagnostic summaries describing compound contributions by DEB axis.

For each axis, computes:
- `n_compounds_contributing_AXIS`: Number of compounds with burden > 0.
- `dominant_compound_AXIS`: Name of the compound with the largest burden (or "none" if total is 0).
- `max_single_compound_fraction_AXIS`: Max single compound fraction on that axis (or 0.0 if total is 0).

The input `compound_axis_burdens` should be an iterable of objects (e.g., NamedTuples) containing:
- `chemical_name`
- `burden_assimilation`
- `burden_maintenance`
- `burden_growth`
- `burden_reproduction`

Ties are broken by preferring the first encountered compound with the maximal burden.
"""
function mixture_contribution_diagnostics(compound_axis_burdens)
    # Totals
    agg = aggregate_deb_axis_burdens(compound_axis_burdens; mixture_method="additive_axis_burden")

    counts = Dict(:assimilation => 0, :maintenance => 0, :growth => 0, :reproduction => 0)
    max_val = Dict(:assimilation => 0.0, :maintenance => 0.0, :growth => 0.0, :reproduction => 0.0)
    dominant = Dict(:assimilation => "none", :maintenance => "none", :growth => "none", :reproduction => "none")

    for r in compound_axis_burdens
        name = r.chemical_name

        A_val = Float64(r.burden_assimilation)
        M_val = Float64(r.burden_maintenance)
        G_val = Float64(r.burden_growth)
        R_val = Float64(r.burden_reproduction)

        # Assimilation
        if A_val > 0
            counts[:assimilation] += 1
            if A_val > max_val[:assimilation]
                max_val[:assimilation] = A_val
                dominant[:assimilation] = name
            end
        end

        # Maintenance
        if M_val > 0
            counts[:maintenance] += 1
            if M_val > max_val[:maintenance]
                max_val[:maintenance] = M_val
                dominant[:maintenance] = name
            end
        end

        # Growth
        if G_val > 0
            counts[:growth] += 1
            if G_val > max_val[:growth]
                max_val[:growth] = G_val
                dominant[:growth] = name
            end
        end

        # Reproduction
        if R_val > 0
            counts[:reproduction] += 1
            if R_val > max_val[:reproduction]
                max_val[:reproduction] = R_val
                dominant[:reproduction] = name
            end
        end
    end

    frac_A = agg.total_burden_assimilation > 0 ? (max_val[:assimilation] / agg.total_burden_assimilation) : 0.0
    frac_M = agg.total_burden_maintenance > 0 ? (max_val[:maintenance] / agg.total_burden_maintenance) : 0.0
    frac_G = agg.total_burden_growth > 0 ? (max_val[:growth] / agg.total_burden_growth) : 0.0
    frac_R = agg.total_burden_reproduction > 0 ? (max_val[:reproduction] / agg.total_burden_reproduction) : 0.0

    return (
        n_compounds_contributing_assimilation = counts[:assimilation],
        n_compounds_contributing_maintenance = counts[:maintenance],
        n_compounds_contributing_growth = counts[:growth],
        n_compounds_contributing_reproduction = counts[:reproduction],

        dominant_compound_assimilation = dominant[:assimilation],
        dominant_compound_maintenance = dominant[:maintenance],
        dominant_compound_growth = dominant[:growth],
        dominant_compound_reproduction = dominant[:reproduction],

        max_single_compound_fraction_assimilation = frac_A,
        max_single_compound_fraction_maintenance = frac_M,
        max_single_compound_fraction_growth = frac_G,
        max_single_compound_fraction_reproduction = frac_R
    )
end
function aggregate_axis_mixture_effects(compound_axis_burdens; mixture_effect_model::String = "axis_toxic_unit_sum")
    X_A, X_M, X_G, X_R = 0.0, 0.0, 0.0, 0.0
    E_A, E_M, E_G, E_R = 0.0, 0.0, 0.0, 0.0

    if mixture_effect_model == "axis_toxic_unit_sum"
        for r in compound_axis_burdens
            x_A = Float64(r.burden_assimilation)
            x_M = Float64(r.burden_maintenance)
            x_G = Float64(r.burden_growth)
            x_R = Float64(r.burden_reproduction)
            
            for x in (x_A, x_M, x_G, x_R)
                if !isfinite(x) || x < 0.0
                    throw(ArgumentError("Burden values must be finite and >= 0. Got: $x"))
                end
            end
            
            X_A += x_A
            X_M += x_M
            X_G += x_G
            X_R += x_R
        end
        
        E_A = X_A / (1.0 + X_A)
        E_M = X_M / (1.0 + X_M)
        E_G = X_G / (1.0 + X_G)
        E_R = X_R / (1.0 + X_R)
        
    elseif mixture_effect_model == "independent_action_axis_effects"
        prod_A, prod_M, prod_G, prod_R = 1.0, 1.0, 1.0, 1.0
        
        for r in compound_axis_burdens
            x_A = Float64(r.burden_assimilation)
            x_M = Float64(r.burden_maintenance)
            x_G = Float64(r.burden_growth)
            x_R = Float64(r.burden_reproduction)
            
            for x in (x_A, x_M, x_G, x_R)
                if !isfinite(x) || x < 0.0
                    throw(ArgumentError("Burden values must be finite and >= 0. Got: $x"))
                end
            end
            
            X_A += x_A
            X_M += x_M
            X_G += x_G
            X_R += x_R
            
            prod_A *= (1.0 - (x_A / (1.0 + x_A)))
            prod_M *= (1.0 - (x_M / (1.0 + x_M)))
            prod_G *= (1.0 - (x_G / (1.0 + x_G)))
            prod_R *= (1.0 - (x_R / (1.0 + x_R)))
        end
        
        E_A = 1.0 - prod_A
        E_M = 1.0 - prod_M
        E_G = 1.0 - prod_G
        E_R = 1.0 - prod_R
        
    elseif mixture_effect_model == "grouped_ca_then_ia_axis_effects"
        # Group by axis and effect_code
        group_X_A = Dict{String, Float64}()
        group_X_M = Dict{String, Float64}()
        group_X_G = Dict{String, Float64}()
        group_X_R = Dict{String, Float64}()
        
        for r in compound_axis_burdens
            x_A = Float64(r.burden_assimilation)
            x_M = Float64(r.burden_maintenance)
            x_G = Float64(r.burden_growth)
            x_R = Float64(r.burden_reproduction)
            
            for x in (x_A, x_M, x_G, x_R)
                if !isfinite(x) || x < 0.0
                    throw(ArgumentError("Burden values must be finite and >= 0. Got: $x"))
                end
            end
            
            X_A += x_A
            X_M += x_M
            X_G += x_G
            X_R += x_R
            
            ec = hasproperty(r, :effect_code) ? r.effect_code : "unknown_effect_code"
            if ismissing(ec) || ec === nothing || ec == ""
                ec = "unknown_effect_code"
            end
            
            group_X_A[ec] = get(group_X_A, ec, 0.0) + x_A
            group_X_M[ec] = get(group_X_M, ec, 0.0) + x_M
            group_X_G[ec] = get(group_X_G, ec, 0.0) + x_G
            group_X_R[ec] = get(group_X_R, ec, 0.0) + x_R
        end
        
        prod_A, prod_M, prod_G, prod_R = 1.0, 1.0, 1.0, 1.0
        
        for (ec, X) in group_X_A
            prod_A *= (1.0 - (X / (1.0 + X)))
        end
        for (ec, X) in group_X_M
            prod_M *= (1.0 - (X / (1.0 + X)))
        end
        for (ec, X) in group_X_G
            prod_G *= (1.0 - (X / (1.0 + X)))
        end
        for (ec, X) in group_X_R
            prod_R *= (1.0 - (X / (1.0 + X)))
        end
        
        E_A = 1.0 - prod_A
        E_M = 1.0 - prod_M
        E_G = 1.0 - prod_G
        E_R = 1.0 - prod_R
        
    else
        throw(ArgumentError("Unknown mixture_effect_model: $mixture_effect_model"))
    end
    
    return (
        mixture_effect_model = mixture_effect_model,
        X_assimilation = X_A,
        X_maintenance = X_M,
        X_growth = X_G,
        X_reproduction = X_R,
        E_assimilation = E_A,
        E_maintenance = E_M,
        E_growth = E_G,
        E_reproduction = E_R
    )
end
