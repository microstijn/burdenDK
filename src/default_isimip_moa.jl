export default_isimip_moa_mapping
export default_isimip_interactions

function default_isimip_interactions()
    interactions = [zeros(7, 7) for _ in 1:8]

    # Mode names: thermal, oxygen, osmotic, immune, eutrophication, toxic, feeding, physical
    # Var indices: WT=1, BOD=2, TDS=3, FC=4, Nutrient=5, Chemical=6, Plastic=7

    # WT x BOD -> oxygen mode, coefficient 0.35
    interactions[2][1, 2] = 0.35

    # WT x TDS -> osmotic mode, coefficient 0.20
    interactions[3][1, 3] = 0.20

    # BOD x FC -> immune mode, coefficient 0.30
    interactions[4][2, 4] = 0.30

    # WT x Nutrient -> eutrophication mode, coefficient 0.30
    interactions[5][1, 5] = 0.30

    # BOD x Chemical -> toxic mode, coefficient 0.20
    interactions[6][2, 6] = 0.20

    return interactions
end

function default_isimip_moa_mapping(; clamp_unit=false)
    U = [
        # WT   BOD  TDS  FC   Nutr Chem Plast
        1.00  0.05 0.05 0.00 0.10 0.05 0.00;  # thermal
        0.25  1.00 0.00 0.10 0.30 0.05 0.00;  # oxygen
        0.00  0.00 1.00 0.00 0.00 0.05 0.00;  # osmotic
        0.00  0.10 0.00 1.00 0.00 0.10 0.00;  # immune
        0.20  0.10 0.00 0.00 1.00 0.00 0.00;  # eutrophication
        0.00  0.05 0.05 0.00 0.00 1.00 0.20;  # toxic
        0.05  0.20 0.00 0.10 0.20 0.20 1.00;  # feeding
        0.00  0.00 0.00 0.00 0.00 0.05 1.00   # physical
    ]

    return ModeOfActionMapping(
        U = U,
        interactions = default_isimip_interactions(),
        clamp_nonnegative = true,
        clamp_unit = clamp_unit
    )
end
