using Test
using TwoTimescaleResilience

@testset "Default ISIMIP MoA" begin
    mapping = default_isimip_moa_mapping()

    # 1. U has size 8 x 7.
    @test size(mapping.U) == (8, 7)

    # 2. interactions has length 8 and each matrix is 7 x 7.
    @test length(mapping.interactions) == 8
    for m in mapping.interactions
        @test size(m) == (7, 7)
    end

    # WT=1, BOD=2, TDS=3, FC=4, Nutr=5, Chem=6, Plast=7

    # 3. WT-only input primarily increases thermal mode.
    wt_only = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    m_wt = mode_of_action(wt_only, mapping)
    @test m_wt.thermal == 1.0
    @test m_wt.thermal > m_wt.oxygen
    @test m_wt.thermal > m_wt.eutrophication

    # 4. BOD-only input primarily increases oxygen mode.
    bod_only = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    m_bod = mode_of_action(bod_only, mapping)
    @test m_bod.oxygen == 1.0
    @test m_bod.oxygen > m_bod.thermal
    @test m_bod.oxygen > m_bod.immune

    # 5. TDS-only input primarily increases osmotic mode.
    tds_only = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    m_tds = mode_of_action(tds_only, mapping)
    @test m_tds.osmotic == 1.0
    @test m_tds.osmotic > m_tds.thermal

    # 6. FC-only input primarily increases immune mode.
    fc_only = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    m_fc = mode_of_action(fc_only, mapping)
    @test m_fc.immune == 1.0
    @test m_fc.immune > m_fc.oxygen

    # 7. Nutrient-only input primarily increases eutrophication mode.
    nutr_only = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    m_nutr = mode_of_action(nutr_only, mapping)
    @test m_nutr.eutrophication == 1.0
    @test m_nutr.eutrophication > m_nutr.thermal

    # 8. Chemical-only input primarily increases toxic mode.
    chem_only = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    m_chem = mode_of_action(chem_only, mapping)
    @test m_chem.toxic == 1.0
    @test m_chem.toxic > m_chem.immune

    # 9. Plastic-only input increases feeding and physical modes.
    plast_only = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    m_plast = mode_of_action(plast_only, mapping)
    @test m_plast.feeding == 1.0
    @test m_plast.physical == 1.0

    # 10. WT x BOD interaction increases oxygen mode above additive-only result.
    wt_bod = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    m_wt_bod = mode_of_action(wt_bod, mapping)
    # expected oxygen = 0.25*1.0 + 1.0*1.0 + 0.35*1.0*1.0 = 1.6
    @test m_wt_bod.oxygen ≈ 1.6
    @test m_wt_bod.oxygen > (m_wt.oxygen + m_bod.oxygen)
end
