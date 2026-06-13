# rubach2010_rate_axis.jl — the rate-axis (toxicokinetic) negative control (synthesis §7b).
#
# Pairs AmP k_M against measured chlorpyrifos elimination rate constants k_out
# (Rubach et al. 2010, ETC 29:2225, Table 2 — 15 freshwater arthropods, radiotracer, non-DEB).
# k_out values + fresh weight (Table 1) are transcribed below (small, fixed dataset).
# Result: AmP overlap = 6 species (+1 genus proxy); k_M <-> k_out weak/wrong-sign, n.s.
# Caveat: k_out is the TOXICOKINETIC elimination rate, not the thesis-relevant TOXICODYNAMIC
# recovery rate k_r — so a null here is the expected off-target result (see scoping memo §5d).
#
# Run:  julia +release --project=. scripts/rubach2010_rate_axis.jl
using Statistics
# species, AmP k_M (/day), Rubach k_out (/day), fresh weight (mg), AmP L_m (cm), EC50 (ug/L), tier
rows = [
("Asellus_aquaticus",   0.11903, 0.185,  11.5,  0.15173, 6.16,  "species"),
("Cloeon_dipterum",     0.030835,0.196,   7.7,  0.29154, 0.763, "species"),
("Culex_pipiens",       0.70975, 0.0240,  3.4,  0.016333,NaN,   "species"),
("Daphnia_magna",       0.27273, 0.546,  13.5,  0.1516,  0.484, "species"),
("Gammarus_pulex",      0.072624,0.398,  19.1,  0.1114,  0.379, "species"),
("Notonecta_maculata",  0.078092,0.152, 133.0,  0.36275, 9.07,  "species"),
("Procambarus_clarkii", 1.0419,  0.0860,2794.2, 1.0486,  1.7,   "genus-proxy"),
]
function report(label, idx)
    kM=[rows[i][2] for i in idx]; kout=[rows[i][3] for i in idx]; fw=[rows[i][4] for i in idx]
    lkM,lko,lfw = log10.(kM),log10.(kout),log10.(fw)
    resid(y,x)= y .- ((mean(y)-cov(x,y)/var(x)*mean(x)) .+ (cov(x,y)/var(x)).*x)
    println("\n=== $label  (n=$(length(idx))) ===")
    println("  cor(log k_M,  log kout) = ", round(cor(lkM,lko),digits=3))
    println("  cor(log k_M,  log fw)   = ", round(cor(lkM,lfw),digits=3))
    println("  PARTIAL cor(k_M,kout|fw)= ", round(cor(resid(lkM,lfw),resid(lko,lfw)),digits=3))
end
report("Solid species-level only", 1:6)
report("Incl. Procambarus genus-proxy", 1:7)
open("data/external/rubach2010_kM_kout_paired.csv","w") do io
    println(io, "species,amp_k_M_perday,rubach_kout_perday,rubach_freshweight_mg,amp_L_m_cm,rubach_EC50_48h_ugL,match_tier")
    for r in rows
        println(io, join([r[1],r[2],r[3],r[4],r[5], isnan(r[6]) ? "" : r[6], r[7]], ","))
    end
end
