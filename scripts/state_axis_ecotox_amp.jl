# state_axis_ecotox_amp.jl — the state-axis negative control (synthesis §7b).
#
# Tests Baas & Kooijman (2015): does AmP specific somatic maintenance ([p_M], or the rate
# k_M = [p_M]/[E_G]) predict acute sensitivity (LC50/EC50) across species, BEYOND body size?
# Pairs per-species acute LC50 (from data/external/ecotox_acute_4chem.csv, produced by
# scripts/extract_ecotox_acute.awk) with AmP k_M / p_M / L_m, per chemical and pooled.
#
# Result (n=310): the raw maintenance->sensitivity signal replicates (rho ~ -0.27, all 4
# chemicals) but NULLS under a body-size control (partial ~ -0.03). A robust null, not underpower.
#
# Run:  julia +release --project=. scripts/state_axis_ecotox_amp.jl
import JSON
using Statistics
lib = JSON.parsefile("data/AmP_Species_Library.json")
amp = Dict{String,NTuple{3,Float64}}()  # latin_underscore -> (k_M, p_M, L_m)
for (sp,v) in lib
    a=get(v,"auxiliary_metrics",Dict())
    haskey(a,"k_M")&&haskey(a,"p_M")&&haskey(a,"L_m") || continue
    amp[sp]=(float(a["k_M"]),float(a["p_M"]),float(a["L_m"]))
end
chemname=Dict("2921882"=>"Chlorpyrifos","121755"=>"Malathion","1563662"=>"Carbofuran","63252"=>"Carbaryl")
grp=Dict{Tuple{String,String},Vector{Float64}}()
for ln in eachline("data/external/ecotox_acute_4chem.csv")
    f=split(ln,","); length(f)>=5||continue
    cas=f[1]; sp=f[3]; val=tryparse(Float64,f[5]); val===nothing&&continue
    push!(get!(grp,(cas,sp),Float64[]), log10(val))
end
rows=[]  # (cas, sp, n, lc50_ugL, kM, pM, Lm)
for ((cas,sp),vs) in grp
    haskey(amp,sp)||continue
    kM,pM,Lm=amp[sp]
    push!(rows,(cas,sp,length(vs),10^median(vs),kM,pM,Lm))
end
println("Matched (chemical×species) pairs to AmP: ", length(rows))
open("data/external/state_axis_ecotox_amp_paired.csv","w") do io
    println(io,"chemical_cas,chemical,species,n_records,LC50_ugL,amp_k_M,amp_p_M,amp_L_m")
    for r in sort(rows;by=x->(x[1],x[2]))
        println(io,join([r[1],chemname[r[1]],r[2],r[3],round(r[4],sigdigits=5),r[5],r[6],r[7]],","))
    end
end
function pcor(y,x,z)  # partial cor of y,x controlling z
    rx=x .- (mean(x).+(cov(z,x)/var(z)).*(z.-mean(z)))
    ry=y .- (mean(y).+(cov(z,y)/var(z)).*(z.-mean(z)))
    cor(rx,ry)
end
function analyze(label, sel)
    sub=[r for r in rows if sel(r)]
    length(sub)<4 && (println("\n[$label] n=",length(sub)," (too few)"); return)
    y =[log10(r[4]) for r in sub]; kM=[log10(r[5]) for r in sub]
    pM=[log10(r[6]) for r in sub]; Lm=[log10(r[7]) for r in sub]
    println("\n=== $label  (n=",length(sub),") ===")
    println("  cor(logLC50, log k_M)        = ",round(cor(y,kM),digits=3))
    println("  cor(logLC50, log p_M [B&K])  = ",round(cor(y,pM),digits=3))
    println("  cor(logLC50, log L_m size)   = ",round(cor(y,Lm),digits=3))
    println("  partial cor(LC50,k_M | size) = ",round(pcor(y,kM,Lm),digits=3))
    println("  partial cor(LC50,p_M | size) = ",round(pcor(y,pM,Lm),digits=3))
end
for (cas,nm) in chemname; analyze(nm, r->r[1]==cas); end
analyze("POOLED (all 4 chemicals)", r->true)
println("\nNote: LC50 high = LESS sensitive. B&K predict maintenance↑ → sensitivity↑ → LC50↓ (negative cor).")
