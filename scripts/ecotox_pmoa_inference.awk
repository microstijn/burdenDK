# ecotox_pmoa_inference.awk — DIAGNOSTIC (negative result, kept on record).
#
# Attempts to infer a stressor's physiological mode of action (pMoA) from the raw ECOTOX
# dump by asking which DEB-axis-relevant effect code responds at the lowest concentration.
# Effect-code -> DEB axis: FDB->assimilation, GRO/DVP/MPH->growth, REP->reproduction,
# PHY->maintenance, MOR->survival.
#
# CONCLUSION (see scripts/ecotox_pmoa_inference.jl output and
# data/external/ecotox_pmoa_inference_results.txt): this does NOT identify pMoA. The DEB
# kappa-rule cascade makes REPRODUCTION the most-sensitive endpoint for almost any pMoA
# (reproduction is funded by what's left after maintenance + growth), and the taxon mix
# differs across effect codes. So "most-sensitive endpoint" != "process taxed". Use FITTED
# DEBtox pMoAs instead (docs/notes/pmoa_evidence_to_gather.md, Tier 1). Kept as a documented
# negative so the route is not rebuilt.
#
# Usage (run from the raw ECOTOX dir; raw dump is gitignored):
#   cd data/ecotox/ecotox_ascii_03_12_2026
#   awk -f ../../../scripts/ecotox_pmoa_inference.awk tests.txt results.txt \
#       > ../../external/ecotox_pmoa_records.csv
# Then: julia +release --project=. scripts/ecotox_pmoa_inference.jl   (from repo root)
#
# CAS set below: Cd, Cu, chlorpyrifos, atrazine, bisphenol-A. Edit `want` to change.
BEGIN{FS="|"; OFS=","; split("7440439 7440508 2921882 1912249 80057", a, " "); for(i in a) want[a[i]]=1}
function uf(u,  s){ s=u; sub(/^AI /,"",s);
  if(s=="ug/L"||s=="ppb"||s=="ug/l") return 1;
  if(s=="mg/L"||s=="ppm"||s=="ug/ml"||s=="mg/l") return 1000;
  if(s=="ng/L"||s=="ng/l") return 0.001;
  if(s=="mg/ml"||s=="g/L"||s=="g/l") return 1000000;
  return -1; }
function axisof(e){ sub(/\/.*/,"",e);
  if(e=="FDB") return "assimilation";
  if(e=="GRO"||e=="DVP"||e=="MPH") return "growth";
  if(e=="REP") return "reproduction";
  if(e=="PHY") return "maintenance";
  if(e=="MOR") return "survival";
  return ""; }
FILENAME~/tests.txt/ { if($3 in want){ tcas[$1]=$3 } next }
FILENAME~/results.txt/ {
  tid=$2; if(!(tid in tcas)) next;
  ax=axisof($22); if(ax=="") next;
  ep=$19; cls="";
  if(ep=="NOEC"||ep=="LOEC"||ep=="EC10"||ep=="NOEL"||ep=="MATC") cls="chronic";
  else if(ep=="EC50"||ep=="LC50") cls="acute";
  else next;
  f=uf($43); if(f<0) next;
  v=$38; if(v=="" || v !~ /^[0-9.]+$/) next; val=v*f; if(val<=0) next;
  print tcas[tid], ax, cls, val
}
