# extract_ecotox_multimoa.awk -- per-species effect concentrations for a MULTI-MoA chemical
# panel, from the raw ECOTOX ASCII dump. Generalises extract_ecotox_acute.awk: instead of a
# hard-coded 4-CAS AChE set, it reads a CAS -> (moa_class, primary_axis) panel from a CSV and
# emits one row per usable (LC50/EC50/IC50/NOEC/LOEC) water-concentration result.
#
# This feeds the across-axis capacity-weighting COVERAGE SCAN (does the same AmP-resident species
# set recur across >=2 DISCRIMINATING MoA classes -- maintenance/growth vs reproduction -- so the
# kappa-driven axis weighting is identifiable beyond body size?).
#
# Usage (file ORDER matters -- panel, then species map, then tests, then results):
#   cd data/ecotox/ecotox_ascii_03_12_2026
#   awk -f ../../../scripts/extract_ecotox_multimoa.awk \
#       ../../ecotox_multimoa_panel.csv validation/species.txt tests.txt results.txt \
#       > ../../external/ecotox_multimoa_extract.csv
#
# Emits: cas, moa_class, primary_axis, species_number, latin, ecotox_group,
#        endpoint, effect, obs_duration, obs_unit, value_ugL
# Raw ECOTOX dump is gitignored; the derived extract is a scan intermediate.
BEGIN{ FS="|"; OFS=",";
  print "cas,moa_class,primary_axis,species_number,latin,ecotox_group,endpoint,effect,obs_duration,obs_unit,value_ugL" }
# unit -> factor to ug/L (water conc only); -1 if not a usable water conc
function uf(u,  s){ s=u; sub(/^AI /,"",s);
  if(s=="ug/L"||s=="ppb"||s=="ug/l") return 1;
  if(s=="mg/L"||s=="ppm"||s=="ug/ml"||s=="mg/l") return 1000;
  if(s=="ng/L"||s=="ng/l") return 0.001;
  if(s=="mg/ml"||s=="g/L"||s=="g/l") return 1000000;
  return -1; }
# file1: panel CSV (comma-delimited; no pipes, so $0 is the whole line)
FILENAME ~ /panel/ {
  n=split($0,a,",");
  if(n<4 || a[1]=="cas" || a[1]=="") next;
  cls[a[1]]=a[3]; axis[a[1]]=a[4];
  next }
# file2: species map  (species_number -> latin_name, ecotox_group)
FILENAME ~ /species/ { lat[$1]=$3; grp[$1]=$15; next }
# file3: tests -> test_id -> cas|species  (only panel CAS)
FILENAME ~ /tests/ {
  if($3 in cls){ tcas[$1]=$3; tsp[$1]=$18 }
  next }
# file4: results
FILENAME ~ /results/ {
  tid=$2; if(!(tid in tcas)) next;
  ep=$19;
  if(!(ep=="LC50"||ep=="EC50"||ep=="IC50"||ep=="NOEC"||ep=="LOEC")) next;
  ef=$22;
  f=uf($43); if(f<0) next;
  v=$38; if(v=="" || v !~ /^[0-9.]+$/) next;
  val=v*f; if(val<=0) next;
  cas=tcas[tid]; sp=tsp[tid]; l=lat[sp]; gsub(/ /,"_",l);
  print cas, cls[cas], axis[cas], sp, l, grp[sp], ep, ef, $12, $17, val
}
