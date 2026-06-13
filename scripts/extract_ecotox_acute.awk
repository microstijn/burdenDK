# extract_ecotox_acute.awk — species-level acute lethality from the raw ECOTOX ASCII dump.
#
# Pulls acute (24-96 h) LC50/MOR + EC50/IMM endpoints in water-concentration units,
# converts to ug/L, for a hard-coded set of chemical CAS (Baas & Kooijman 2015's four
# AChE inhibitors: chlorpyrifos 2921882, malathion 121755, carbofuran 1563662,
# carbaryl 63252). Emits: cas, species_number, latin_underscore, ecotox_group, value_ugL.
#
# Usage (file ORDER matters — species map, then tests, then results):
#   cd data/ecotox/ecotox_ascii_03_12_2026
#   awk -f scripts/extract_ecotox_acute.awk validation/species.txt tests.txt results.txt \
#       > data/external/ecotox_acute_4chem.csv
#
# To target other chemicals, edit the CAS set in the tests.txt and results.txt rules.
# Raw ECOTOX dump is gitignored; the derived CSV is committed.
BEGIN{FS="|"; OFS=","}
# unit -> factor to ug/L (water conc only); returns -1 if not a usable water conc
function uf(u,  s){ s=u; sub(/^AI /,"",s);
  if(s=="ug/L"||s=="ppb"||s=="ug/l") return 1;
  if(s=="mg/L"||s=="ppm"||s=="ug/ml"||s=="mg/l") return 1000;
  if(s=="ng/L"||s=="ng/l") return 0.001;
  if(s=="mg/ml"||s=="g/L"||s=="g/l") return 1000000;
  return -1; }
# file1: species map  (species_number -> latin_name, ecotox_group)
FNR==NR && FILENAME~/species.txt/ { lat[$1]=$3; grp[$1]=$15; next }
# file2: tests -> test_id -> cas|species  (only our 4 CAS)
FILENAME~/tests.txt/ {
  if($3==2921882||$3==121755||$3==1563662||$3==63252){ tcas[$1]=$3; tsp[$1]=$18 }
  next }
# file3: results
FILENAME~/results.txt/ {
  tid=$2; if(!(tid in tcas)) next;
  ep=$19; ef=$22;
  ok = ((ep=="LC50"&&ef=="MOR")||(ep=="EC50"&&(ef=="IMM"||ef=="MOR")));
  if(!ok) next;
  # acute duration filter: 24-96 h, or 1-4 d
  dv=$12+0; ac=0;
  if($17=="h" && dv>=24 && dv<=96) ac=1;
  else if($17=="d" && dv>=1 && dv<=4) ac=1;
  if(!ac) next;
  f=uf($43); if(f<0) next;
  v=$38; if(v=="" || v !~ /^[0-9.]+$/) next;
  val=v*f; if(val<=0) next;
  sp=tsp[tid]; l=lat[sp]; gsub(/ /,"_",l);
  print tcas[tid], sp, l, grp[sp], val
}
