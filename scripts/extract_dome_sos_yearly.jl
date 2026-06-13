# Station-YEAR panel for the Stress-on-Stress temporal analysis (power-boost + erosion-over-
# time test). Base Julia only -> julia +release scripts/extract_dome_sos_yearly.jl
# Reads the gitignored raw DOME biota download, writes data/external/sos_dome_ukcemp_yearly.csv
# (one row per SURVT station-year). Contaminant/size values are matched to the NEAREST year
# with data at the same station within +/-2 years (tissue burden changes slowly), so all SoS
# station-years are usable; a `contam_lag` column records the year offset used.
#
# Companion to the station-level test (scripts/extract_dome_sos.jl). Same source/units.

using Statistics

const RAW = joinpath(@__DIR__, "..", "data", "external", "ospar_cemp2024_biota_data.csv")
const OUT = joinpath(@__DIR__, "..", "data", "external", "sos_dome_ukcemp_yearly.csv")
const WINDOW = 2   # nearest-year match window (years)

function splitcsv(line)
    out = String[]; buf = IOBuffer(); inq = false
    for c in line
        if c == '"'; inq = !inq
        elseif c == ',' && !inq; push!(out, String(take!(buf)))
        else; print(buf, c); end
    end
    push!(out, String(take!(buf))); return out
end
unq(s) = strip(s, ['"', ' '])
num(s) = (s = unq(s); isempty(s) ? NaN : something(tryparse(Float64, s), NaN))

const C_STATION = 5; const C_NAME = 6; const C_YEAR = 12; const C_SPECIES = 17
const C_DET = 22; const C_CONC = 26; const C_LNMEA = 29; const C_DRYWT = 30
const METALS = ["CD","CU","HG","PB","ZN"]; const PAH23 = ["NAP","ACNLE","ACNE","FLE","PA","ANT"]
const CONG = ["CB28","CB52","CB101","CB118","CB138","CB153","CB180"]
const CAPTURE = Set(vcat(["SURVT","SCB7"], METALS, PAH23, CONG))

# (station,year) => determinand => [values]; plus size keyed the same
acc = Dict{Tuple{String,Int}, Dict{String, Vector{Float64}}}()
name = Dict{String,String}()
survt_st = Set{String}()

open(RAW) do io
    readline(io)
    for line in eachline(io)
        f = splitcsv(line); length(f) < C_DRYWT && continue
        unq(f[C_DET]) == "SURVT" && push!(survt_st, unq(f[C_STATION]))
    end
end
open(RAW) do io
    readline(io)
    for line in eachline(io)
        f = splitcsv(line); length(f) < C_DRYWT && continue
        st = unq(f[C_STATION]); (st in survt_st) || continue
        unq(f[C_SPECIES]) == "Mytilus edulis" || continue
        yr = tryparse(Int, unq(f[C_YEAR])); yr === nothing && continue
        get!(name, st, unq(f[C_NAME]))
        d = get!(acc, (st, yr), Dict{String, Vector{Float64}}())
        det = unq(f[C_DET])
        det in CAPTURE && (v = num(f[C_CONC]); isfinite(v) && push!(get!(d, det, Float64[]), v))
        l = num(f[C_LNMEA]); isfinite(l) && push!(get!(d, "LNMEA", Float64[]), l)
        w = num(f[C_DRYWT]); isfinite(w) && push!(get!(d, "DRYWT", Float64[]), w)
    end
end

med(d, k) = (haskey(d, k) && !isempty(d[k])) ? median(d[k]) : NaN
# nearest-year median of determinand `k` at station `st` to year `y0` within WINDOW
function nearest(st, y0, k)
    best = NaN; bestlag = 99
    for ((s, y), d) in acc
        s == st || continue
        abs(y - y0) <= WINDOW || continue
        m = med(d, k); isfinite(m) || continue
        if abs(y - y0) < bestlag; bestlag = abs(y - y0); best = m; end
    end
    return best, (isfinite(best) ? (bestlag == 99 ? 0 : bestlag) : 99)
end
fmt(x) = isfinite(x) ? string(round(x, sigdigits = 5)) : "nd"

# PCB per (station,year): reported SCB7 else sum of congener medians
function pcb(d)
    s = med(d, "SCB7"); isfinite(s) && return s
    c = [med(d, g) for g in CONG]; all(isfinite, c) ? sum(c) : NaN
end

cols = vcat(["station_code","station_name","year","survt_median_d","survt_n","contam_lag",
             "LNMEA","DRYWT"], METALS, PAH23, ["SCB7"])
rows = String[]
for ((st, yr), d) in acc
    haskey(d, "SURVT") && !isempty(d["SURVT"]) || continue   # only SoS station-years
    # nearest-year contaminants + size + the lag actually used (max over fields)
    getn(k) = nearest(st, yr, k)
    metvals = [getn(m) for m in METALS]; pahvals = [getn(p) for p in PAH23]
    # PCB nearest-year: search nearest year that yields a finite pcb()
    pcbval = NaN; pcblag = 99
    for ((s, y), dd) in acc
        s == st && abs(y - yr) <= WINDOW || continue
        p = pcb(dd); isfinite(p) || continue
        if abs(y - yr) < pcblag; pcblag = abs(y - yr); pcbval = p; end
    end
    lags = filter(l -> l < 99, vcat([v[2] for v in metvals], [v[2] for v in pahvals], pcblag == 99 ? Int[] : [pcblag]))
    lag = isempty(lags) ? 99 : maximum(lags)
    ln, _ = getn("LNMEA"); dw, _ = getn("DRYWT")
    row = vcat([st, "\"$(name[st])\"", string(yr), fmt(med(d, "SURVT")), string(length(d["SURVT"])),
                string(lag), fmt(ln), fmt(dw)],
               [fmt(v[1]) for v in metvals], [fmt(v[1]) for v in pahvals], [fmt(pcbval)])
    push!(rows, join(row, ","))
end
sort!(rows)

open(OUT, "w") do io
    println(io, "# Station-YEAR Stress-on-Stress panel, Mytilus edulis, ICES DOME 2024 OSPAR CEMP.")
    println(io, "# One row per SURVT station-year; contaminants/size matched to nearest year with")
    println(io, "# data at the same station within +/-$(WINDOW) yr (contam_lag = yr offset used, 99=none).")
    println(io, "# survt_median_d = median survival-in-air days. Generated by scripts/extract_dome_sos_yearly.jl.")
    println(io, join(cols, ","))
    for r in rows; println(io, r); end
end
println("wrote ", OUT, "  (", length(rows), " SURVT station-years)")
for line in eachline(OUT); startswith(line, "#") || println(line); end
