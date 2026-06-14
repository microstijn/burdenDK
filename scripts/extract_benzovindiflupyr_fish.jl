# ===========================================================================
# Extract the 5-fish benzovindiflupyr survival data from the Nickisch Born Gericke et al. 2022
# (ETC 41(7):1732-1741, doi:10.1002/etc.5348) Supporting Information (etc5348-sup-0004) into a
# committed long-format CSV, and compute the per-species acute LC50 (the cross-species sensitivity
# metric for the capacity probe). Raw SI is gitignored; this writes the derived CSV.
#
# The .xlsx was unzipped to a temp dir (an xlsx is a zip of XML). Sheets 4-8 are the per-species
# survival tables (cols: SampleID, replicate, time, conc, Nsurv). We use the ACUTE tests.
#
#   julia +release --project=. scripts/extract_benzovindiflupyr_fish.jl
# Writes data/external/benzovindiflupyr_fish_survival.csv (+ prints per-species acute LC50).
# ===========================================================================

using Statistics, Printf

const DUMP = "C:/Users/peete074/AppData/Local/Temp/xlsx_dump"
const SHEETS = ["sheet4.xml" => "Cyprinodon_variegatus", "sheet5.xml" => "Cyprinus_carpio",
                "sheet6.xml" => "Lepomis_macrochirus",   "sheet7.xml" => "Oncorhynchus_mykiss",
                "sheet8.xml" => "Pimephales_promelas"]

# --- shared strings (column A = SampleID is a string index) ---
function shared_strings(path)
    xml = read(path, String); out = String[]
    for m in eachmatch(r"<si>(.*?)</si>"s, xml)
        push!(out, join([t.captures[1] for t in eachmatch(r"<t[^>]*>([^<]*)</t>", m.captures[1])], ""))
    end
    return out
end

# --- parse one worksheet into rows of (sample, replicate, time, conc, nsurv) ---
function parse_sheet(path, strings)
    xml = read(path, String)
    cells = Dict{Tuple{Int,Char},Any}()
    for m in eachmatch(r"<c r=\"([A-Z]+)(\d+)\"([^>]*)>\s*<v>([^<]*)</v>", xml)
        col = m.captures[1][1]; row = parse(Int, m.captures[2]); attrs = m.captures[3]; v = m.captures[4]
        cells[(row, col)] = occursin("t=\"s\"", attrs) ? strings[parse(Int, v) + 1] : parse(Float64, v)
    end
    rows = NamedTuple[]
    maxrow = maximum(r for (r, _) in keys(cells))
    for r in 2:maxrow
        haskey(cells, (r, 'A')) || continue
        push!(rows, (sample = String(cells[(r, 'A')]), rep = cells[(r, 'B')],
                     time = cells[(r, 'C')], conc = cells[(r, 'D')], nsurv = cells[(r, 'E')]))
    end
    return rows
end

# --- per-species acute LC50 at the final acute timepoint (log-linear interp on conc-response) ---
function acute_lc50(rows)
    ac = filter(r -> startswith(r.sample, "acute"), rows)
    isempty(ac) && return (NaN, NaN, 0)
    tmax = maximum(r.time for r in ac)
    concs = sort(unique(r.conc for r in ac))
    frac = Float64[]
    for c in concs
        n0 = sum(r.nsurv for r in ac if r.conc == c && r.time == 0.0; init = 0.0)
        nf = sum(r.nsurv for r in ac if r.conc == c && r.time == tmax; init = 0.0)
        push!(frac, n0 > 0 ? nf / n0 : NaN)
    end
    # interpolate concentration at 50% survival (log conc), between the bracketing points
    lc50 = NaN
    for i in 1:length(concs)-1
        if isfinite(frac[i]) && isfinite(frac[i+1]) && frac[i] >= 0.5 >= frac[i+1] && concs[i] > 0
            t = (frac[i] - 0.5) / (frac[i] - frac[i+1])
            lc50 = exp(log(concs[i]) + t * (log(concs[i+1]) - log(concs[i])))
            break
        end
    end
    return (lc50, tmax, length(concs))
end

function main()
    strings = shared_strings(joinpath(DUMP, "xl/sharedStrings.xml"))
    out = open(joinpath(@__DIR__, "..", "data", "external", "benzovindiflupyr_fish_survival.csv"), "w")
    println(out, "# Benzovindiflupyr fish survival (long), from Nickisch Born Gericke et al. 2022,")
    println(out, "# ETC 41(7):1732-1741 (doi:10.1002/etc.5348) SI etc5348-sup-0004. Underlying acute")
    println(out, "# data: Ashauer et al. 2013. conc_ugL = ug a.s./L; time_d days; nsurv = # alive.")
    println(out, "species,sample,replicate,time_d,conc_ugL,nsurv")
    println("\n=== per-species ACUTE LC50 (benzovindiflupyr, final acute timepoint) ===")
    @printf("  %-24s %10s  %6s  %5s\n", "species", "LC50 (ug/L)", "t_d", "nconc")
    for (file, sp) in SHEETS
        rows = parse_sheet(joinpath(DUMP, "xl/worksheets", file), strings)
        for r in rows
            @printf(out, "%s,%s,%g,%g,%g,%g\n", sp, replace(r.sample, "," => ";"), r.rep, r.time, r.conc, r.nsurv)
        end
        lc50, tmax, nc = acute_lc50(rows)
        @printf("  %-24s %10.3f  %6.1f  %5d\n", sp, lc50, tmax, nc)
    end
    close(out)
    println("\nwrote data/external/benzovindiflupyr_fish_survival.csv")
end

main()
