function write_ascii_grid(filename::String, grid::Matrix{Float64}; xllcorner=0.0, yllcorner=0.0, cellsize=1.0, nodata=-9999.0)
    nrows, ncols = size(grid)
    open(filename, "w") do f
        println(f, "ncols         ", ncols)
        println(f, "nrows         ", nrows)
        println(f, "xllcorner     ", xllcorner)
        println(f, "yllcorner     ", yllcorner)
        println(f, "cellsize      ", cellsize)
        println(f, "NODATA_value  ", nodata)

        for i in 1:nrows
            for j in 1:ncols
                val = isnan(grid[i, j]) ? nodata : grid[i, j]
                print(f, val, (j == ncols ? "" : " "))
            end
            println(f)
        end
    end
end

function read_ascii_grid(filename::String)
    lines = readlines(filename)

    ncols = parse(Int, split(lines[1])[2])
    nrows = parse(Int, split(lines[2])[2])
    xllcorner = parse(Float64, split(lines[3])[2])
    yllcorner = parse(Float64, split(lines[4])[2])
    cellsize = parse(Float64, split(lines[5])[2])
    nodata = parse(Float64, split(lines[6])[2])

    grid = zeros(Float64, nrows, ncols)
    for i in 1:nrows
        row_vals = split(lines[6+i])
        for j in 1:ncols
            val = parse(Float64, row_vals[j])
            if val == nodata
                grid[i, j] = NaN
            else
                grid[i, j] = val
            end
        end
    end

    return grid
end
