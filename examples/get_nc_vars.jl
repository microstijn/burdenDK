using NCDatasets

files = [
    raw"C:\Users\peete074\Downloads\pathogen_monthlyAvg_1980_2019.nc",
    raw"C:\Users\peete074\Downloads\organic_monthlyAvg_1980_2019.nc"
]

for file in files
    println("\n==============================")
    println("File: ", file)
    println("==============================")

    ds = NCDataset(file, "r")

    println("\nDimensions:")
    for (name, dim) in ds.dim
        println("  ", name, " => ", length(dim))
    end

    println("\nVariables:")
    for (name, var) in ds
        println("  ", name)
        println("    size: ", size(var))
        println("    dims: ", dimnames(var))
        if haskey(var.attrib, "units")
            println("    units: ", var.attrib["units"])
        end
        if haskey(var.attrib, "long_name")
            println("    long_name: ", var.attrib["long_name"])
        end
    end

    close(ds)
end