files = [
    "src/deb_axes.jl",
    "src/mode_of_action.jl",
    "src/exposure_filters.jl",
    "src/moa_deb_mapping.jl",
]

function extract_struct_blocks(file)
    println("\n# ", file)

    lines = readlines(file)
    inside_struct = false
    depth = 0
    buffer = String[]
    start_line = 0

    for (i, line) in enumerate(lines)
        s = strip(line)

        if !inside_struct &&
           occursin(r"^(@kwdef\s+)?(mutable\s+)?struct\s+", s)

            inside_struct = true
            depth = 1
            start_line = i
            empty!(buffer)
            push!(buffer, line)
            continue
        end

        if inside_struct
            if i != start_line
                push!(buffer, line)
            end

            # Struct blocks can contain inner constructors with function/end,
            # so we count major block keywords conservatively.
            if occursin(r"^(@kwdef\s+)?(mutable\s+)?struct\s+", s) && i != start_line
                depth += 1
            elseif s == "end"
                depth -= 1
            end

            if depth == 0
                println("\nL$start_line:")
                println(join(buffer, "\n"))
                inside_struct = false
            end
        end
    end
end

for file in files
    extract_struct_blocks(file)
end