
 

function parse_file(filepath::String)

    operations = Vector{Tuple{String, Vector{Int}}}()
    for raw in eachline(filepath)
        line = strip(split(raw, "#")[1])
        isempty(line) && continue
        parts = split(line)
        cmd = parts[1]
        qubits = parse.(Int, parts[2:end])
        push!(operations, (cmd, qubits))
    end
    return operations
end

function infer_nqubits(path::String)
    max_idx = -1
    for raw in eachline(path)
        
        line = strip(split(raw, "#")[1])
        isempty(line) && continue

        parts = split(line)
       
        qs = parse.(Int, parts[2:end])
        if !isempty(qs)
            max_idx = max(max_idx, maximum(qs))
        end
    end
    return max_idx + 1
end


const STIM2QC = Dict(
    "H"    => (1,qs -> sHadamard(qs[1])),
    "X"    => (1,qs -> sX(qs[1])),
    "Z"    => (1,qs -> sZ(qs[1])),
    "CNOT" => (2,qs -> sCNOT(qs[1], qs[2])),
    "CX"   => (2,qs -> sCNOT(qs[1], qs[2])),
    "S"    => (1,qs -> sPhase(qs[1])),
)

function converter(filepath::String)
    ops = parse_file(filepath)
    n   = infer_nqubits(filepath)
    st  = one(Stabilizer, n, basis=:Z)

    

    for (cmd, qs0) in ops
        if haskey(STIM2QC, cmd)
            num, build = STIM2QC[cmd]

            # Ensure total qubits is a multiple of arity
            if length(qs0) % num != 0
                error("Invalid qubit list for $cmd: length $(length(qs0)) isn’t a multiple of arity $num")
            end

            # Chunk and apply
            for i in 1:num:length(qs0)
                slice = qs0[i : i+num-1]
                qs    = slice .+ 1           # Stim is 0-based → 1-based indexing
                op    = build(qs)
                apply!(st, op)
            end
        end
    end

    return st
end

filepath = "/Users/rohan/Desktop/CS/lab_work/stim_test/my_circuit.stim"
converter(filepath)
