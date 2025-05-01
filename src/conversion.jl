using QuantumClifford
  
struct ClassicalParity{N} <: AbstractOperation
    "The indices of the classical bits to be xor-ed"
    bits::NTuple{N,Int}
    "The index of the classical bit that will store the results"
    store::Int
end

struct ClassicalParityRelative{N}
    bits::NTuple{N, Int}
end
struct Repeat
    body::Vector
    times:: Int
end

struct sMRZAbsolute <: AbstractOperation
    qubit::Int
    bit::Int
    sMRZ(q, args...) = if q<=0 throw(NoZeroQubit) else new(q,args...) end
end

struct PauliFrame
    measurements::Matrix{Int}
end


function apply!(frame::PauliFrame, xor::ClassicalParity)
    for row in 1:size(frame.measurements, 1)
        value = frame.measurements[row, xor.bits[1] + 1]
        for i in xor.bits[2:end]
            value ⊻= frame.measurements[row, i + 1]
        end
        frame.measurements[row, xor.store + 1] = value
    end
    return frame
end


function apply!(frame::PauliFrame, r::Repeat)
    for i in 1:r.times
        for instructions in r.body
            abs_bits = ntuple(j -> i - 1 + instr.bits[j], length(instructions.bits))
            abs_instr = ClassicalParity(abs_bits, i)
            apply!(frame, abs_instr)
        end
    end
end


function convert(r::Repeat)::Vector{ClassicalParity}
    absolute_instructions = ClassicalParity[]
    for i in r.times
        for instructions in r.body
            if instructions isa ClassicalParity
                absolute_bits = ntuple(j->i -1+instructions.bits[j], length(instructions.bits))#translates relative position of the bits into an absolute index
                push!(absolute_instructions, ClassicalParity(absolute_bits, i-1+100))#creates a classicalparity check and stores the the parity check to the list absolute_instructions list
            end
        end
    end
    return absolute_instructions
end

