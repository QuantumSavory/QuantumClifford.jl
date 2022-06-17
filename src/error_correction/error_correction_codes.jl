#currently in the works

module Codes
abstract type Code end

struct Shor_code <: Code
end

struct Toric_code <: Code
    l::Int
end

struct CSSfromClassical <: Code
    classicalcode::Matrix{Bool}
end

"""documents"""

function rate end

function logicalqubits end

function physicalqubits end

function codedistance end

function H(code::Shor_code) end

function QuantumClifford.MixedDestabilizer(code::Shor_code)
    ...
end

...

function rate(code::Shor_code) return 1//9 end

... decoding algorithms and how to structure them
end # module

