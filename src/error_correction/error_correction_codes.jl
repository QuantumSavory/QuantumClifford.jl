#currently in the works

module Codes
abstract type Code end

struct Shor9 <: Code
end

struct Toric <: Code
    l::Int
end

struct CSSfromClassical <: Code
    classicalcode::Matrix{Bool}
end

"""docs"""
function rate end

function logicalqubits end

function physicalqubits end

function codedistance end

function H(code::Shor9) end

function QuantumClifford.MixedDestabilizer(code::Shor9)
    ...
end

...

function rate(code::Shor9) return 1//9 end

... decoding algorithms and how to structure them
end # module

