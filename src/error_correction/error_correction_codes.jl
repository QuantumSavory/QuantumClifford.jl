#currently in the works

module Codes
abstract type Code end

struct Shor_code <: Code
end

"""documents"""

function rate end

function logicalqubits end

function physicalqubits end

function codedistance end

function H(code::Shor_code) end

function QuantumClifford.MixedDestabilizer(code::Shor_code)

    Stabilizer([P"X",P"Y",P"Z"])

    Stabilizer([0x2, 0x0],
                Bool[0 1;
                    1 0],
                Bool[0 -i;
                    i 1],
                Bool[1 0;
                    0 -1])
end

...

function rate(code::Shor_code) return 1//9 end

... decoding algorithms and how to structure them
end # module

