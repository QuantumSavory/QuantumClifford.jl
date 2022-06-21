
module Codes # requires end
abstract type Code end

struct Shor_code <: Code
end

"""documents"""

function rate end

function logicalqubits end

function physicalqubits end

function codedistance end