# todo is there anyway to do this automatically so that if the code in QuantumClifford changes, we don't have to change this?!
# especially with UInt8 types

const TableauGPU{T} = QuantumClifford.Tableau{Tzv, Tm} where {T <: Unsigned, Tzv <: CuArray{UInt8, 1}, Tm <: CuArray{T, 2}}
const StabilizerGPU{T} = QuantumClifford.Stabilizer{<:TableauGPU{T}} where {T <: Unsigned}
const PauliOperatorGPU{T} = QuantumClifford.PauliOperator{Tz, Tv} where {T <: Unsigned, Tz<:CuArray{UInt8,0}, Tv<:CuArray{T, 1}}

# todo. type definition here is stronger than the code in pauliframes.jl  this will cause serious problems
# especially because its not obvious whether TFrame is Tableau or Stabilizer in pauliframes.jl
# and we are assuming that TMeasurement is made of booleans
const PauliFrameGPU{T} = QuantumClifford.PauliFrame{TFrame, TMeasurement} where {TFrame <: StabilizerGPU{T}, TMeasurement <: CuArray{Bool, 2}}

