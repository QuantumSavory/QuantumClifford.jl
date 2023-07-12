TableauGPU{T} = QuantumClifford.Tableau{Tzv, Tm} where {T <: Unsigned, Tzv <: CuArray{UInt8}, Tm <: CuArray{T, 2}}
StabilizerGPU{T} = QuantumClifford.Stabilizer{Tab} where {T <: Unsigned, Tzv <: CuArray{UInt8}, Tm <: CuArray{T, 2}, Tab <: QuantumClifford.Tableau{Tzv, Tm}}
