using LinearAlgebra: Adjoint


const PhaseInnerType = UInt8
const MeasurementInnerType = Bool

CUDAValue{T} = CuArray{T, 0, CUDA.Mem.DeviceBuffer} where {T}
CUDAVector{T} = CuArray{T, 1, CUDA.Mem.DeviceBuffer} where {T}
CUDAMatrix{T} = CuArray{T, 2, CUDA.Mem.DeviceBuffer} where {T}
CUDAParams = [CUDAValue, CUDAVector, CUDAMatrix]

AdjCUDAValue{T} = CuArray{T, 0, CUDA.Mem.DeviceBuffer} where {T} # do not adjoint value
AdjCUDAVector{T} = CuArray{T, 1, CUDA.Mem.DeviceBuffer} where {T} # do not adjoint vector
AdjCUDAMatrix{T} = Adjoint{T, CuArray{T, 2, CUDA.Mem.DeviceBuffer}} where {T}
AdjCUDAParams = [AdjCUDAValue, AdjCUDAVector, AdjCUDAMatrix]

DeviceValue{T} = Union{CuDeviceArray{T, 0, 1}, Adjoint{T, CuDeviceArray{T, 0, 1}}} where {T}
DeviceVector{T} = Union{CuDeviceArray{T, 1, 1}, Adjoint{T, CuDeviceArray{T, 1, 1}}} where {T}
DeviceMatrix{T} = Union{CuDeviceArray{T, 2, 1}, Adjoint{T, CuDeviceArray{T, 2, 1}}} where {T}


function getTableauGPU(GPUValue, GPUVector, GPUMatrix)
    TableauGPU{T} = QuantumClifford.Tableau{Tₚᵥ, Tₘ} where {T <: Unsigned, Tₚᵥ <: GPUVector{PhaseInnerType}, Tₘ <: GPUMatrix{T}}
end
function getStabilizerGPU(GPUValue, GPUVector, GPUMatrix, GPUTableau)
    StabilizerGPU{T} = QuantumClifford.Stabilizer{<:GPUTableau{T}} where {T <: Unsigned}
end
function getPauliOperatorGPU(GPUValue, GPUVector, GPUMatrix, GPUTableau)
    PauliOperatorGPU{T} = QuantumClifford.PauliOperator{Tₚ, Tᵥ} where {T <: Unsigned, Tₚ<:GPUValue{PhaseInnerType}, Tᵥ<:GPUVector{T}}
end

# todo. type definition here is stronger than the code in pauliframes.jl  this will cause serious problems
# especially because its not obvious whether TFrame is Tableau or Stabilizer in pauliframes.jl
# and we are assuming that TMeasurement is made of booleans
function getPauliFrameGPU(GPUValue, GPUVector, GPUMatrix, GPUTableau)
    StabilizerGPU{T} = getStabilizerGPU(GPUValue, GPUVector, GPUMatrix, GPUTableau){T} where {T <: Unsigned}
    PauliFrameGPU{T} = QuantumClifford.PauliFrame{TFrame, TMeasurement} where {T <: Unsigned, TFrame <: StabilizerGPU{T}, TMeasurement <: GPUMatrix{MeasurementInnerType}}
end

const TableauCUDA{T} = getTableauGPU(CUDAParams...){T} where {T <: Unsigned}
const StabilizerCUDA{T} = getStabilizerGPU(CUDAParams..., TableauCUDA){T} where {T <: Unsigned}
const PauliOperatorCUDA{T} = getPauliOperatorGPU(CUDAParams..., TableauCUDA){T} where {T <: Unsigned}
const PauliFrameCUDA{T} = getPauliFrameGPU(CUDAParams..., TableauCUDA){T} where {T <: Unsigned}

const TableauAdj{T} = getTableauGPU(AdjCUDAParams...){T} where {T <: Unsigned}
const StabilizerAdj{T} = getStabilizerGPU(CUDAParams..., TableauAdj){T} where {T <: Unsigned}
const PauliOperatorAdj{T} = getPauliOperatorGPU(CUDAParams..., TableauAdj){T} where {T <: Unsigned}
const PauliFrameAdj{T} = getPauliFrameGPU(CUDAParams..., TableauAdj){T} where {T <: Unsigned}

const TableauGPU{T} = Union{TableauCUDA{T}, TableauAdj{T}} where {T <: Unsigned}
const StabilizerGPU{T} = Union{StabilizerCUDA{T}, StabilizerAdj{T}} where {T <: Unsigned}
const PauliOperatorGPU{T} = Union{PauliOperatorCUDA{T}, PauliOperatorAdj{T}} where {T <: Unsigned}
const PauliFrameGPU{T} = Union{PauliFrameCUDA{T}, PauliFrameAdj{T}} where {T <: Unsigned}
