const InnerGPUType = UInt32;
const InnerCPUType = UInt64;


to_gpu(array::AbstractArray) = CuArray(array)
to_cpu(array::AbstractArray) = Array(array)

to_gpu(array::AbstractArray, tp::Type{T}) where {T <: Unsigned} = 
    CuArray(reinterpret(T, collect(array)))
to_cpu(array::AbstractArray, tp::Type{T}) where {T <: Unsigned} = 
    Array(reinterpret(T, collect(array)))

# todo later add some type checking to avoid copying (or throw error) if the data is already on gpu/cpu
to_gpu(tab::QuantumClifford.Tableau, tp::Type{T}=InnerGPUType) where {T <: Unsigned} =
    QuantumClifford.Tableau(to_gpu(tab.phases), tab.nqubits, to_gpu(tab.xzs, tp))

to_cpu(tab::QuantumClifford.Tableau, tp::Type{T}=InnerCPUType) where {T <: Unsigned} =
    QuantumClifford.Tableau(to_cpu(tab.phases), tab.nqubits, to_cpu(tab.xzs, tp))

to_gpu(pauli::QuantumClifford.PauliOperator, tp::Type{T}=InnerGPUType) where {T <: Unsigned} =
    QuantumClifford.PauliOperator(to_gpu(pauli.phase), pauli.nqubits, to_gpu(pauli.xz, tp))

to_cpu(pauli::QuantumClifford.PauliOperator, tp::Type{T}=InnerCPUType) where {T <: Unsigned} =
    QuantumClifford.PauliOperator(to_cpu(pauli.phase), pauli.nqubits, to_cpu(pauli.xz, tp))

to_gpu(stabilizer::QuantumClifford.Stabilizer, tp::Type{T}=InnerGPUType) where {T <: Unsigned} =
    Stabilizer(to_gpu(tab(stabilizer), tp))

to_cpu(stabilizer::QuantumClifford.Stabilizer, tp::Type{T}=InnerCPUType) where {T <: Unsigned} = 
    Stabilizer(to_cpu(tab(stabilizer), tp))

to_gpu(pauli_frame::QuantumClifford.PauliFrame, tp::Type{T}=InnerGPUType) where {T <: Unsigned} =
    QuantumClifford.PauliFrame(to_gpu(pauli_frame.frame, tp), to_gpu(pauli_frame.measurements))

to_cpu(pauli_frame::QuantumClifford.PauliFrame, tp::Type{T}=InnerCPUType) where {T <: Unsigned} =
    QuantumClifford.PauliFrame(to_cpu(pauli_frame.frame, tp), to_cpu(pauli_frame.measurements))
