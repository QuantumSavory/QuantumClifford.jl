const InnerGPUType = UInt32;
const InnerCPUType = UInt64;
# todo. make the conversion types consiting with this...
# also define a cpu default type?!

to_gpu(array::AbstractArray) = CuArray(array)
to_cpu(array::AbstractArray) = Array(array);

# changes the type to InnerGPUType or InnerCPUType as well as copying to cpu/gpu
to_gpu_reinterpret(array::AbstractArray) = CuArray(reinterpret(InnerGPUType, collect(array)))
to_cpu_reinterpret(array::AbstractArray) = Array(reinterpret(InnerCPUType, collect(array)))

# maybe change the format of storing the data in gpu array 
# so that it is more convinient to work with them on gpu?
# todo later add some type checking to avoid copying (or throw error) if the data is already on gpu/cpu
to_gpu(tab::QuantumClifford.Tableau) =
    QuantumClifford.Tableau(to_gpu(tab.phases), tab.nqubits, to_gpu_reinterpret(tab.xzs))

to_cpu(tab::QuantumClifford.Tableau) =
    QuantumClifford.Tableau(to_cpu(tab.phases), tab.nqubits, to_cpu_reinterpret(tab.xzs))

to_gpu(pauli::QuantumClifford.PauliOperator) =
    QuantumClifford.PauliOperator(to_gpu(pauli.phase), pauli.nqubits, to_gpu_reinterpret(pauli.xz))

to_cpu(pauli::QuantumClifford.PauliOperator) =
    QuantumClifford.PauliOperator(to_cpu(pauli.phase), pauli.nqubits, to_cpu_reinterpret(pauli.xz))

to_gpu(stabilizer::QuantumClifford.Stabilizer) =
    Stabilizer(to_gpu(tab(stabilizer)))

to_cpu(stabilizer::QuantumClifford.Stabilizer) = 
    Stabilizer(to_cpu(tab(stabilizer)))

to_gpu(pauli_frame::QuantumClifford.PauliFrame) =
    QuantumClifford.PauliFrame(to_gpu(pauli_frame.frame), to_gpu(pauli_frame.measurements))

to_cpu(pauli_frame::QuantumClifford.PauliFrame) =
    QuantumClifford.PauliFrame(to_cpu(pauli_frame.frame), to_cpu(pauli_frame.measurements))
