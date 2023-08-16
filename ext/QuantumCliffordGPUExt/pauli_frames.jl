function apply!(f::PauliFrameGPU{T}, op::QuantumClifford.AbstractCliffordOperator) where {T <: Unsigned}
    _apply!(f.frame, op; phases=Val(false))
    return f
end

function apply_sMZ_kernel!(xzs::DeviceMatrix{Tme},
                          measurements::DeviceMatrix{Bool},
                          op::sMZ,
                          ibig::Int,
                          ismallm::Tme,
                          rows::Int) where {Tme <: Unsigned} 
    f = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    if f > rows
        return nothing
    end    
    should_flip = !iszero(xzs[ibig,f] & ismallm)
    measurements[f,op.bit] = should_flip
    return nothing
end

function apply!(frame::PauliFrameGPU{T}, op::QuantumClifford.sMZ) where {T <: Unsigned} # TODO sMX, sMY
    op.bit == 0 && return frame
    i = op.qubit
    xzs = frame.frame.tab.xzs
    lowbit = T(1)
    ibig = QuantumClifford._div(T,i-1)+1
    ismall = QuantumClifford._mod(T,i-1)
    ismallm = lowbit<<(ismall)
    (@run_cuda apply_sMZ_kernel!(xzs, frame.measurements, op, ibig, ismallm, length(frame)) length(frame))
    return frame
end

function apply_sMRZ_kernel!(xzs::DeviceMatrix{Tme},
                          measurements::DeviceMatrix{Bool},
                          op::QuantumClifford.sMRZ,
                          ibig::Int, # todo change to Int
                          ismallm::Tme,
                          rows::Int) where {Tme <: Unsigned} 
    f = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    if f > rows
        return nothing
    end
    if op.bit != 0 # not good practice. how to replace 0?
        should_flip = !iszero(xzs[ibig,f] & ismallm)
        measurements[f,op.bit] = should_flip
    end
    xzs[ibig,f] &= ~ismallm
    rand(Bool) && (xzs[end÷2+ibig,f] ⊻= ismallm)
    return nothing
end

function apply!(frame::PauliFrameGPU{T}, op::QuantumClifford.sMRZ) where {T <: Unsigned} # TODO sMRX, sMRY
    i = op.qubit
    xzs = frame.frame.tab.xzs
    lowbit = T(1)
    ibig = QuantumClifford._div(T,i-1)+1
    ismall = QuantumClifford._mod(T,i-1)
    ismallm = lowbit<<(ismall)
    (@run_cuda apply_sMRZ_kernel!(xzs, frame.measurements, op, ibig, ismallm, length(frame)) length(frame))
    return frame
end



##################################################################################

# function partition_circuit(::Type{T}, circ) where {T <: Unsigned}
#     # time complexity: O(gates * qubits/64) for UInt64

#     circ = if eltype(circ) <: QuantumClifford.CompactifiedGate
#         circ
#     else
#         compactify_circuit(circ)
#     end

#     qmax = maximum((maximum(affectedqubits(g)) for g in circ))
#     columns_cnt = QuantumClifford._div(T, qmax) + 1

#     partitions = []
#     partition_idx = []
#     last_column_occupier = [0 for _ in 1:columns_cnt] # last_column_occupier[i] = group id of the last gate that occupies indice i

#     for g in circ
#         column_indices = [QuantumClifford.getbigindex(T, qbit) for qbit in affectedqubits(g)]
#         my_partition_idx = maximum(last_column_occupier[column_indices]) + 1
#         push!(partition_idx, my_partition_idx)
#         if my_partition_idx > length(partitions)
#             push!(partitions, [g])
#         else
#             push!(partitions[my_partition_idx], g)
#         end
#         last_column_occupier[column_indices] .= my_partition_idx
#     end
#     partitions
# end

# function partition_circuit(circ)
#     # time complexity: O(gates * qubits/64) for UInt64

#     circ = if eltype(circ) <: QuantumClifford.CompactifiedGate
#         circ
#     else
#         compactify_circuit(circ)
#     end

#     qmax = maximum((maximum(affectedqubits(g)) for g in circ))

#     partitions = []
#     partition_idx = []
#     last_column_occupier = [0 for _ in 1:qmax] # last_column_occupier[i] = group id of the last gate that occupies indice i

#     for g in circ
#         indices = [x for x in affectedqubits(g)]
#         my_partition_idx = maximum(last_column_occupier[indices]) + 1
#         push!(partition_idx, my_partition_idx)
#         if my_partition_idx > length(partitions)
#             push!(partitions, [g])
#         else
#             push!(partitions[my_partition_idx], g)
#         end
#         last_column_occupier[indices] .= my_partition_idx
#     end
#     partitions
# end

# function pftrajectories(state::PauliFrameGPU{T}, circuit) where {T <: Unsigned}
#     for gates in partition_circuit(circuit)
#         # for gate in gates
#         #     CUDA.@sync apply!(state, gate)
#         # end

#         @sync begin
#             for gate in gates
#                 @async apply!(state, gate)
#             end
#         end
#     end
#     return state
# end
