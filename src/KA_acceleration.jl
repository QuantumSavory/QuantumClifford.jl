using GPUArraysCore: AbstractGPUArray
import KernelAbstractions as KA
# Resolves issue due to KA comparing against the literal Symbol("@Const").
using KernelAbstractions: @Const
import Atomix

const DeviceUnsigned = UInt32
# Cannot allocate shared memory if this is not known at compile time.
const default_block_size = 256

# For use whenever KA.@index(Global, NTuple) is unavailable.
@inline global_index(block_index, block_dim, thread_index) =
    (block_index .- one(eltype(block_index))) .* block_dim .+ thread_index

# CAUTION: Requires block_size == length(buffer) == prod(KA.@groupsize()).
# CAUTION: Requires unsafe_indices = true if num_active_threads < block_size.
# TODO: Fix function signature once support is either confirmed or denied.
# TODO: Revisit once warp level primitives are supported.
@inline function shared_memory_reduce!(
    f, buffer::AbstractArray{T}, value::T, index, ::Val{block_size}
    ) where {T, block_size}

    @inbounds buffer[index] = value

    # This branch is a power of 2, take the quick route.
    if count_ones(block_size) == one(block_size)

    # This is absolutely hideous but only for-loop unrolling is supported.
    # TODO: Confirm whether this unrolls as desired.
    KA.Extras.@unroll for bit = one(block_size) : trailing_zeros(block_size)
        stride = KA.@uniform ((2)^(trailing_zeros(block_size) - bit))
        # The call to KA.@synchronize is ALL or NOTHING.
        KA.@synchronize()
        if index <= stride
            @inbounds value = convert(T, f(value, buffer[index + stride]))
            @inbounds buffer[index] = value
        end
    end

    else

    current = KA.@uniform (block_size)
    # TODO: Unroll this branch should it become possible.
    while current > one(current)
        # The call to KA.@synchronize is ALL or NOTHING.
        KA.@synchronize()
        # This splits into even/odd steps. Should target SALU.
        if iseven(current)
            current = KA.@uniform (current >> one(current))
            if index <= current
                @inbounds value = convert(T, f(value, buffer[index + current]))
                @inbounds buffer[index] = value
            end
        else
            current = KA.@uniform ((current >> one(current)) + one(current))
            # The strict inequality is intentional.
            if index < current
                @inbounds value = convert(T, f(value, buffer[index + current]))
                @inbounds buffer[index] = value
            end
        end
    end

    end
    return value

end

# CAUTION: Keep in mind that mutable = order_right_left ? right : left.
KA.@kernel inbounds = true unsafe_indices = true function kernel_mul_ordered!(
    mutable_xzs, mutable_phases, @Const(const_xzs), @Const(const_phases),
    ::Val{order_right_left}, ::Val{phases}, ::Val{block_size}
    ) where {order_right_left, block_size, phases}

    # unsafe_indices is required for shared_memory_reduce, do this manually.
    global_pos = global_index(
        KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
        )
    i = global_pos[1]
    j_mutable = size(mutable_xzs, 2) > 1 ? global_pos[2] : 1
    j_const = size(const_xzs, 2) > 1 ? global_pos[2] : 1
    z_offset = KA.@uniform (size(const_xzs, 1) >> 1)
    if phases
        low = zero(eltype(const_xzs))
        high = zero(eltype(const_xzs))
    end

    if i <= z_offset
        if order_right_left
            x_left = const_xzs[i, j_const]
            z_left = const_xzs[i + z_offset, j_const]
            x_right = mutable_xzs[i, j_mutable]
            z_right = mutable_xzs[i + z_offset, j_mutable]
        else
            x_left = mutable_xzs[i, j_mutable]
            z_left = mutable_xzs[i + z_offset, j_mutable]
            x_right = const_xzs[i, j_const]
            z_right = const_xzs[i + z_offset, j_const]
        end

        x_new = xor(x_left, x_right)
        z_new = xor(z_left, z_right)
        if phases
            xl_zr = x_left & z_right
            zl_xr = z_left & x_right
            low = xor(xl_zr, zl_xr)
            high = xor(x_new, z_new, xl_zr) & low
        end

        mutable_xzs[i, j_mutable] = x_new
        mutable_xzs[i + z_offset, j_mutable] = z_new
    end

    if phases
        value::DeviceUnsigned = (count_ones(high) << 1) & 0x3
        value += count_ones(low) & 0x3
        buffer = KA.@localmem DeviceUnsigned (block_size,)
        index = KA.@index(Local, Linear)
        value = shared_memory_reduce!(
            +, buffer, value, index, Val(block_size)
            )
        if index == one(index)
            if i == one(i)
                value += const_phases[j_const] & 0x3
            end
            # This assumes bits(mutable_phases) >= bits(DeviceUnsigned).
            Atomix.@atomic mutable_phases[j_mutable] +=
                convert(eltype(mutable_phases), value)
            Atomix.@atomic mutable_phases[j_mutable] &=
                convert(eltype(mutable_phases), 0x3)
        end
    end

end

function mul_ordered!(
    mutable_xzs::AbstractGPUArray{T},
    mutable_phases::AbstractGPUArray{<: Unsigned},
    const_xzs::AbstractGPUArray{T},
    const_phases::AbstractGPUArray{<: Unsigned};
    order_right_left::Val{RL}, phases::Val{B} = Val(true)
    ) where {T <: Unsigned, RL, B}

    dim_x = size(const_xzs, 1) >> 1
    dim_y = max(size(mutable_xzs, 2), size(const_xzs, 2))
    # This resolves a race condition if doing mul(::Tableau, n, n).
    length(const_phases) == 1 && (const_phases = copy(const_phases))
    kernel = kernel_mul_ordered!(KA.get_backend(const_xzs), (default_block_size, 1))
    kernel(
        mutable_xzs, mutable_phases, const_xzs, const_phases,
        order_right_left, phases, Val(default_block_size);
        workgroupsize = (default_block_size, 1), ndrange = (dim_x, dim_y)
        )

end
