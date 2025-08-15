
#=============================================================================#
@inline function create_mutex(backend::KA.Backend)
    output = KA.allocate(backend, DeviceUnsigned)
    fill!(output, mutex_state_unlocked)
    return output
end

@inline function reset_mutex!(mutex::AbstractMutex)
    fill!(mutex, mutex_state_unlocked)
    return mutex
end

# WARNING: Atomix has an API for atomic load/store, yet support is missing.
# WARNING: KernelAbstractions does not provide any sort of memory fence.

#==============================================================================
SINFUL IMPLEMENTATION WROUGHT ABOUT BY THE FOLLY OF MANKIND
==============================================================================#

# TODO: Overhaul this entirely once support becomes available.
@inline @generated function lock_mutex!(
    mutex::M, data::AbstractArray...
    )::Nothing where {M <: AbstractMutex}

    if length(data) > 0x0
        clause = :(
            @atomicreplace :release :acquire mutex[0x1] mutex_exchange_lock
            )
        # CAUTION: Atomics forces the necessary memory synchronisation barrier.
        @inbounds nil = zero(eltype(data[0x1]))
        fence = :(
            # This is just a fancy atomic NOOP.
            @atomicreplace :release :acquire data[0x1][0x1] $nil => $nil;
            )
    end
    for n in 0x2 : length(data)
        @inbounds nil = zero(eltype(data[n]))
        fence = :(
            $fence;
            @atomicreplace :release :acquire data[$n][0x1] $nil => $nil;
            )
    end
    return :(
        @inbounds begin;
            while true;
                ($clause).success && break;
            end;
            $fence;
        end;
        return nothing;
        )

end

# TODO: Overhaul this entirely once support becomes available.
@inline @generated function unlock_mutex!(
    mutex::M, data::AbstractArray...
    )::Nothing where {M <: AbstractMutex}

    if length(data) > 0x0
        # CAUTION: Atomics forces the necessary memory synchronisation barrier.
        fence = :(
            temp_1 = data[0x1][0x1];
            # This is just a fancy atomic NOOP. Always succeeds and releases.
            @atomicreplace :release :acquire data[0x1][0x1] temp_1 => temp_1;
            )
    end
    for n in 0x2 : length(data)
        sym = Symbol(:temp_, n)
        fence = :(
            $fence;
            $sym = data[$n][0x1];
            @atomicreplace :release :acquire data[$n][0x1] $sym => $sym;
            )
    end
    return :(
        @inbounds begin;
            $fence;
            # This will always succeed.
            @atomicreplace :release :acquire mutex[0x1] mutex_exchange_unlock;
        end;
        return nothing;
        )

end
#=============================================================================#
