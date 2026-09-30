
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
    mutex::AbstractMutex, data::AbstractArray...
    )::Nothing

    if length(data) > 0
        swap = :(@atomicreplace :release :acquire mutex[1] mutex_exchange_lock)
        # CAUTION: Atomics force the necessary memory synchronisation barrier.
        @inbounds nil = zero(eltype(data[1]))
        # This is just a fancy atomic NOOP.
        fence = :(@atomicreplace :release :acquire data[1][1] $nil => $nil;)
    end
    for n in 2 : length(data)
        @inbounds nil = zero(eltype(data[n]))
        fence = :(
            $fence;
            @atomicreplace :release :acquire data[$n][1] $nil => $nil;
            )
    end
    return :(
        @inbounds begin;
            while true;
                ($swap).success && break;
            end;
            $fence;
        end;
        return nothing;
        )

end

# TODO: Overhaul this entirely once support becomes available.
@inline @generated function unlock_mutex!(
    mutex::AbstractMutex, data::AbstractArray...
    )::Nothing

    if length(data) > 0
        # CAUTION: Atomics force the necessary memory synchronisation barrier.
        fence = :(
            temp_1 = data[1][1];
            # This is just a fancy atomic NOOP. Always succeeds and releases.
            @atomicreplace :release :acquire data[1][1] temp_1 => temp_1;
            )
    end
    for n in 2 : length(data)
        sym = Symbol(:temp_, n)
        fence = :(
            $fence;
            $sym = data[$n][1];
            @atomicreplace :release :acquire data[$n][1] $sym => $sym;
            )
    end
    return :(
        @inbounds begin;
            $fence;
            # This will always succeed.
            @atomicreplace :release :acquire mutex[1] mutex_exchange_unlock;
        end;
        return nothing;
        )

end
#=============================================================================#
