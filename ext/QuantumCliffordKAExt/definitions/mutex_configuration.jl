
#=============================================================================#
# TODO: Eliminate these once KernelAbstractions becomes more feature complete.
AbstractMutex = AbstractArray{DeviceUnsigned, 0}

# These should be confined to an enum but atomics require type conversion.
const mutex_state_locked = zero(DeviceUnsigned)
const mutex_state_unlocked = one(DeviceUnsigned)

# Defines the mappings for the compare-and-swap(CAS)/@atomicreplace call.
const mutex_exchange_lock = mutex_state_unlocked => mutex_state_locked
const mutex_exchange_unlock = mutex_state_locked => mutex_state_unlocked
#=============================================================================#
