
#=============================================================================#
# Maintains compatibility with the main package unless explicitly specified.
const default_multiplication_order = multiplication_order_left

# Customary choice to position X pauli operators before their Z counterparts.
const default_pauli_preferance = pauli_preferance_x

# Potentially boosts cache hits and reduces atomic contention.
const default_primary_axis = primary_axis_rows

# Strict correctness unless there is an explicit opt-out.
const default_phases = true

# Reasonable size that is generally ideal for most vendors and use cases.
# TODO: Modify this to DeviceUnsigned once CUDA bugs are resolved.
const default_block_size = 256

# Ameliorate overhead and enhance performance by doing more work per thread.
const default_batch_size = UInt8(32)
#=============================================================================#
