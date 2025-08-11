
#=============================================================================#
# Maintains compatibility with the main package unless explicitly specified.
const default_multiplication_order = multiplication_order_left
# Strict correctness unless there is an explicit opt-out.
const default_phases = true
# Potentially boosts cache hits and reduces atomic contention.
const default_primary_axis = primary_axis_rows
# Reasonable size that is generally ideal for most vendors and use cases.
const default_block_size = 256
# Ameliorate overhead and enhance performance by doing more work per thread.
const default_batch_size = 32
# TODO: Eliminate this in favour of complete asynchronicity.
const default_scheduling_limit = 64
#=============================================================================#
