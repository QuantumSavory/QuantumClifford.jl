# Small sizes for encoding issues, large sizes for race conditions.
const test_sizes = [
    31, 32, 33, 63, 64, 65, 127, 128, 129,
    64 * 1023, 64 * 1024, 64 * 1025, 64 * 2047, 64 * 2048, 64 * 2049
    ]
# The tests are for correctness, not for device memory limits.
const max_rows = 1024
# Keep it reasonable so that local testing remains accessible.
const round_count = 16
# Correctness should be independent of parameter values.
# The omission of the const specifier is intentional, overriden in OpenCL.
block_sizes = rand(1:256, round_count)
const batch_sizes = rand(1:256, round_count)
