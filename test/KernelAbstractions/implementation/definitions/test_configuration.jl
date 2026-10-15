
#=============================================================================#
# Small counts for encoding issues, large counts for race conditions.
const qubit_counts = [
    31, 32, 33, 63, 64, 65, 127, 128, 129,
    64 * 1023, 64 * 1024, 64 * 1025, 64 * 2047, 64 * 2048, 64 * 2049
    ]
# The tests are for correctness, not for device memory limits.
const max_rows = 1024

# Keep it reasonable so that local testing remains accessible.
const max_rounds = 16
# Certain operations are quite expensive and so should run fewer iterations.
const min_rounds = min(4, max_rounds)
#=============================================================================#
