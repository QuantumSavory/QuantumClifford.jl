# Small sizes for encoding issues, large sizes for race conditions.
const test_sizes = [
	31, 32, 33, 63, 64, 65, 127, 128, 129,
	64 * 1023, 64 * 1024, 64 * 1025, 64 * 2047, 64 * 2048, 64 * 2049
	]
# Keep it reasonable so that local testing remains accessible.
const cycle_range = 1:16
