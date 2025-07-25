# (La)TeX hates SVG but the Plots package has issues with transparent PDFs.
const format = "svg"

# BenchmarkTools parameters.
# Evaluations per sample point.
const evals = 5
# Maximum number of samples.
const samples = 10_000
# Maximum runtime for each sample group.
const seconds = 300

# By definition, (unsigned) char is the smallest addressable unit of memory.
const MiB = 1024 * 1024 * count_zeros(zero(Cuchar))
# Avoid consuming too many resources, 1 GiB is plenty.
const n_MiB = [2^i for i = 1:10]
# TODO: Keep these or remove them now that a good default has been set?
const batch_sizes = [1, 4, 8, 16, 32, 64]

# These values are inaccessible since they originate from a package extension.
const default_block_size = 256
const default_batch_size = 32
