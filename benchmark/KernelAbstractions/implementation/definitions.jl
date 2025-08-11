# (La)TeX hates SVG but the Plots package has issues with transparent PDFs.
const format = "svg"

# BenchmarkTools parameters.
# Evaluations per sample point.
const evals = 16
# Maximum number of samples.
const samples = 2^10
# Maximum runtime for each sample group.
const seconds = 60

# By definition, (unsigned) char is the smallest addressable unit of memory.
const MiB = 1024 * 1024 * count_zeros(zero(Cuchar))
# Avoid consuming too many resources, 1 GiB is plenty.
const n_MiB = [2^i for i = 1:10]
# TODO: Keep these or remove them now that a good default has been set?
const batch_sizes = [1, 4, 8, 16, 32, 64]

# These values originate from a package extension, hence the query.
const KAExt = Base.get_extension(QuantumClifford, :QuantumCliffordKAExt)
const default_phases = KAExt.default_phases
const default_primary_axis = KAExt.default_primary_axis
const default_block_size = KAExt.default_block_size
const default_batch_size = KAExt.default_batch_size
