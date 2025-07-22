# By definition, (unsigned) char is the smallest addressable unit of memory.
const MiB = 1024 * 1024 * count_zeros(zero(Cuchar))
# Avoid consuming too many resources, 1 GiB is plenty.
const n_MiB = [2^i for i = 1:10]

const batch_sizes = [1, 4, 8, 16, 32, 64]

# These values are inaccessible since they originate from a package extension.
const default_block_size = 256
const default_batch_size = 32
