
#=============================================================================#

#==============================================================================
BENCHMARK TOOLS
==============================================================================#

# CAUTION: Functions mutate their arguments, induces disparity between runs.
# Evaluations per sample point.
const evals = 1

# Maximum number of sample points.
const samples = 2^14

# Maximum runtime for each trial.
const seconds = 60

#==============================================================================
SAMPLE EXTRAPOLATION
==============================================================================#

# Maximum sampling period before before being aborted and extrapolated instead.
const extrapolation_threshold = seconds << 1

# Whether to include the aborted run in the data set used for extrapolation.
const include_threshold_point = false

# In the absence of sufficient data points for a fit, perform O(n^k) scaling.
const host_permit_simple_scaling = true
const device_permit_simple_scaling = true

#==============================================================================
PROBLEM SIZE
==============================================================================#

# By definition, (unsigned) char is the smallest addressable unit of memory.
const MiB = 1024 * 1024 * count_zeros(zero(Cuchar))

# Avoid consuming too many resources, 1 GiB is plenty.
const sizes_MiB = [2^i for i in 1 : 10]

#==============================================================================
TUNING PARAMETERS
==============================================================================#

const benchmark_primary_axis = true

const benchmark_phases = true

# TODO: Enable this by default once the POCL code generation bugs are fixed.
const benchmark_block_size = false

const benchmark_batch_size = true
#=============================================================================#
