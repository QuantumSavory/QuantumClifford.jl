
#=============================================================================#
# Correctness should be independent of the tuning parameter values.
# They originate from a package extension, hence requiring this query.
const KAExt = Base.get_extension(QuantumClifford, :QuantumCliffordKAExt)

const primary_axes = rand(instances(KAExt.PrimaryAxis), max_rounds)
# The omission of the const specifier is intentional, overridden in OpenCL.
# TODO: Revisit this once the POCL code generation issues are resolved.
block_sizes = rand(Base.OneTo(KAExt.default_block_size), max_rounds)
const batch_sizes = rand(Base.OneTo(KAExt.default_batch_size), max_rounds)
#=============================================================================#
