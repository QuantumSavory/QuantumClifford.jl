
#=============================================================================#
# These values originate from a package extension, hence requiring this query.
const KAExt = Base.get_extension(QuantumClifford, :QuantumCliffordKAExt)

if benchmark_primary_axis
    const primary_axes = [axis for axis in instances(KAExt.PrimaryAxis)]
else
    const primary_axes = [KAExt.default_primary_axis]
end

if benchmark_phases
    const phases = [true, false]
else
    const phases = [KAExt.default_phases]
end

if benchmark_block_size
    const block_sizes = [32, 64, 128, 256, 512]
else
    const block_sizes = [KAExt.default_block_size]
end

if benchmark_batch_size
    const batch_sizes = [1, 4, 8, 16, 32, 64]
else
    const batch_sizes = [KAExt.default_batch_size]
end
#=============================================================================#
