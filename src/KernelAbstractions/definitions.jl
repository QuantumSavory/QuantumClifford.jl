
#=============================================================================#
# Reasonable size that is generally ideal for most vendors and usecases.
const default_block_size = 256
# Ameliorate overhead and enhance performance by doing more work per thread.
const default_batch_size = 32

# Most ideal type for shared reductions due to avoiding bank conflicts.
const DeviceUnsigned = UInt32

# Keeps function definitions succinct.
const DevicePauliOperator = PauliOperator{T_P, T_XZ} where {
	T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray
	}
const DeviceTableau = Tableau{T_P, T_XZ} where {
	T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray
	}
# This is a bit redundant but it keeps the REPL expansions more readable.
const DeviceUnionStabilizer = Union{
	Stabilizer{T}, MixedStabilizer{T}
	} where {
		T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
		T <: Tableau{T_P, T_XZ}
		}
const DeviceUnionDestabilizer = Union{
	Destabilizer{T}, MixedDestabilizer{T}
	} where {
		T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
		T <: Tableau{T_P, T_XZ}
		}
const DeviceAbstractStabilizer = Union{
	Stabilizer{T}, MixedStabilizer{T}, Destabilizer{T}, MixedDestabilizer{T}
	} where {
		T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
		T <: Tableau{T_P, T_XZ}
		}
#=============================================================================#
