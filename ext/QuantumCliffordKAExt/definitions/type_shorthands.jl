
#=============================================================================#
# Keeps the function definitions succinct.
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
