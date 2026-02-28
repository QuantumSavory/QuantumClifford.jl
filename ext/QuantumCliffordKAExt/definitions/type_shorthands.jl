
#=============================================================================#
const DevicePauliOperator = PauliOperator{T_P, T_XZ} where {
    T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray
    }

const DeviceTableau = Tableau{T_P, T_XZ} where {
    T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray
    }

const DeviceUnionTableau = Union{
    T, Stabilizer{T}, MixedStabilizer{T}, Destabilizer{T}, MixedDestabilizer{T}
    } where {
        T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
        T <: Tableau{T_P, T_XZ}
        }

# Utilised to specialise dispatch.
const DeviceUnionDestabilizer = Union{
    Destabilizer{T}, MixedDestabilizer{T}
    } where {
        T_P <: AbstractGPUArray, T_XZ <: AbstractGPUArray,
        T <: Tableau{T_P, T_XZ}
        }
#=============================================================================#
