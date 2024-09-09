function LPCode(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `LPCode` depends on the package `AbstractAlgebra` and `Hecke` but you have not installed or imported them yet. Immediately after you import `AbstractAlgebra` and `Hecke`, the `LPCode` will be available.")
    end
    return ext.LPCode(args...; kwargs...)
end

function two_block_group_algebra_codes(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `two_block_group_algebra_codes` depends on the package `AbstractAlgebra` and `Hecke` but you have not installed or imported them yet. Immediately after you import `AbstractAlgebra` and `Hecke`, the `two_block_group_algebra_codes` will be available.")
    end
    return ext.two_block_group_algebra_codes(args...; kwargs...)
end

function generalized_bicycle_codes(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `generalized_bicycle_codes` depends on the package `AbstractAlgebra` and `Hecke` but you have not installed or imported them yet. Immediately after you import `AbstractAlgebra` and `Hecke`, the `generalized_bicycle_codes` will be available.")
    end
    return ext.generalized_bicycle_codes(args...; kwargs...)
end

function bicycle_codes(args...; kwargs...)
    ext = Base.get_extension(QuantumClifford, :QuantumCliffordHeckeExt)
    if isnothing(ext)
        throw("The `bicycle_codes` depends on the package `AbstractAlgebra` and `Hecke` but you have not installed or imported them yet. Immediately after you import `AbstractAlgebra` and `Hecke`, the `bicycle_codes` will be available.")
    end
    return ext.bicycle_codes(args...; kwargs...)
end
