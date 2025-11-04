
#=============================================================================#
# TODO: Remove these once the main package establishes them.
@inline function copy_to!(
    target::PauliOperator, source::PauliOperator
    )

    target.nqubits == source.nqubits || throw(ArgumentError("BAD COPY_TO!"))
    copyto!(target.phase, source.phase)
    copyto!(target.xz, source.xz)
    return target

end

@inline function copy_to!(
    target::Tableau, source::Tableau
    )

    target.nqubits == source.nqubits || throw(ArgumentError("BAD COPY_TO!"))
    copyto!(target.phases, source.phases)
    copyto!(target.xzs, source.xzs)
    return target

end

@inline function copy_to!(
    target::AbstractStabilizer, source::AbstractStabilizer
    )

    copy_to!(tab(target), tab(source))
    return target

end
#=============================================================================#
