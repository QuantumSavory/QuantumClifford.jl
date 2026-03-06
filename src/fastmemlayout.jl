"""Convert a tableau to a memory layout that is fast for row operations.

In this layout a Pauli string (a row of the tableau) is stored contiguously in memory.

See also: [`fastcolumn`](@ref)"""
function fastrow end

"""Convert a tableau to a memory layout that is fast for column operations.

In this layout a column of the tableau is stored (mostly) contiguously in memory.
Due to bitpacking, e.g., packing 64 bits into a single `UInt64`,
the memory layout is not perfectly contiguous,
but it is still optimal given that some bitwrangling is required to extract a given bit.

See also: [`fastrow`](@ref)"""
function fastcolumn end

fastrow(t::Tableau{Tₚᵥ,Tₘ}) where {Tₚᵥ, Tₘ} = t
fastrow(t::Tableau{Tₚᵥ,Tₘ}) where {Tₚᵥ, Tₘ<:Adjoint} = Tableau(t.phases, t.nqubits, collect(t.xzs))
fastcolumn(t::Tableau{Tₚᵥ,Tₘ}) where {Tₚᵥ, Tₘ} = Tableau(t.phases, t.nqubits, collect(t.xzs')')
fastcolumn(t::Tableau{Tₚᵥ,Tₘ}) where {Tₚᵥ, Tₘ<:Adjoint} = t

fastrow(s::Stabilizer) = Stabilizer(fastrow(tab(s)))
fastcolumn(s::Stabilizer) = Stabilizer(fastcolumn(tab(s)))

fastrow(s::Destabilizer) = Destabilizer(fastrow(tab(s)))
fastcolumn(s::Destabilizer) = Destabilizer(fastcolumn(tab(s)))

fastrow(s::MixedStabilizer) = MixedStabilizer(fastrow(tab(s)), rank(s))
fastcolumn(s::MixedStabilizer) = MixedStabilizer(fastcolumn(tab(s)), rank(s))

fastrow(s::MixedDestabilizer) = MixedDestabilizer(fastrow(tab(s)), rank(s))
fastcolumn(s::MixedDestabilizer) = MixedDestabilizer(fastcolumn(tab(s)), rank(s))
