"""Convert a tableau to a memory layout that is fast for row operations.

In this layout a Pauli string (a row of the tableau) is stored contiguously in memory.

See also: [`fastrow`](@ref)"""
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

fastrow(s::Stabilizer) = Stabilizer(fastrow(s.tab))
fastcolumn(s::Stabilizer) = Stabilizer(fastcolumn(s.tab))

fastrow(s::Destabilizer) = Destabilizer(fastrow(s.tab))
fastcolumn(s::Destabilizer) = Destabilizer(fastcolumn(s.tab))

fastrow(s::MixedStabilizer) = MixedStabilizer(fastrow(s.tab), s.rank)
fastcolumn(s::MixedStabilizer) = MixedStabilizer(fastcolumn(s.tab), s.rank)

fastrow(s::MixedDestabilizer) = MixedDestabilizer(fastrow(s.tab), s.rank)
fastcolumn(s::MixedDestabilizer) = MixedDestabilizer(fastcolumn(s.tab), s.rank)
