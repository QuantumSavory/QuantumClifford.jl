"""Convert a tableau to a memory layout that is fast for row operations.

In this layout a Pauli string (a row of the tableau) is stored contiguously in memory.
This corresponds to column-major storage of the underlying `xzs` matrix.

**This is the default layout** for [`Stabilizer`](@ref), [`Destabilizer`](@ref), 
[`MixedStabilizer`](@ref), and [`MixedDestabilizer`](@ref).

Optimal for:
- Row operations (e.g., `mul_left!`)
- Canonicalization (`canonicalize!`)
- Projective measurements (`project!`)
- Measuring large Pauli operators
- Applying large dense n-qubit Clifford operations

# Indexing Convention
Regardless of memory layout, indexing is always `xzs[tableau_column_index, tableau_row_index]`.
The X component of qubit `i` in stabilizer `j` is stored at `xzs[i_big, j]` and 
the Z component at `xzs[i_big + end÷2, j]`, where `i_big` accounts for bit packing.

See also: [`fastcolumn`](@ref)"""
function fastrow end

"""Convert a tableau to a memory layout that is fast for column operations.

In this layout a column of the tableau (the bits of a given qubit) is stored 
(mostly) contiguously in memory. This corresponds to row-major storage of the 
underlying `xzs` matrix, implemented via `transpose(collect(transpose(xzs)))`.

Due to bitpacking (e.g., packing 64 bits into a single `UInt64`),
the memory layout is not perfectly contiguous,
but it is still optimal given that some bitwrangling is required to extract a given bit.

**This is the default layout** for [`PauliFrame`](@ref).

Optimal for:
- Column operations (applying single-qubit and two-qubit gates)
- Sparse gate applications like `apply!(s, sCNOT(i, j))`
- Operations that only touch a few qubits

# Indexing Convention
Regardless of memory layout, indexing is always `xzs[tableau_column_index, tableau_row_index]`.
The X component of qubit `i` in stabilizer `j` is stored at `xzs[i_big, j]` and 
the Z component at `xzs[i_big + end÷2, j]`, where `i_big` accounts for bit packing.

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
