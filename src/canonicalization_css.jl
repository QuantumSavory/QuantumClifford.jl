"""
    lc_to_css(stabilizer::Stabilizer)

Find a local Clifford transformation to a Calderbank–Shor–Steane (CSS) code.
Return `(css, circuit)`, or `nothing` if no such transformation exists.
The input is not modified.

`css` is an independent generating set whose rows contain only X and I, or
only Z and I, including their signs. `circuit` is a vector of
[`SingleQubitOperator`](@ref)s: applying it to the input produces the same
signed stabilizer group as `css`, up to generator changes. There is at most
one gate per qubit. The output is a CSS representative, not a unique form
shared by all locally Clifford-equivalent inputs. If the input already
generates a CSS group, the circuit is empty.

Redundant generators, identity rows, unused qubits, and zero qubits are
supported. The generators must commute, be Hermitian, and generate a group
that excludes `-I`; otherwise an `ArgumentError` is thrown.

Uses the preserving-projector algorithm of Mingyuan Wang and Stefan Krastanov,
*Local Clifford equivalence to CSS codes: Polynomial-time classification and
conversion*. On each indecomposable qubit component, at most two affine
binary systems determine the local Clifford gates. For `m` input generators
on `n` qubits, the cost is O(m n² + n⁴) binary operations and O(m n + n³)
storage. Signed row reduction extracts the final pure checks.

```jldoctest
julia> s = S"XZ
             ZX";

julia> css, circuit = lc_to_css(s);

julia> css
+ XX
+ ZZ

julia> converted = copy(s);

julia> foreach(gate -> apply!(converted, gate), circuit);

julia> canonicalize!(converted) == css
true
```

See also: [`canonicalize!`](@ref), [`graphstate`](@ref).
"""
function lc_to_css(stabilizer::Stabilizer)
    all(iseven, phases(stabilizer)) || throw(ArgumentError("Stabilizer generators must be Hermitian."))
    reduced, _, r = canonicalize!(copy(stabilizer); ranks=true)
    all(iszero, @view phases(reduced)[r+1:end]) || throw(ArgumentError("Stabilizer generators must not generate -I."))
    basis = reduced[1:r]
    all(iseven, phases(basis)) && all(iszero(comm(basis, i, j)) for i in 1:r for j in 1:i-1) ||
        throw(ArgumentError("Stabilizer generators must commute."))
    circuit = SingleQubitOperator[]
    all(i -> all(q -> !basis[i,q][1], 1:nqubits(basis)) ||
             all(q -> !basis[i,q][2], 1:nqubits(basis)), 1:r) && return basis, circuit

    # canonicalize! gives RREF in grouped coordinates (x₁,…,xₙ,z₁,…,zₙ).
    # Its row supports identify the indecomposable physical factors, even
    # when the supplied generators mix otherwise independent code blocks.
    n = nqubits(basis)
    G = stab_to_gf2(basis)
    pivots = [findfirst(@view G[i,:])::Int for i in 1:r]
    support = Graphs.SimpleGraph(n)
    for i in 1:r, q in 1:n
        if G[i,q] || G[i,n+q]
            Graphs.add_edge!(support, mod1(pivots[i], n), q)
        end
    end
    for qubits in Graphs.connected_components(support)
        rows = findall(p -> mod1(p, n) in qubits, pivots)
        isempty(rows) && continue
        component = G[rows, vcat(qubits, n .+ qubits)]
        system = _css_preservation_system(component)
        # Pᵢ = [aᵢ bᵢ; cᵢ 1+aᵢ]. Preservation makes bᵢcᵢ constant
        # on a component, so b₁=0 or c₁=0 enforces Pᵢ²=Pᵢ everywhere.
        solution = _css_solve_affine(system, 2)
        if isnothing(solution)
            solution = _css_solve_affine(system, 3)
        end
        isnothing(solution) && return nothing
        for (i, q) in enumerate(qubits)
            a, b, c = solution[3i-2:3i]
            gate = _css_local_clifford(q, a, b, c)
            gate == SingleQubitOperator(sId1(q)) || push!(circuit, gate)
        end
    end
    # A CSS row space has pure rows in RREF. Using the tableau operations
    # tracks signs of products such as (+YY)(+ZZ) = -XX automatically.
    for gate in circuit
        apply!(basis, gate)
    end
    canonicalize!(basis), circuit
end

# Build G P Hᵀ = 0, where H spans the ordinary dot-product annihilator of
# the RREF basis G. Each row of H corresponds to one nonpivot column.
# Columns of the augmented system are (a₁,b₁,c₁,…,aₙ,bₙ,cₙ | rhs).
function _css_preservation_system(G)
    r, columns = size(G)
    n = columns ÷ 2
    pivots = [findfirst(@view G[i,:])::Int for i in 1:r]
    nonpivots = setdiff(1:columns, pivots)
    H = falses(length(nonpivots), columns)
    for (j, column) in enumerate(nonpivots)
        H[j,column] = true
        H[j,pivots] .= @view G[:,column]
    end
    # Reserve the final row for either anchor constraint.
    system = falses(r * length(nonpivots) + 1, 3n + 1)
    row = 0
    for α in 1:r, β in eachindex(nonpivots)
        row += 1
        for i in 1:n
            x, z = G[α,i], G[α,n+i]
            u, v = H[β,i], H[β,n+i]
            system[row,3i-2] = (x & u) ⊻ (z & v)
            system[row,3i-1] = x & v
            system[row,3i] = z & u
            system[row,end] ⊻= z & v
        end
    end
    system
end

# Set free variables to zero. A pivot in the augmented column certifies
# inconsistency; otherwise each pivot variable equals its row's RHS.
function _css_solve_affine(system, anchor)
    augmented = copy(system)
    augmented[end,anchor] = true
    gf2_gausselim!(augmented)
    variables = size(augmented, 2) - 1
    solution = falses(variables)
    for row in eachrow(augmented)
        pivot = findfirst(row)
        isnothing(pivot) && continue
        pivot > variables && return nothing
        solution[pivot] = row[end]
    end
    solution
end

function _css_local_clifford(q, a, b, c)
    directions = ((true, false), (false, true), (true, true))
    image = only(v for v in directions if ((a & v[1]) ⊻ (b & v[2]), (c & v[1]) ⊻ (!a & v[2])) == v)
    kernel = only(v for v in directions if !((a & v[1]) ⊻ (b & v[2])) && !((c & v[1]) ⊻ (!a & v[2])))
    # Q has the column image/kernel vectors: P Q = Q diag(1,0).
    # Row labels transform as (x,z) Q. Positive images of X and Z fix
    # a Clifford lift; apply! handles its induced sign on Y.
    SingleQubitOperator(q, image[1], kernel[1], image[2], kernel[2], false, false)
end
