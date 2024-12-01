"""
TYPEDSIGNATURES

Computes the minimum Hamming weight of a binary vector `x` by solving an **mixed
integer program (MIP)** that satisfies the following constraints:

- \$\\text{stab} \\cdot x \\equiv 0 \\pmod{2}\$: The binary vector x must have an
even overlap with each X-check of the stabilizer binary representation `stab`.
- \$\\text{logicOp} \\cdot x \\equiv 1 \\pmod{2}\$: The binary vector x must have
an odd overlap with logical-X operator `logicOp` on the i-th logical qubit.

Specifically, it calculates the minimum Hamming weight \$d_{Z}\$ for the Z-type
logical operator, and the minimum distance for X-type logical operators is the same.

### Background on Minimum Distance

For *classical* codes, the minimum distance, which measures a code's error-correcting
capability, is equivalent to its minimum weight. This can be computed by generating
all possible codewords from combinations of the generator matrix rows, calculating
their weights, and finding the smallest. While accurate, this method takes exponential
time. Vardy [vardy1997intractability](@cite) demonstrated that computing the minimum
distance is *NP-hard*, and the corresponding decision problem is *NP-complete*,
making polynomial-time algorithms unlikely.

In the case of *quantum* codes, classical intuition does not always apply. For instance,
the Steane code has a minimum distance of three, even though all its elements have weight
four [Sabo:2022smk](@cite).

```jldoctest example1
julia> using QuantumClifford.ECC: Steane7, distance;

julia> c = parity_checks(Steane7());

julia> stab_to_gf2(c)
6×14 Matrix{Bool}:
 0  0  0  1  1  1  1  0  0  0  0  0  0  0
 0  1  1  0  0  1  1  0  0  0  0  0  0  0
 1  0  1  0  1  0  1  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  1  1  1  1
 0  0  0  0  0  0  0  0  1  1  0  0  1  1
 0  0  0  0  0  0  0  1  0  1  0  1  0  1

julia> distance(Steane7())
3
```

Such discrepancies arise because stabilizer codes are defined by parity-check matrices,
but their minimum distances are determined by the dual, specifically the minimum weight
of non-trivial logical operators.

```jldoctest example1
julia> lz = stab_to_gf2(logicalzview(canonicalize!(MixedDestabilizer(c))))
1×14 Matrix{Bool}:
 0  0  1  0  1  1  0  0  0  0  0  0  0  0

julia> lz = stab_to_gf2(logicalzview(canonicalize!(MixedDestabilizer(c))))
1×14 Matrix{Bool}:
 0  0  0  0  0  0  0  0  1  0  1  0  1  0
```

Brute-force methods remain viable but inefficient. White and Grassl [white2006new](@cite)
proposed mapping quantum codes to higher-dimensional classical linear codes. This mapping
allows the minimum distance of the quantum additive code to be inferred from that of the
classical linear code but increases parameters from `n` to `3n` and `d` to `2d`, adding
complexity.

The minimum distance problem for quantum codes is *NP-hard*, and this hardness extends
to multiplicative and additive approximations, even when restricted to stabilizer or
CSS codes, with the result established through a reduction from classical problems in
the CWS framework using a 4-cycle free graph [kapshikar2023hardness](@cite). Despite
this, methods that improve on brute-force approaches are actively explored.

For a more in-depth background on minimum distance, see Chapter 3 of [Sabo:2022smk](@cite).

### Mixed Integer Programming

The MIP minimizes the Hamming weight `w(x)`, defined as the number of nonzero
elements in `x`, while satisfying the constraints:

\$\$
\\begin{aligned}
    \\text{Minimize} \\quad & w(x) = \\sum_{i=1}^n x_i, \\\\
    \\text{subject to} \\quad & \\text{stab} \\cdot x \\equiv 0 \\pmod{2}, \\\\
                           & \\text{logicOp} \\cdot x \\equiv 1 \\pmod{2}, \\\\
                           & x_i \\in \\{0, 1\\} \\quad \\text{for all } i.
\\end{aligned}
\$\$

Here:
- \$\\text{stab}\$ is the binary matrix representing the stabilizer group.
- \$\\text{logicOp}\$ is the binary vector representing the logical-X operator.
- `x` is the binary vector (decision variable) being optimized.

The optimal solution \$w_{i}\$ for each logical-X operator corresponds to the
minimum weight of a Pauli Z-type operator satisfying the above conditions.
The Z-type distance is given by:

\$\$
\\begin{aligned}
    d_Z = \\min(w_1, w_2, \\dots, w_k),
\\end{aligned}
\$\$

where k is the number of logical qubits.

# Example

A [[40, 8, 5]] 2BGA code with the minimum distance of 5 from
Table 2 of [lin2024quantum](@cite).

```jldoctest examples
julia> import Hecke: group_algebra, GF, abelian_group, gens; import GLPK; import JuMP;

julia> using QuantumClifford.ECC: two_block_group_algebra_codes, generalized_bicycle_codes; # hide

julia> l = 10; m = 2;

julia> GA = group_algebra(GF(2), abelian_group([l,m]));

julia> x, s = gens(GA);

julia> A = 1 + x^6;

julia> B = 1 + x^5 + s + x^6 + x + s*x^2;

julia> c = two_block_group_algebra_codes(A,B);

julia> code_n(c), code_k(c), distance(c)
(40, 8, 5)
```

A [[48, 6, 8]] GB code with the minimum distance of 8 from (A3)
in Appendix B of [panteleev2021degenerate](@cite).

```jldoctest examples
julia> l = 24;

julia> c1 = generalized_bicycle_codes([0, 2, 8, 15], [0, 2, 12, 17], l);

julia> code_n(c1), code_k(c1), distance(c1)
(48, 6, 8)
```

### Applications

- The first usecase of the MIP approach was the code capacity Most Likely
Error (MLE) decoder for color codes introduced in [landahl2011fault](@cite).
- For all quantum LDPC codes presented in [panteleev2021degenerate](@cite),
the lower and upper bounds on the minimum distance was obtained by reduction
to a mixed integer linear program and using the GNU Linear Programming Kit.
- For all the Bivariate Bicycle (BB) codes presented in [bravyi2024high](@cite),
the code distance was calculated using the mixed integer programming approach.

"""
function distance(c::Stabilizer; num_logical_qubits=code_k(c), upper_bound=false)
    lx, _ = get_lx_lz(c)
    H = stab_to_gf2(parity_checks(c))
    H = SparseMatrixCSC{Int, Int}(H)
    hx = get_stab_hx(H)
    1 <= num_logical_qubits < code_k(c) && return _minimum_distance(hx, lx[num_logical_qubits, :]) # for large instances
    1 <= num_logical_qubits <= code_k(c) || throw(ArgumentError("The number of logical qubits must be between 1 and $(code_k(c)) inclusive"))
    if num_logical_qubits == code_k(c)
        weights = []
        for i in 1:num_logical_qubits
            w = _minimum_distance(hx, lx[num_logical_qubits, :])
            push!(weights, w)
        end
        if upper_bound
            return minimum(weights), maximum(weights) # return upper bound of minimum distance if required.
        else
            return minimum(weights)
        end
    end
end

# Computing minimum distance for quantum LDPC codes is a NP-Hard problem,
# but we can solve mixed integer program (MIP) for small instance sizes.
# inspired from [bravyi2024high](@cite) and [landahl2011fault](@cite).
function _minimum_distance(hx, lx)
    n = size(hx, 2) # number of qubits
    m = size(hx, 1) # number of stabilizers
    # Maximum stabilizer weight
    whx = maximum(sum(hx[i, :] for i in 1:m)) # Maximum row sum of hx
    # Weight of the logical operator
    wlx = count(!iszero, lx) # Count non-zero elements in lx
    # Number of slack variables needed to express orthogonality constraints modulo two
    num_anc_hx = ceil(Int, log2(whx))
    num_anc_logical = ceil(Int, log2(wlx))
    num_var = n + m * num_anc_hx + num_anc_logical
    model = Model(GLPK.Optimizer; add_bridges = false) # model
    set_silent(model)
    @variable(model, x[1:num_var], Bin) # binary variables
    @objective(model, Min, sum(x[i] for i in 1:n)) # objective function
    # Caching weights for orthogonality constraints
    orthogonality_weights = sparsevec(1:m, [Vector{Int}() for _ in 1:m], m)
    # Orthogonality to rows of hx constraints
    for row in 1:m
        weight = spzeros(Int, num_var)
        supp = findnz(hx[row, :])[1] # Non-zero indices in hx[row, :]
             for q in supp
                 weight[q] = 1
             end
             cnt = 1
             for q in 1:num_anc_hx
                 weight[n + (row - 1) * num_anc_hx + q] = -(1 << cnt)
             cnt += 1
         end
         orthogonality_weights[row] = weight # Cache the weight vector
    end
    # Add orthogonality constraints using cached weights
    for row in 1:m
        @constraint(model, sum(orthogonality_weights[row][i] * x[i] for i in 1:num_var) == 0)
    end
    # Odd overlap with lx constraint
    supp = findnz(lx)[1] # Non-zero indices in lx
    weight = spzeros(Int, num_var)
    for q in supp
        weight[q] = 1
    end
    cnt = 1
    for q in 1:num_anc_logical
        weight[n + m * num_anc_hx + q] = -(1 << cnt)
        cnt += 1
    end
    @constraint(model, sum(weight[i] * x[i] for i in 1:num_var) == 1)
    optimize!(model)
    opt_val = sum(value(x[i]) for i in 1:n)
    return Int(opt_val)
end
