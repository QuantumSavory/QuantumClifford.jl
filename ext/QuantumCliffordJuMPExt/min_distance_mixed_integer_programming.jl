"""
TYPEDSIGNATURES

Computes the minimum Hamming weight of a binary vector `x` by solving an **mixed
integer program (MIP)** that satisfies the following constraints:

- \$\\text{stab} \\cdot x \\equiv 0 \\pmod{2}\$: The binary vector x must have an
even overlap with each X-check of the stabilizer binary representation `stab`.
- \$\\text{logicOp} \\cdot x \\equiv 1 \\pmod{2}\$: The binary vector x must have
an odd overlap with logical-X operator `logicOp` on the i-th logical qubit.

Specifically, it calculates the minimum Hamming weight \$d_{Z}\$ for the Z-type
logical operator. The minimum distance for X-type logical operators is the same.

### Background on Minimum Distance

For *classical* codes, the minimum distance, which measures a code's error-correcting
capability, is equivalent to its minimum weight. This can be computed by generating
all possible codewords from combinations of the generator matrix rows, calculating
their weights, and finding the smallest. While accurate, this method takes exponential
time. Vardy [vardy1997intractability](@cite) demonstrated that computing the minimum
distance is *NP-hard*, and the corresponding decision problem is *NP-complete*,
making polynomial-time algorithms unlikely.

For *quantum* codes, classical intuition does not always apply. The minimum distance
is given by the minimum weight of a non-trivial logical operator. This is generally
unrelated to the minimum distance of the corresponding stabilizer code when viewed as
a classical, additive code. White and Grassl [white2006new](@cite) proposed mapping
quantum codes to higher-dimensional classical linear codes. This mapping allows the
minimum distance of the quantum additive code to be inferred from that of the classical
linear code but increases parameters from `n` to `3n` and `d` to `2d`, adding complexity.
Furthermore, once a minimal weight vector is identified, it is essential to verify
whether it belongs to the Pauli group `𝒫ₙ` over `n` qubits [Sabo:2022smk](@cite).

Additionally, to illustrate how classical intuition can be misleading in this context,
consider that the [[7, 1, 3]] Steane code has a minimum distance of three, despite all
its elements having a weight of four. This discrepancy occurs because stabilizer codes
are defined by parity-check matrices, while their minimum distances are determined by
the dual [Sabo:2022smk](@cite).

```jldoctest example
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

julia> minimum(sum(stab_to_gf2(c), dims=2))
4

julia> distance(Steane7())
3
```

!!! note The minimum distance problem for quantum codes is *NP-hard*, and this hardness
extends to multiplicative and additive approximations, even when restricted to stabilizer
or CSS codes, with the result established through a reduction from classical problems in
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

!!!note Since the [[48, 6, 8]] GB code does not have specific lower and upper
bounds (e.g., consider [[48, 6, 5 ≤ d ≤ 8]]), the minimum distance for all
`Z`-type and `X`-type logical operators remains the same. Here, the *exact*
minimum distance of `8` is provided.

```jldoctest examples
julia> distance(c1, all_logical_qubits=true)
Dict{Int64, Int64} with 6 entries:
  5 => 8
  4 => 8
  6 => 8
  2 => 8
  3 => 8
  1 => 8
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
function distance(c::Stabilizer; upper_bound=false, logical_qubit=code_k(c), all_logical_qubits=false, logical_operator_type=:X)
    1 <= logical_qubit <= code_k(c) || throw(ArgumentError("The number of logical qubit must be between 1 and $(code_k(c)) inclusive"))
    logical_operator_type == :X || logical_operator_type == :Z || throw(ArgumentError("Invalid type of logical operator: Use :X or :Z"))
    l = get_lx_lz(c)[1]
    H = SparseMatrixCSC{Int, Int}(stab_to_gf2(parity_checks(c)))
    h = get_stab(H, :X)
    if logical_operator_type == :Z
        l = get_lx_lz(c)[2]
        H = SparseMatrixCSC{Int, Int}(stab_to_gf2(parity_checks(c)))
        h = get_stab(H, :Z)
    end
    weights = Dict{Int, Int}()
    for i in 1:logical_qubit
        w = _minimum_distance(h, l[i, :])
        weights[i] = w
    end
    return upper_bound ? maximum(values(weights)) : (all_logical_qubits ? weights : minimum(values(weights)))
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
