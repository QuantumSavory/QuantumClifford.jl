"""
$TYPEDSIGNATURES

Compute the distance of a code using mixed integer programming.
See [`QuantumClifford.ECC.DistanceMIPAlgorithm`](@ref) for configuration options.

Computes the minimum Hamming weight of a binary vector `x` by solving an **mixed
integer program (MIP)** that satisfies the following constraints:

- ``\\text{stab} \\cdot x \\equiv 0 \\pmod{2}``: The binary vector `x` must have an
even overlap with each `X`-check of the stabilizer binary representation `stab`.
- ``\\text{logicOp} \\cdot x \\equiv 1 \\pmod{2}``: The binary vector `x` must have
an odd overlap with logical-`X` operator `logicOp` on the `i`-th logical qubit.

Specifically, it calculates the minimum Hamming weight ``d_{Z}`` for the `Z`-type
logical operator. The minimum distance for `X`-type logical operators is the same.

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
consider that the `[[7, 1, 3]]` Steane code has a minimum distance of three, despite all
its elements having a weight of four. This discrepancy occurs because stabilizer codes
are defined by parity-check matrices, while their minimum distances are determined by
the dual [Sabo:2022smk](@cite).

```jldoctest
julia> using QuantumClifford, QuantumClifford.ECC

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

!!! note
    The minimum distance problem for quantum codes is *NP-hard*, and this hardness extends
    to multiplicative and additive approximations, even when restricted to stabilizer or
    CSS codes, with the result established through a reduction from classical problems in
    the CWS framework using a 4-cycle free graph [kapshikar2023hardness](@cite). Despite
    this, methods that improve on brute-force approaches are actively explored.

For a more in-depth background on minimum distance, see Chapter 3 of [Sabo:2022smk](@cite).

### Mixed Integer Programming

The MIP minimizes the Hamming weight `w(x)`, defined as the number of nonzero
elements in `x`, while satisfying the constraints:

```math
\\begin{aligned}
    \\text{Minimize} \\quad & w(x) = \\sum_{i=1}^n x_i, \\\\
    \\text{subject to} \\quad & \\text{stab} \\cdot x \\equiv 0 \\pmod{2}, \\\\
                            & \\text{logicOp} \\cdot x \\equiv 1 \\pmod{2}, \\\\
                            & x_i \\in \\{0, 1\\} \\quad \\text{for all } i.
\\end{aligned}
```

Here:
- ``\\text{stab}`` is the binary matrix representing the stabilizer group.
- ``\\text{logicOp}`` is the binary vector representing the logical-`X` operator.
- `x` is the binary vector (decision variable) being optimized.

The optimal solution ``w_{i}`` for each logical-`X` operator corresponds to the
minimum weight of a Pauli `Z`-type operator satisfying the above conditions.
The `Z`-type distance is given by:

```math
\\begin{aligned}
    d_Z = \\min(w_1, w_2, \\dots, w_k),
\\end{aligned}
```

where `k` is the number of logical qubits.

# Example

A [[40, 8, 5]] 2BGA code with the minimum distance of 5 from
Table 2 of [lin2024quantum](@cite).

```jldoctest jumpexamples
julia> import Hecke: group_algebra, GF, abelian_group, gens; import HiGHS; import JuMP;

julia> using QuantumClifford, QuantumClifford.ECC

julia> l = 10; m = 2;

julia> GA = group_algebra(GF(2), abelian_group([l,m]));

julia> x, s = gens(GA);

julia> A = 1 + x^6;

julia> B = 1 + x^5 + s + x^6 + x + s*x^2;

julia> c = two_block_group_algebra_codes(A,B);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(40, 8, 5)
```

A [[48, 6, 8]] GB code with the minimum distance of 8 from (A3)
in Appendix B of [panteleev2021degenerate](@cite).

```jldoctest jumpexamples
julia> l = 24;

julia> c1 = generalized_bicycle_codes([0, 2, 8, 15], [0, 2, 12, 17], l);

julia> code_n(c1), code_k(c1), distance(c1, DistanceMIPAlgorithm(solver=HiGHS))
(48, 6, 8)
```

### Applications

Mixed-integer programming (MIP) is applied in quantum error correction,
notably for decoding and minimum distance computation. Some applications
are as follows:

- The first usecase of the MIP approach was the code capacity Most Likely Error (MLE) decoder for color codes introduced in [landahl2011color](@cite).

- For all quantum LDPC codes presented in [panteleev2021degenerate](@cite), the lower and upper bounds on the minimum distance was obtained by reduction to a mixed integer linear program and using the GNU Linear Programming Kit ([makhorin2008glpk](@cite)).

- For all the Bivariate Bicycle (BB) codes presented in [bravyi2024high](@cite), the code distance was calculated using the mixed integer programming approach.

- [lacroix2024scaling](@cite) developed a MLE decoder that finds the most likely chain of Pauli errors given the observed error syndrome by solving a mixed-integer program using `HiGHS` package ([huangfu2018parallelizing](@cite)).

- [cain2025correlateddecodinglogicalalgorithms](@cite) formulate maximum-likelihood decoding as a mixed-integer program maximizing ``\\prod_{j=1}^M p_j^{E_j}(1-p_j)^{1-E_j}`` (where binary variables ``E_j \\in {0,1}`` indicate error occurrence) subject to syndrome constraints, solved optimally via MIP solvers despite its NP-hard complexity.
"""
function distance(code::AbstractECC, alg::DistanceMIPAlgorithm)
    logical_qubits = isnothing(alg.logical_qubit) ? (1:code_k(code)) : (alg.logical_qubit:alg.logical_qubit)
    isnothing(alg.logical_qubit) || (1 <= alg.logical_qubit <= code_k(code)) || throw(ArgumentError(THROW_BOUNDS))
    # Get the appropriate logical operators and matrices based on operator type
    logical_operator_type = alg.logical_operator_type
    l, H, h = if logical_operator_type == :X
        l_val = SparseMatrixCSC{Int, Int}(stab_to_gf2(logx_ops(code)))
        H_val = SparseMatrixCSC{Int, Int}(stab_to_gf2(parity_checks(code)))
        px = parity_matrix_x(code)
        h_val = cat(px, spzeros(Int, size(px, 1), nqubits(code)); dims=2)
        (l_val, H_val, h_val)
    else
        l_val = SparseMatrixCSC{Int, Int}(stab_to_gf2(logz_ops(code)))
        H_val = SparseMatrixCSC{Int, Int}(stab_to_gf2(parity_checks(code)))
        pz = parity_matrix_z(code)
        h_val = cat(spzeros(Int, size(pz, 1), nqubits(code)), pz; dims=2)
        (l_val, H_val, h_val)
    end
    # Calculate weights for each logical qubit
    if alg.opt_summary
        weights = Dict(
            i => (weight=w, summary=s)
            for (i, (w, s)) in (
                (i, _minimum_distance(h, l[i, :], alg.solver.Optimizer, true, alg.time_limit))
                for i in logical_qubits
            )
        )
    else
        weights = Dict(
            i => first(_minimum_distance(h, l[i, :], alg.solver.Optimizer, false, alg.time_limit))
            for i in logical_qubits
        )
    end
    # Return appropriate result based on configuration
    if alg.opt_summary
        min_entry = argmin(x -> x.weight, values(weights))
        return (min_entry.weight, min_entry.summary)
    else
        return minimum(values(weights))
    end
end

"""Computing minimum distance for quantum LDPC codes is a NP-Hard problem,
but we can solve mixed integer program (MIP) for small instance sizes.
inspired from [bravyi2024high](@cite) and [landahl2011color](@cite)."""
function _minimum_distance(hx, lx, opt, opt_summary, time_limit)
    n = size(hx, 2) # number of qubits
    m = size(hx, 1) # number of stabilizers
    # Maximum stabilizer weight
    whx = maximum(sum(hx[i, :] for i in 1:m)) # Maximum row sum of hx
    # Weight of the logical operator
    wlx = count(!iszero, lx) # Count non-zero elements in lx
    # Number of slack variables needed to express orthogonality constraints modulo two
    num_anc_hx = ilog2(whx, RoundUp)
    num_anc_logical = ilog2(wlx, RoundUp)
    num_var = n + m * num_anc_hx + num_anc_logical
    # Let the user choose which MIP solver to use
    model = Model(opt; add_bridges = false) # model
    set_silent(model)
    set_time_limit_sec(model, time_limit)
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
    # Ensure the model is solved and feasible
    if !is_solved_and_feasible(model)
        if termination_status(model) == MOI.MEMORY_LIMIT
            throw(ErrorException(THROW_MODEL_MEMORY_LIMIT))
        elseif termination_status(model) == MOI.TIME_LIMIT
            throw(ErrorException(THROW_MODEL_TIME_LIMIT))
        else
            throw(ErrorException("Model failed to solve :$(termination_status(model))"))
        end
    end
    opt_val = sum(value(x[i]) for i in 1:n)
    opt_summary ? (round(Int, opt_val), solution_summary(model)) : round(Int, opt_val)
end
