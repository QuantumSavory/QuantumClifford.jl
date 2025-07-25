"""
$TYPEDSIGNATURES

Compute the distance of a code using mixed integer programming.
See [`QuantumClifford.ECC.DistanceMIPAlgorithm`](@ref) for configuration options.

Specifically, it computes the minimum distance of a Calderbank-Shor-Steane code by
solving two independent **Mixed Integer Programs** for `X`-type (``d_X``) and `Z`-type
(``d_Z``) distances. The code distance is ``d = \\min(d_X, d_Z)``.

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

### Quantum Error Correction and Code Distance

The foundation of quantum error correction lies in protecting logical quantum information
from physical errors by encoding it across multiple qubits. A quantum code's performance is
fundamentally characterized by its distance (`d`), which quantifies the code's ability to detect and correct errors. Practically, the distance represents the minimum number of physical
qubit errors required to cause an undetectable logical error - one that corrupts encoded
information while evading the code's error detection mechanisms.

### Fundamentals of Quantum Code Distance

The distance `d` of a quantum error-correcting code represents its robustness against
physical errors and is defined as:

```math
\\begin{align*}
\\begin{equation}
d = \\min_{P \\in N(S)\\setminus S} \\mathrm{wt}(P)
\\end{equation}
\\end{align*}
```

where:

- ``N(S)`` denotes the normalizer of the stabilizer group `S`
- ``\\mathrm{wt}(P)`` represents the weight of Pauli operator `P` (number of non-identity components)
- The minimization is taken over all logical operators that are not stabilizers.

The distance reveals the essential property that it equals the smallest number of qubits that must be affected to produce a logical error undetectable by stabilizer measurements. The normalizer condition ``P \\in N(S)`` ensures the operator commutes with all stabilizers, while ``P \\notin S`` guarantees it performs a non-trivial logical operation.

### Mixed Integer Programming

We compute the minimum code distance for CSS (Calderbank-Shor-Steane) codes by solving MIPs. 

The distance is computed separately for `X`-type (``d_X``) and `Z`-type (``d_Z``) logical
operators, then combined to give the true code distance: ``d = \\min(d_X, d_Z)``.

#### `X-type` Distance (``d_X``)

It is defined as the minimum number of `Z`-errors (phase flips) required to implement a non-trivial `X`-logical operator, where the errors must both commute with all `X`-stabilizers and anti-commute with at least one `X`-logical operator.


```math
\\begin{align*}
\\text{Minimize} \\quad & \\sum_{i=1}^n e_{Z,i} \\quad \\text{(Hamming weight of Z-errors)} \\\\
\\text{Subject to} \\quad & \\mathbf{H_X} \\cdot \\mathbf{e}_Z \\equiv \\mathbf{0} \\pmod{2} \\quad \\text{(Commutes with X-stabilizers)} \\\\
                      & \\mathbf{L_X} \\cdot \\mathbf{e}_Z \\equiv 1 \\pmod{2} \\quad \\text{(Anti-commutes with X-logical)} \\\\
                      & e_{Z,i} \\in \\{0,1\\} \\quad \text{(Binary error variables)}
\\end{align*}
```

Here:
- ``\\mathbf{H_X}`` represent `X`-stabilizer matrix)
- ``\\mathbf{L_X}`` represent `X`-logical operators).
- ``\\mathbf{e_Z}`` are binary vector representing `Z`-error locations.

#### `Z-type` Distance (``d_Z``)

It is defined as the minimum number of `X`-errors (bit flips) required to implement a non-trivial `Z`-logical operator, where the errors must both commute with all `Z`-stabilizers and anti-commute with at least one `Z`-logical operator.

```math
\\begin{align*}
\\text{Minimize} \\quad & \\sum_{i=1}^n e_{X,i} \\quad \\text{(Hamming weight of X-errors)} \\\\
\\text{Subject to} \\quad & \\mathbf{H_Z} \\cdot \\mathbf{e}_X \\equiv \\mathbf{0} \\pmod{2} \\quad \\text{(Commutes with Z-stabilizers)} \\\\
                      & \\mathbf{L_Z} \\cdot \\mathbf{e}_X \\equiv 1 \\pmod{2} \\quad \\text{(Anti-commutes with Z-logical)} \\\\
                      & e_{X,i} \\in \\{0,1\\} \\quad \\text{(Binary error variables)}
\\end{align*}
```

Here:
- ``\\mathbf{H_Z}`` represent `Z`-stabilizer matrix)
- ``\\mathbf{L_Z}`` represent `Z`-logical operators).
- ``\\mathbf{e_X}`` are binary vector representing `X`-error locations.

# Example

A `[[40, 8, 5]]` 2BGA code with the minimum distance of `5` from
Table `2` of [lin2024quantum](@cite).

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

A `[[48, 6, 8]]` GB code with the minimum distance of `8` from (A3)
in Appendix `B` of [panteleev2021degenerate](@cite).

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
    isnothing(alg.logical_qubit) || (1 <= alg.logical_qubit <= code_k(code)) || throw(ArgumentError("Logical qubit out of range"))
    if alg.logical_operator_type == :minXZ
        dX = _compute_distance(code, logical_qubits, :X, alg)
        dZ = _compute_distance(code, logical_qubits, :Z, alg)
        return min(dX, dZ)
    else
        return _compute_distance(code, logical_qubits, alg.logical_operator_type, alg)
    end
end

function _compute_distance(code, logical_qubits, logical_operator_type, alg)
    # Get the appropriate logical operators and matrices based on operator type
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
            throw(ErrorException("Model exceeded memory limits"))
        elseif termination_status(model) == MOI.TIME_LIMIT
            throw(ErrorException("Model exceeded time limit"))
        else
            throw(ErrorException("Model failed to solve: $(termination_status(model))"))
        end
    end
    opt_val = sum(value(x[i]) for i in 1:n)
    opt_summary ? (round(Int, opt_val), solution_summary(model)) : round(Int, opt_val)
end
