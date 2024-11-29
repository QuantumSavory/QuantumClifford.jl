"""
TYPEDSIGNATURES

Computes the minimum Hamming weight of a binary vector `x` by solving an **mixed
integer program (MIP)** that satisfies the following constraints:

- \$\\text{stab} \\cdot x \\equiv 0 \\pmod{2}\$: The binary vector x must have an
even overlap with each X-check of the stabilizer binary representation `stab`.
- \$\\text{logicOp} \\cdot x \\equiv 1 \\pmod{2}\$: The binary vector x must have
an odd overlap with logical-X operator `logicOp` on the i-th logical qubit.

It computes the minimum Hamming weight \$d_{Z}\$ for the Z-type logical
operator. The distance for X-type logical operators is the same.

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

A [[40, 8, 5]] 2BGA code with the minimum distance of 5.

```jldoctest
julia> import Hecke: group_algebra, GF, abelian_group, gens; import GLPK; import JuMP;

julia> using QuantumClifford.ECC: minimum_distance, two_block_group_algebra_codes;

julia> l = 10; m = 2;

julia> GA = group_algebra(GF(2), abelian_group([l,m]));

julia> x, s = gens(GA);

julia> A = 1 + x^6;

julia> B = 1 + x^5 + s + x^6 + x + s*x^2;

julia> c = two_block_group_algebra_codes(A,B);

julia> minimum_distance(c)
5
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
function minimum_distance(c::Stabilizer; num_logical_qubits=code_k(c))
    lx, _ = get_lx_lz(c)
    H = stab_to_gf2(parity_checks(c))
    H = sparse(H)
    hx =  Matrix{Int}(get_stab_hx(H))
    num_logical_qubits == 1 && return _minimum_distance(hx, lx[num_logical_qubits, :])
    if 2 <= num_logical_qubits && num_logical_qubits <= code_k(c)
        w_values = []
        for i in 1:num_logical_qubits
            w = _minimum_distance(hx, lx[num_logical_qubits, :])
            push!(w_values, w)
        end
        return round(Int, sum(w_values)/num_logical_qubits)
    else
        throw(ArgumentError("The number of logical qubits must be between 1 to $(code_k(c)) inclusive"))
    end
end

# Computing minimum distance for quantum LDPC codes is a NP-Hard problem,
# but we can solve mixed integer program (MIP) for small instance sizes.
# inspired from [bravyi2024high](@cite) and [landahl2011fault](@cite).
function _minimum_distance(hx, lx)
    n = size(hx, 2) # number of qubits
    m = size(hx, 1) # number of stabilizers
    # Maximum stabilizer weight
    whx = maximum(sum(hx[i, :] for i in 1:m))  # Maximum row sum of hx
    # Weight of the logical operator
    wlx = count(!iszero, lx)  # Count non-zero elements in lx
    # Number of slack variables needed to express orthogonality constraints modulo two
    num_anc_hx = ceil(Int, log2(whx))
    num_anc_logical = ceil(Int, log2(wlx))
    num_var = n + m * num_anc_hx + num_anc_logical
    model = Model(GLPK.Optimizer) # model
    set_silent(model)
    @variable(model, x[1:num_var], Bin) # binary variables
    @objective(model, Min, sum(x[i] for i in 1:n)) # objective function
    # Orthogonality to rows of hx constraints
    for row in 1:m
        weight = zeros(Int, num_var)
        supp = findall(hx[row, :] .!= 0)  # Non-zero indices in hx[row, :]
             for q in supp
                 weight[q] = 1
             end
             cnt = 1
             for q in 1:num_anc_hx
                 weight[n + (row - 1) * num_anc_hx + q] = -(1 << cnt)
             cnt += 1
         end
        @constraint(model, sum(weight[i] * x[i] for i in 1:num_var) == 0)
    end
    # Odd overlap with lx constraint
    supp = findall(lx .!= 0)  # Non-zero indices in lx
    weight = zeros(Int, num_var)
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
