"""
Construct a regular LDPC code parity-check matrix `H` using Gallager's
original construction method [1057683](@cite).

# Mathematical Formulation

The parity-check matrix `H` is constructed as:

```math
\\begin{aligned}
H = \\begin{bmatrix} H_1 \\\\ H_2 \\\\ \\vdots \\\\ H_{w_c} \\end{bmatrix}
\\end{aligned}
```

where each submatrix ``H_d`` is a ``μ × μw_r`` binary matrix with:
- Row weight ``w_r`` (number of 1s per row)
- Column weight 1 (each column has exactly one 1)

The first submatrix ``H_1`` has a specific structure where row ``i`` (``1 ≤ i ≤ μ``) contains 
its ``w_r`` 1s in columns ``(i-1)w_r + 1`` to ``iw_r``. Subsequent submatrices ``H_2,...,H_{w_c}`` 
are *column permutations* of ``H_1``.

### Fields
- ``μ`` : Integer determining submatrix size (controls code dimensions).
- ``w_c`` : Column weight of final H matrix (number of submatrices,` ≥ 3` recommended).
- ``w_r`` : Row weight of final H matrix (must satisfy ``w_r > w_c``).
- `rng` : Random number generator.
"""
function random_gallager_ldpc_code(μ::Int, w_c::Int, w_r::Int; rng::AbstractRNG)
    μ ≤ 0 && throw(ArgumentError("μ must be positive"))
    w_c < 3 && throw(ArgumentError("w_c must be ≥ 3 (Gallager's condition)"))
    w_r ≤ w_c && throw(ArgumentError("w_r must be > w_c"))
    n = μ * w_r
    m = μ * w_c
    I = Int[]
    J = Int[]
    V = Bool[]
    sizehint!(I, m * w_c)
    sizehint!(J, m * w_c)
    sizehint!(V, m * w_c)
    template_cols = [((i-1)*w_r+1):(i*w_r) for i in 1:μ]
    for d in 1:w_c
        rows = ((d-1)*μ+1):(d*μ)
        if d == 1
            for (i, cols) in enumerate(template_cols)
                for j in cols
                    push!(I, rows[i])
                    push!(J, j)
                    push!(V, true)
                end
            end
        else
            perm = randperm(rng, n)
            for (i, cols) in enumerate(template_cols)
                for j in cols
                    push!(I, rows[i])
                    push!(J, perm[j])
                    push!(V, true)
                end
            end
        end
    end
    return sparse(I, J, V, m, n)
end
