"""
    $TYPEDEF

Construct a regular LDPC code parity-check matrix `H` using Gallager's
original construction method [gallager1962ldpc](@cite).

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

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/gallager).

# Example

Here is an example of `(3,4)`-LDPC code

```jldoctest examples
julia> using QECCore

julia> μ = 3; wc = 3; wr = 4;

julia> H = random_Gallager_ldpc(μ, wc, wr);
```

The Gallager's classical code exhibits regular structure with all columns
having weight 3 and all rows having weight 4, ensuring it maintains constant
column and row weights throughout the parity-check matrix.

```jldoctest examples
julia> all(sum(Matrix(H), dims=1) .== 3)
true

julia> all(sum(Matrix(H), dims=2) .== 4)
true
```

### Arguments
- `block_rows::Int` Number of block rows called "submatrices" in the Gallager's construction [gallager1962ldpc](@cite).
- `col_weight::Int`: Column weight
- `row_weight::Int`: Row weight which must be greater column weight.
"""
function random_Gallager_ldpc(rng::AbstractRNG, block_rows::Int, col_weight::Int, row_weight::Int)
    block_rows ≤ 1 && throw(ArgumentError("Number of block rows must be > 1, got $block_rows"))
    col_weight < 3 && throw(ArgumentError("Column weight must be ≥ 3 (Gallager's condition)"))
    row_weight ≤ 1 && throw(ArgumentError("Row weight must be > 1, got $row_weight"))
    row_weight ≤ col_weight && throw(ArgumentError("Row weight must be > column weight"))
    n = block_rows*row_weight
    m = block_rows*col_weight
    I = Int[]
    J = Int[]
    V = Bool[]
    sizehint!(I, m*col_weight)
    sizehint!(J, m*col_weight)
    sizehint!(V, m*col_weight)
    template_cols = [((i-1)*row_weight+1):(i*row_weight) for i in 1:block_rows]
    perms = [randperm(rng, n) for _ in 2:col_weight]
    for d in 1:col_weight
        rows = ((d-1)*block_rows+1):(d*block_rows)
        if d == 1
            for (i, cols) in enumerate(template_cols)
                for j in cols
                    push!(I, rows[i])
                    push!(J, j)
                    push!(V, true)
                end
            end
        else
            perm = perms[d-1]
            for (i, cols) in enumerate(template_cols)
                for j in cols
                    push!(I, rows[i])
                    push!(J, perm[j])
                    push!(V, true)
                end
            end
        end
    end
    H = sparse(I, J, V, m, n)
    return H
end

random_Gallager_ldpc(block_rows::Int, col_weight::Int, row_weight::Int) = random_Gallager_ldpc(GLOBAL_RNG, block_rows, col_weight, row_weight)
