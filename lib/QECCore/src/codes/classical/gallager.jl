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

```jldoctest
julia> using QECCore; using Nemo

julia> μ = 3; wc = 3; wr = 4;

julia> c = GallagerLDPC(μ, wc, wr);

julia> H = parity_matrix(c)
9×12 SparseArrays.SparseMatrixCSC{Bool, Int64} with 36 stored entries:
 1  1  1  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  1  1  1  1  ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  1  1  1
 1  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  1  ⋅  ⋅  ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  1  ⋅  ⋅  1  1
 ⋅  1  ⋅  1  1  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅
 ⋅  1  1  ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1
 ⋅  ⋅  ⋅  1  ⋅  1  1  1  ⋅  ⋅  ⋅  ⋅
 1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  1  1  ⋅

julia> code_n(c), code_k(c)
(12, 5)

julia> all(sum(Matrix(H), dims=1) .== 3)
true

julia> all(sum(Matrix(H), dims=2) .== 4)
true
```

### Fields
    $TYPEDFIELDS
"""
struct GallagerLDPC <: AbstractCECC # TODO this should not be a constructor but a function, given that it just randomly generates a matrix
    """Number of block rows called "submatrices" in the Gallager's construction [gallager1962ldpc](@cite)."""
    block_rows::Int
    """Column weight"""
    col_weight::Int
    """Row weight which must be greater column weight."""
    row_weight::Int
    """Random seed for reproducible parity check matrix generation (`default=42`)."""
    seed::Int

    function GallagerLDPC(block_rows, col_weight, row_weight, seed::Union{Int,Nothing}=42)
        block_rows ≤ 1 && throw(ArgumentError("Number of block rows must be > 1, got $block_rows"))
        col_weight < 3 && throw(ArgumentError("Column weight must be ≥ 3 (Gallager's condition)"))
        row_weight ≤ 1 && throw(ArgumentError("Row weight must be > 1, got $row_weight"))
        row_weight ≤ col_weight && throw(ArgumentError("Row weight must be > column weight"))
        if isnothing(seed)
            seed = rand(UInt32)
        end
        new(block_rows, col_weight, row_weight, seed)
    end
end

function _gallager_ldpc_code(μ::Int, w_c::Int, w_r::Int, seed::Int)
    n = μ*w_r
    m = μ*w_c
    I = Int[]
    J = Int[]
    V = Bool[]
    sizehint!(I, m*w_c)
    sizehint!(J, m*w_c)
    sizehint!(V, m*w_c)
    template_cols = [((i-1)*w_r+1):(i*w_r) for i in 1:μ]
    rng = MersenneTwister(seed)
    perms = [randperm(rng, n) for _ in 2:w_c]
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
    return sparse(I, J, V, m, n)
end

code_n(c::GallagerLDPC) = c.block_rows*c.row_weight

parity_matrix(c::GallagerLDPC) = _gallager_ldpc_code(c.block_rows, c.col_weight, c.row_weight, c.seed)
