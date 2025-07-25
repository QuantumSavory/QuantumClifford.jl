"""
For circulant *bivariate bicycle codes*, the parity-check matrices ``H_X`` and
``H_Z`` satisfy the rank equality ``\\text{rank}(H_X) = \\text{rank}(H_Z)``. This
follows from a symmetry in the code structure induced by a self-inverse permutation
matrix. We define ``C_\\ell`` as an `\\ell \\times \\ell`` permutation matrix such
that its `i`-th column has a `1` in row `j = -i \\mod \\ell``. Similarly, we define
``C_m`` for size ``m \\times m``, and let ``C = C_\\ell \\otimes C_m``. Using the
property that this matrix satisfies:

```math
\\begin{aligned}
C_\\ell S_\\ell C_\\ell = S_\\ell^T, \\quad C_m S_m C_m = S_m^T,
\\end{aligned}
```

it follows that the transposes of the circulant blocks are related by:

```math
\\begin{aligned}
A^T = C A C, \\quad B^T = C B C.
\\end{aligned}
```

Hence, the matrix ``H_Z`` can be expressed in terms of ``H_X`` as:

```math
\\begin{aligned}
H_Z = [B^T | A^T] = [C B C | C A C] = C [B | A] \\begin{bmatrix} 0 & C \\\\ C & 0 \\end{bmatrix} = C H_X \\begin{bmatrix} 0 & C \\\\ C & 0 \\end{bmatrix}.
\\end{aligned}
```

Since both left and right multiplications are with invertible matrices (`C` is
self-inverse and permutation), it follows that:

```math
\\begin{aligned}
\text{rank}(H_X) = \\text{rank}(H_Z)
\\end{aligned}
```

This yields a new expression for the number of encoded qubits:

```math
\\begin{aligned}
k = n - \\text{rank}(H_X) - \\text{rank}(H_Z) = n - 2 \\cdot \\text{rank}(H_Z).
\\end{aligned}
```

And further, using:

```math
\\begin{aligned}
\\ker(H_Z^T) = \\ker(A) \\cap \\ker(B),
\\end{aligned}
```

we get:

```math
\\begin{aligned}
k = 2 \\cdot \\dim(\\ker(A) \\cap \\ker(B)).
\\end{aligned}
```
"""
function bivariate_bicycle_code_k(c::AbstractCSSCode)
    Hx = parity_matrix_x(c)
    n = size(Hx,2)÷2
    A = matrix(GF(2), Hx[:,1:n])
    B = matrix(GF(2), Hx[:,n+1:end])
    F = GF(2)
    hA = hom(free_module(F, size(A, 1)), free_module(F, size(A, 2)), A)
    hB = hom(free_module(F, size(B, 1)), free_module(F, size(B, 2)), B)
    ans = kernel(hA)[1] ∩ kernel(hB)[1]
    k = 2*size(map(domain(hA), gens(ans[1])), 1)
    return k
end
