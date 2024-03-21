"""The family of Reed-Muller codes, as discovered by Muller in his 1954 paper [muller1954application](@cite) and Reed who proposed the first efficient decoding algorithm [reed1954class](@cite).

You might be interested in consulting [raaphorst2003reed](@cite), [abbe2020reed](@cite), and [djordjevic2021quantum](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/reed_muller)
"""

abstract type ClassicalCode end

struct ReedMuller <: ClassicalCode
    r::Int
    m::Int

    function ReedMuller(r, m)
        if r < 0 || m < 1
            throw(ArgumentError("Invalid parameters: r must be non-negative and m must be positive"))
        end
        new(r, m)
    end
end

function variables_xi(m, i)
    return repeat([fill(1, 2^(m - i - 1)); fill(0, 2^(m - i - 1))], outer = 2^i)
end

function vmult(vecs...)
    return [reduce(*, a, init=1) for a in zip(vecs...)]
end

function parity_checks(c::ReedMuller)
    r=c.r
    m=c.m
    xi = [variables_xi(m, i) for i in 0:m - 1]
    row_matrices = [reduce(vmult, [xi[i + 1] for i in S], init = ones(Int, 2^m)) for s in 0:r for S in combinations(0:m - 1, s)]
    rows = length(row_matrices)
    cols = length(row_matrices[1])
    H = reshape(vcat(row_matrices...), cols, rows)'
end