"""The family of Reed-Solomon codes, as discovered by Reed and Solomon in their 1960 paper [reed1960polynomial](@cite). 

You might be interested in consulting [wicker1999reed](@cite), [sklar2001reed](@cite), [berlekamp1978readable](cite) and [https://youtu.be/K26Ssr8H3ec?si=QOeohq_6I0Oyd8qu](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/reed_solomon)
"""
struct ReedSolomon <: ClassicalCode
    n::Int
    k::Int
    j::Int

    function ReedSolomon(n, k, j)
        if n < 0 || k < 0 || j < 0 || n >= 500
            throw(ArgumentError("Invalid parameters: n, k, j must be non-negative and n >= 500 in order to obtain a valid code and to remain tractable"))
        end 
            new(n, k, j)
    end
end

function generator_matrix(rs::ReedSolomon)
    r = ceil(Int, log2(rs.n + 1))
    t = div(rs.n - rs.k, 2)
    GF2ͬ, a = finite_field(2, r, "a")
    P, x = GF2ͬ[:x]
    polynomial = x - a^rs.j
    for i in 1:(2*t - 1)
       polynomial *= (x - a^(rs.j + i))
    end
    return polynomial
end
