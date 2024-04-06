"""The family of Reed-Solomon codes, as discovered by Reed and Solomon in their 1960 paper [reed1960polynomial](@cite). 

You might be interested in consulting [geisel1990tutorial](cite), [wicker1999reed](@cite), [sklar2001reed](@cite), [berlekamp1978readable](cite) and [https://youtu.be/K26Ssr8H3ec?si=QOeohq_6I0Oyd8qu](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/reed_solomon)
"""
struct ReedSolomon <: ClassicalCode
    n::Int
    k::Int

    function ReedSolomon(n, k)
        if n < 0 || k < 0 || n > 500
            throw(ArgumentError("Invalid parameters: n and k must be non-negative and n > 500 in order to obtain a valid code and to remain tractable"))
        end
            new(n, k)
    end
end

#section 3.3.2 [geisel1990tutorial](cite) for constructing custom message polynomial m(x)
function _message_polynomial_rs(k, a, an, x, positions)
    message = 0*x
    for pos in positions
        if pos <= k
            message += a^(an)*x^(pos)
        else
            throw(DomainError("Invalid bit positions: The number of bit positions [0th,..., kth] to assign coefficient values should be <= k"))
        end
    end
    return message
end

#section 3.3.2 [geisel1990tutorial](cite) for constructing parity_check polynomial ck(x) given m(x) and g(x)
function _parity_check_polynomial_rs(n, k, x, message, gx)
    x_nminusk = x^(n - k)
    ck = mod(x_nminusk*message, gx)
    return ck
end

#section 3.3.2 [geisel1990tutorial](cite) for constructing codeword polynomial c(x) given m(x) and ck(x)
function _codeword_polynomial_rs(message, ck)
    cx = message + ck
    return cx
end

function generator_polynomial(rs::ReedSolomon)
    r = ceil(Int, log2(rs.n + 1))
    t = div(rs.n - rs.k, 2)
    GF2ͬ, a = finite_field(2, r, "a")
    P, x = GF2ͬ[:x]
    poly_zeros = 2*t
    gx = x - a^poly_zeros
    for i in poly_zeros:(poly_zeros + 2*t - 1)
       gx *= (x - a^(poly_zeros + i))
    end
    return gx
end
