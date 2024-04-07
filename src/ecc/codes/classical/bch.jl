"""The family of Bose–Chaudhuri–Hocquenghem (BCH) codes, as discovered in 1959 by Alexis Hocquenghem [hocquenghem1959codes](@cite), and independently in 1960 by Raj Chandra Bose and D.K. Ray-Chaudhuri [bose1960class](@cite).

You might be interested in consulting [bose1960further](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/q-ary_bch)
"""
struct BCH <: ClassicalCode
    n::Int
    t::Int

    function BCH(n, t)
        if n < 0 || t < 1 || n > 500 || t > 500 
            throw(ArgumentError("Invalid parameters: n and t must be non-negative and < 500 in order to obtain a valid code and to remain tractable"))
        end
            new(n, t)
    end
end

function generator_polynomial(rs::BCH)
    r = ceil(Int, log2(rs.n + 1))
    GF2ͬ, a = finite_field(2, r, "a")
    GF2x, x = GF(2)["x"]
    minimal_poly = FqPolyRingElem[]
    for i in 1:(2*rs.t - 1)
        if i % 2 != 0
             minimal_poly = [minimal_poly; minpoly(GF2x, a^i)]
        end 
    end
    gx = lcm(minimal_poly)
    return gx
end
