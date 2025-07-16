"""
    $TYPEDEF

The family of `[[2ʲ, 2ʲ - j - 2, 3]]` Gottesman codes, also known as quantum Hamming codes, as described in [Gottesman's 1997 PhD thesis](@cite gottesman1997stabilizer) and in [gottesman1996class](@cite).

You might be interested in consulting [yu2013all](@cite) and [chao2018quantum](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_hamming)

### Fields
    $TYPEDFIELDS
"""
struct Gottesman <: AbstractQECC
    """The parameter `j` of the code."""
    j::Int

    function Gottesman(j)
        (j >= 3 && j < 21) || error("In `Gottesman(j)`, `j` must be ≥  3 in order to obtain a valid code and < 21 to remain tractable")
        new(j)
    end
end

code_n(c::Gottesman) = 2^c.j

function parity_matrix(c::Gottesman)
    j = c.j
    s = j+2
    n = 2^j

    H = zeros(Bool, s, 2*n)
    for i in 1:n
        H[1, i] = true
        H[2, i+n] = true
    end
    for i in 0:n-1 # column of H, corresponds to a single qubit error that is detectable)
        Xⁱ = i
        Zⁱ = i÷2
        jeven = j%2 == 0
        ieven = i%2 == 0
        if (jeven && ieven) || (!jeven && ieven && i < n÷2) || (!jeven && !ieven && i ≥ n÷2)
            Zⁱ = ~Zⁱ
        end
        for b in 0:j-1 # which check to consider (row of H), also which bit to extract
            H[s-b,i+1] = isone((Zⁱ>>b)&0x1)
            H[s-b,i+n+1] = isone((Xⁱ>>b)&0x1)
        end
    end
    H
end
