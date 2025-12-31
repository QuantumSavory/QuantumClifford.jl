"""
    $TYPEDEF

Repetition code is a classical error correction code that repeats the same bit `n` times.

### Fields
    $TYPEDFIELDS
"""
struct RepCode <: AbstractCECC
    """The number of times to repeat the bit."""
    n::Int
end

function parity_matrix(c::RepCode)
    n = c.n
    I = [i for i in 1:n for δ in (0,1)]
    J = [(i+δ-1)%n+1 for i in 1:n for δ in (0,1)]
    V = true
    sparse(I,J,V,n,n)
end

code_n(c::RepCode) = c.n
code_s(c::RepCode) = c.n
