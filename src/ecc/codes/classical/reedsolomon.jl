import Base: length, iterate, ==, <<, +, *, ÷, %, copy, zero, eltype

# Define a structure for Reed-Solomon code

abstract type ClassicalCode end

struct ReedSolomon{T<:Integer} <: ClassicalCode
    m::Int
    e::Int
end

# Constructor for Reed-Solomon
function ReedSolomon(m::Int, e::Int)
    return ReedSolomon{Int}(m, e)
end

# Define a structure for polynomials
struct Poly{T<:Integer}
    coeff::Vector{T}
end

# Precompute powers of 2 in GF(2^8)
const GF256 = let
    tab = ones(UInt8, 255)
    v = 1
    for i in 2:255
        v <<= 1
        if v > 255
            v = xor(v, 285)
        end
        tab[i] = v
    end
    tab
end

# Precompute inverse of powers of 2 in GF(2^8)
const INV_GF256 = let
    tab = Vector{UInt8}(undef, 255)
    tab[GF256] .= 0:254
    tab
end

# Compute power of 2 in GF(2^8)
GF256_POW2(n::Integer) = GF256[mod(n, 255) + 1]
GF256_POW2(n::UInt8) = GF256[n + 0x1 + (n == 0xff)]

# Compute logarithm of a number in GF(2^8)
GF256_LOG2(n::Integer) = INV_GF256[n]

# Compute inverse in GF(2^8)
G256INV(a::Integer) = GF256_POW2(0xff-GF256_LOG2(a))

# Precompute multiplication table in GF(2^8)
const GF256_MULT = [iszero(i * j) ? 0x0 : GF256_POW2(Int(GF256_LOG2(i)) + GF256_LOG2(j)) for i in 0:255, j in 0:255] 

# Precompute division table in GF(2^8)
const GF256_DIV = [iszero(i) ? 0x0 : GF256_POW2(Int(GF256_LOG2(i)) - GF256_LOG2(j)) for i in 0:255, j in 1:255]

# Multiplication in GF(2^8)
function G256MULT(a::Integer, b::Integer)
    return GF256_MULT[a + 1, b + 1]
end

# Division in GF(2^8)
function G256DIV(a::Integer, b::Integer)
    iszero(b) && throw(DivideError())
    return GF256_DIV[a + 1, b]
end

# Check if polynomial is zero
iszeropoly(p::Poly) = all(iszero, p)

# Remove trailing zeros from polynomial
rtrailingzeros(p::Poly) = rtrailingzeros!(copy(p))
function rtrailingzeros!(p::Poly{T}) where T
    iszeropoly(p) && return zero(Poly{T})
    deleteat!(p.coeff, findlast(!iszero, p.coeff) + 1:length(p))
    return p
end

# Pad polynomial with zeros to a specified length
rpadzeros(p::Poly, n::Int) = rpadzeros!(copy(p), n)
function rpadzeros!(p::Poly{T}, n::Int) where T
    length(p) > n && throw("rpadzeros: length(p) > n")
    append!(p.coeff, zeros(T, n - length(p)))
    return p
end

# Get degree of polynomial
function deg(p::Poly)
    iszeropoly(p) && return -1
    return findlast(!iszero, p.coeff) - 1
end

# Zero polynomial constructor
zero(::Type{Poly{T}}) where T = Poly{T}(zeros(T, 1))

# Unit polynomial constructor
unit(::Type{Poly{T}}) where T = Poly{T}(ones(T, 1))

# Copy polynomial
copy(p::Poly{T}) where T = Poly{T}(copy(p.coeff))

# Get length of polynomial
length(p::Poly) = length(p.coeff)

# Iterator for polynomial
iterate(p::Poly) = iterate(p.coeff)
iterate(p::Poly, i) = iterate(p.coeff, i)

# Element type of polynomial
eltype(::Poly{T}) where T = T

# Check equality of polynomials
==(a::Poly, b::Poly)::Bool = iszeropoly(a + b)

# Left shift for polynomial
<<(p::Poly{T}, n::Int) where T = Poly{T}(vcat(zeros(T, n), p.coeff))

# Addition of polynomials
function +(a::Poly{T}, b::Poly{T}) where T
    l, o = max(length(a), length(b)), zero(T)
    return Poly{T}([xor(get(a.coeff, i, o), get(b.coeff, i, o)) for i in 1:l])
end

# Scalar multiplication with polynomial
*(a::Integer, p::Poly{T}) where T = Poly{T}(G256MULT.(a, p.coeff))
*(p::Poly{T}, a::Integer) where T = Poly{T}(G256MULT.(a, p.coeff))

# Multiplication of polynomials
function *(a::Poly{T}, b::Poly{T}) where T
    polyproduct = Poly(zeros(T, length(a) + length(b) - 1))
    for (i, c1) in enumerate(a.coeff), (j, c2) in enumerate(b.coeff)
        @inbounds polyproduct.coeff[i + j - 1] ⊻= G256MULT(c2, c1)
    end
    return polyproduct
end

# Construct polynomial x^n
polynomial_xn(::Type{T}, n::Int) where T = Poly{T}(push!(zeros(T, n), one(T)))
polynomial_xn(n::Int) = polynomial_xn(UInt8, n)

# Construct polynomial x^n padded to a specific length
function polynomial_xn(::Type{T}, n::Int, pad::Int) where T
    pad > n || throw(DomainError("pad length $pad should be greater than the power $n"))
    coef = zeros(T, pad)
    coef[n+1] = one(T)
    return Poly{T}(coef)
end
polynomial_xn(n::Int, pad::Int) = polynomial_xn(UInt8, n, pad)

# Euclidean division of polynomials
euclidean_division(f::Poly, g::Poly) = euclidean_division!(copy(f), copy(g))
function euclidean_division!(f::Poly{T}, g::Poly{T}) where T
    g, f = rtrailingzeros!(g), rtrailingzeros!(f)
    iszeropoly(g) && throw(DivideError())
    fcoef, gcoef, lf, lg = f.coeff, g.coeff, length(f), length(g)
    gn = last(gcoef)
    if isone(lg)
        fcoef .= G256DIV.(fcoef, gn)
        return f, zero(Poly{T})
    end
    quodeg = lf - lg
    quodeg < 0 && return zero(Poly{T}), f
    @inbounds for i in 0:quodeg
        leadterm = G256DIV(fcoef[end-i], gn)
        for (j, c) in enumerate(gcoef)
            fcoef[quodeg - i + j] ⊻= G256MULT(leadterm, c)
        end
        fcoef[end-i] = leadterm
    end
    quo = Poly{T}(fcoef[end-quodeg:end])
    deleteat!(fcoef, lf-quodeg:lf)
    return quo, f
end

# Division operator for polynomials
÷(f::Poly, g::Poly) = first(euclidean_division(f, g))

# Modulus operator for polynomials
%(f::Poly, g::Poly) = last(euclidean_division(f, g))

# Generate generator polynomial for Reed-Solomon code
generator(n::Int) = generator(UInt8, n)
generator(::Type{T}, n::Int) where T = prod([Poly{T}([GF256[i], one(T)]) for i in 1:n])

# Encode a message polynomial for Reed-Solomon code
function encodepoly(msgpoly::Poly{T}, n::Int) where T
    f = msgpoly << n
    f + f % generator(T, n)
end

# Get error code from received polynomial for Reed-Solomon code
geterrcode(f::Poly{T}, n::Int) where T = rpadzeros!(f << n % generator(T, n), n)

# Multiply matrices in GF(2^8)
function G256multmat(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    (n, m), (p, q) = size(A), size(B)
    m == p || throw(DimensionMismatch("A and B have different size"))
    C = Array{T}(undef, n, q)
    @inbounds for i in 1:n, j in 1:q
        C[i, j] = @views reduce(xor, G256MULT.(A[i, :], B[:, j]))
    end
    return C
end

# Multiply matrix with vector in GF(2^8)
function G256multmat(A::AbstractMatrix{T}, b::AbstractVector{T}) where T
    (n, m) = size(A)
    m == length(b) || throw(DimensionMismatch("A and b have different size"))
    c = Vector{T}(undef, n)
    @inbounds for i in 1:n
        c[i] = @views reduce(xor, G256MULT.(A[i, :], b))
    end
    return c
end

# Generate generator matrix for Reed-Solomon code
function generator_matrix(rs::ReedSolomon)
    G = Array{UInt8}(undef, rs.m + rs.e, rs.m)
    g = generator(UInt8, rs.e)
    @inbounds for i in 1:rs.m
        xi = polynomial_xn(UInt8, i-1, rs.m) << rs.e
        G[:, i] = (xi + xi % g).coeff
    end
    G = G[end:-1:1, end:-1:1]
    G = map(x -> x == 0 ? 0 : 1, G)
    return G
end
