"""
A bivariate bicycle (BB) quantum LDPC code was introduced by Bravyi et al. in their 2024 paper [bravyi2024high](@cite). This code uses identity and cyclic shift matrices. Define `Iₗ` as the `l × l` identity matrix and `Sₗ` as the cyclic shift matrix of the same size, where each row of `Sₗ` has a single '1' at the column `(i + 1) mod l`.

The matrices `x = Sₗ ⊗ Iₘ` and `y = Iₗ ⊗ Sₘ` are used. The BB code is represented by matrices `A` and `B`, defined as: `A = A₁ + A₂ + A₃` and `B = B₁ + B₂ + B₃`. The addition and multiplication operations on binary matrices are performed modulo 2. The check matrices are: `Hx = [A|B]` and `Hz = [B'|A']`. Both `Hx` and `Hz` are `(n/2)×n` matrices.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/qcga).
"""
struct BBQLDPC <: AbstractECC
    l::Int
    m::Int
    A::Vector{Int}
    B::Vector{Int}
end

function iscss(::Type{BBQLDPC})
    return true
end

function _AB(c::BBQLDPC)
    a₁,a₂,a₃ = c.A[1],c.A[2],c.A[3]
    b₁,b₂,b₃ = c.B[1],c.B[2],c.B[3]
    Iₗ  = Matrix{Bool}(LinearAlgebra.I,c.l,c.l)
    Iₘ = Matrix{Bool}(LinearAlgebra.I,c.m,c.m)
    x  = Dict{Bool, Matrix{Bool}}()
    y  = Dict{Bool, Matrix{Bool}}()
    x  = Dict(i => kron(circshift(Iₗ,(0,i)),Iₘ) for i in 0:(c.l))
    y  = Dict(i => kron(Iₗ,circshift(Iₘ,(0,i))) for i in 0:(c.m))
    A  = mod.(x[a₁]+y[a₂]+y[a₃],2)
    B  = mod.(y[b₁]+x[b₂]+x[b₃],2)
    return A, B
end

function parity_checks(c::BBQLDPC)
    A, B = _AB(c)
    Hx = hcat(A,B)
    Hz = hcat(B',A')
    H  = CSS(Hx,Hz)
    Stabilizer(H)
end

code_n(c::BBQLDPC) = 2*c.l*c.m
code_k(c::BBQLDPC) = code_n(c) - LinearAlgebra.rank(matrix(GF(2), parity_checks_x(c))) - LinearAlgebra.rank(matrix(GF(2), parity_checks_z(c)))

parity_checks_x(c::BBQLDPC) = hcat(_AB(c)...)

parity_checks_z(c::BBQLDPC) = hcat(_AB(c)[2]', _AB(c)[1]')
