using Test
using QuantumClifford.ECC.QECCore
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: BCH
using InteractiveUtils
using SparseArrays

import Nemo: GF, matrix, rank, finite_field, polynomial_ring
import Hecke: group_algebra, abelian_group, gens, representation_matrix
import LinearAlgebra

include("test_ecc_base.jl")

# generate instances of all implemented codes to make sure nothing skips being checked

# Goppa Codes
m₁ = 4
t₁ = 2
F, α = finite_field(2, m₁, :α)
R, x = polynomial_ring(F, :x)
g₁ = x^2 + x + α^3
L₁ = [α^i for i in 2:13]

m₂ = 6
F, α = finite_field(2, m₂, :α)
R, x = polynomial_ring(F, :x)
t₂ = 9
g₂ = x^t₂ + 1

m₃ = 3
F, α = finite_field(2, m₃, :α)
R, x = polynomial_ring(F, :x)
t₃ = 2
g₃ = x^t₃ + x + 1

m₄ = 4
t₄ = 2
F, α = finite_field(2, m₄, :α)
R, z = polynomial_ring(F, :z)
g₄ = z^2 + α^7*z + 1
L₄ = [α^i for i in 2:13]

# Lifted Codes
l₁ = 12; GA₁ = group_algebra(GF(2), abelian_group(l)); x = gens(GA)[]
B₁ = reshape([1 + x + x^3 + x^6], (1, 1))
c = LiftedCode(B, repr = representation_matrix)

base_matrix₂ = [0 0 0 0; 0 1 2 5; 0 6 3 1]; l₂ = 3;

const classical_code_instance_args = Dict(
    :RepCode => [3, 4, 5, 6, 7, 8, 9, 10],
    :BCH => [(3, 1), (3, 2), (4, 1), (4, 1), (5, 1), (5, 2), (6, 1), (6, 2)],
    :ReedMuller => [(1, 3), (1, 4), (2, 3), (1, 5), (2, 4), (2, 5), (3, 5)],
    :RecursiveReedMuller => [(1, 3), (1, 4), (2, 3), (1, 5), (2, 4), (2, 5), (3, 5)],
    :Golay => [(23), (24)],
    :Hamming => [2, 3, 4, 5, 6, 7, 8],
    :GallagerLDPC => [(3, 3, 4), (3, 4, 5), (4, 5, 7), (4, 6, 7)],
    :GoppaCode => [(m₁, t₁, g₁, L₁), (m₂, t₂, g₂), (m₃, t₃, g₃), (m₄, t₄, g₄, L₄)],
    :LiftedCode => [(B₁, repr = representation_matrix), (base_matrix₂, l₂), (base_matrix₂, 5), (base_matrix₂, 7)]
)

function all_testable_classical_code_instances(; maxn=nothing)
    codeinstances = []
    i = 1
    classical_code_args = copy(classical_code_instance_args)
    for t in concretesubtypes(AbstractCECC)
        for args in pop!(classical_code_args, nameof(t), [])
            codeinstance = t(args...)
            !isnothing(maxn) && code_n(codeinstance) > maxn && continue
            push!(codeinstances, codeinstance)
            #@show i, t, code_n(codeinstance), code_k(codeinstance), code_n(codeinstance)-code_k(codeinstance)
            i += 1
        end
    end
    @test isempty(classical_code_args) # if this fails, then some code instances were not tested
    return codeinstances
end
