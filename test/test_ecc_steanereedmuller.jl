using Test
using Combinatorics
using LinearAlgebra
using QuantumClifford
using QuantumClifford: mul_left!
using QuantumClifford.ECC
using QuantumClifford.ECC: AbstractECC, SteaneReedMuller, ReedMuller

function split(H::Matrix{Bool})
    n = size(H, 2)
    m = div(n, 2)
    Hx = H[:, 1:m]
    Hz = H[:, m+1:end]
    return Hx, Hz
end

function min_distance(Hx::Matrix{Bool}, Hz::Matrix{Bool})
    w = typemax(Int)
    for (x, z) in zip(eachrow(Hx), eachrow(Hz))
        r = map(|, x, z)
        w = min(w, sum(r))
    end
    return w 
end

@testset "Test QRM(t, r) Matrix Minimum Distance" begin
        stab = parity_checks(SteaneReedMuller(2, 4))
        H = stab_to_gf2(stab)
        Hx, Hz = split(H)
        computed_distance = min_distance(Hx, Hz)
        expected_distance = 2^2 + 2^(2 - 1)
        @test computed_distance == expected_distance
end

@testset "Gottesman codes should correct all single-qubit errors" begin
    for j in 3:7
        H = parity_checks(SteaneReedMuller(1, j))
        syndromes = Set([]) 
        for error_type in (single_x, single_y, single_z)
            for bit_index in 1:nqubits(H)
                syndrome = comm(H, error_type(nqubits(H), bit_index))
                @test any(==(0x1), syndrome)
                push!(syndromes, syndrome)
            end
        end
        @test length(syndromes) == 3*nqubits(H)
    end
end

@testset "Test @t =1, r =3  Gotteman family in Steane's matrix convention [steane1999quantum](@cite)" begin
    stab = parity_checks(SteaneReedMuller(1, 3))
    H = stab_to_gf2(stab)

    #Gottesman [[8,3,3]], Parity Check Tableau is taken from 'Table 15' of Steane's paper [steane1999quantum](@cite)
    @test H  ==  [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
                  0 0 0 0 1 1 1 1 0 0 1 1 0 0 1 1;
                  0 0 1 1 0 0 1 1 0 1 0 1 0 1 0 1;
                  0 1 0 1 0 1 0 1 0 0 1 1 1 1 0 0] 
end
