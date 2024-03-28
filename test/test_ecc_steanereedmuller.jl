using Test
using Combinatorics
using LinearAlgebra
using QuantumClifford
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
    for (t, r) in [(1, 3), (2, 4)]
        stab = parity_checks(SteaneReedMuller(t, r))
        H = stab_to_gf2(stab)
        Hx, Hz = split(H)
        computed_distance = min_distance(Hx, Hz)
        expected_distance = 2^t + 2^(t - 1)
        @test computed_distance == expected_distance
    end
end

function generate_gH(t, r, Hx, Hz)

   gr = size(parity_checks(ReedMuller(t, r)), 1)
   gHx = Hx[1:1, :]
   gHx = vcat(gHx, zeros(1, size(Hx, 2)))
   gHx = vcat(gHx, Hx[2:gr, :])
   gHx = convert(Matrix{Int}, gHx)
   sa = div(length(Hx[2, :]), 4)
   sr = circshift(Hx[2, :], -sa)
   sr = reshape(sr, 1, length(sr))
   gHz = vcat(Hz[gr:1, :], zeros(1, size(Hz, 2)))
   gHz = vcat(gHz, Hz[gr+1:gr+1, :]) 
   gHz = vcat(gHz, Hz[gr+3:2*gr, :]) 
   gHz = vcat(gHz, sr)
   gHz = convert(Matrix{Int}, gHz)
   gH = hcat(gHx, gHz)
   return gH
end

@testset "Test @t = 1, r = 3, 4, 5 Gotteman family in Steane's matrix convention [steane1999quantum](@cite)" begin
    
    stab = parity_checks(SteaneReedMuller(1, 3))
    H = stab_to_gf2(stab)
    Hx, Hz = split(H)
    gH = generate_gH(1, 3, Hx, Hz)

    #Gottesman [[8, 3, 3]]
    @test gH  == [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
                  0 0 0 0 1 1 1 1 0 0 1 1 0 0 1 1;
                  0 0 1 1 0 0 1 1 0 1 0 1 0 1 0 1;
                  0 1 0 1 0 1 0 1 0 0 1 1 1 1 0 0] 


    stab = parity_checks(SteaneReedMuller(1, 4))
    H = stab_to_gf2(stab)
    Hx, Hz = split(H)
    gH = generate_gH(1, 4, Hx, Hz)

    #Gottesman [[16, 10, 3]]
    @test gH == [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 
                 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
                 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1; 
                 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1;
                 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0]

    stab = parity_checks(SteaneReedMuller(1, 5))
    H = stab_to_gf2(stab)
    Hx, Hz = split(H)
    gH = generate_gH(1, 5, Hx, Hz)
   
    #Gottesman [[32, 25, 3]]
    @test gH == [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; 
                 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
                 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
                 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1;
                 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1; 
                 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0]
end