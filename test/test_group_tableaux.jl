using Test

using Random
using QuantumClifford

# Including sizes that would test off-by-one errors in the bit encoding.
test_sizes = [1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17] 
# Zero function(in groupify) slows down around 2^30(n=30),eventually breaks
small_test_sizes = [1, 2, 3, 4, 5, 7] # Pauligroup slows around n = 8

@testset "group_tableaux" begin
    #Test groupify
    for n in [1, test_sizes...]
        s = random_stabilizer(n)
        s_test = copy(s)
        group = groupify(s)
        @test length(group) == 2^n
        unchanged = true
        for stabilizer in group
            apply!(s, stabilizer)
            if !(s == s_test)
                unchanged = false
            end
            @test unchanged == true
        end
    end
    #Test minimal_generating_set
    for n in [1, small_test_sizes...]
        s = random_stabilizer(n)
        group = groupify(s)
        gen_set = minimal_generating_set(Stabilizer(group))
        new_group = groupify(gen_set)
        canonicalize!(Stabilizer(group))
        canonicalize!(Stabilizer(new_group))
        @test group == new_group
        s = zero(Stabilizer, rand(1:(2*n)), n)
        for i in 1:length(s)
            s[i] = random_pauli(n)
        end
        gen_set = minimal_generating_set(s)
        new_group = groupify(s)
        for operator in s
            @test operator in new_group
        end
    end
    #Test pauligroup
    for n in [1, small_test_sizes...] 
        @test length(QuantumClifford.pauligroup(n, phases=false)) == 4^n
        @test length(QuantumClifford.pauligroup(n, phases=true)) == 4^(n+1)
    end
    #Test normalizer
    for n in [1, small_test_sizes...] # pauligroup is very slow at n=14
        t = zero(QuantumClifford.Tableau, rand(1:(2*n)), n)
        for i in eachindex(t) t[i] = random_pauli(n) end
        normalized = normalizer(t)
        paulis = QuantumClifford.pauligroup(n, phases=true)
        for n_pauli in normalized
            for t_pauli in t
                @test comm(n_pauli, t_pauli) == 0x0
            end
        end
        for pauli in paulis
            commutes = true
            for t_pauli in t
                if comm(t_pauli, pauli) == 0x01
                    commutes = false
                end
            end
           @test (!commutes) || (pauli in normalized)
        end
    end
    #Test centralizer
    for n in [1, test_sizes...]
        t = zero(QuantumClifford.Tableau, rand(1:(2*n)), n)
        for i in eachindex(t) t[i] = random_pauli(n) end
        c = centralizer(t)
        for c_pauli in c
            for t_pauli in t
                @test comm(c_pauli, t_pauli) == 0x0
            end
        end
        for pauli in t
            commutes = true
            for t_pauli in t
                if comm(t_pauli, pauli)==0x01 commutes = false end
            end
            @test !commutes || pauli in c
        end

    end
    #Test contract
    for n in [1, test_sizes...]
        s = random_stabilizer(n)
        g = groupify(s)
        subset = []
        for i in 1:nqubits(s) #create a random subset
            if rand(1:2) == 1 push!(subset, i) end
        end
        c = contracter(s, subset)
        count = 0
        for stabilizer in g 
            contractable = true
            for i in subset
                if stabilizer[i] != (false, false) contractable = false end
            end
            if contractable count+=1 end
        end
        @test count == length(c)
        for contracted in c
            p = zero(PauliOperator, nqubits(s))
            index = 0
            for i in 1:nqubits(g)
                if !(i in subset) 
                    index+=1 
                    p[i] = contracted[index] 
                end
            end
            @test p in g || -1* p in g || 1im * p in g || -1im * p in g
        end
    end
end
