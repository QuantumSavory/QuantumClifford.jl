@testitem "Classical" begin
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
        for stabilizer in group
            apply!(s, stabilizer)
            @test s == s_test
        end
    end
    #Test minimal_generating_set
    for n in [1, test_sizes...]
        s = random_stabilizer(n)
        group = groupify(s)
        gen_set = minimal_generating_set(Stabilizer(group))
        @test length(group) == 2^(length(gen_set))
        new_group = groupify(gen_set)
        s1, _, r = canonicalize!(Stabilizer(group), ranks = true)
        s2, _, r = canonicalize!(Stabilizer(new_group), ranks=true)
        @test group[1:r, :] == new_group[1:r, :]
    end
    #Test logical_operator_canonicalize
    for n in [1, small_test_sizes...]
        t = zero(QuantumClifford.Tableau, rand(1:(2*n)), n)
        for i in eachindex(t) t[i] = random_pauli(n) end 
        loc = logical_operator_canonicalize(t).tab 
        index = 0
        for i in range(1, stop=length(loc)+1, step=2)
            if i + 1 > length(loc)|| comm(loc[i], loc[i+1]) == 0x0 
                index = i # index to split loc into non-commuting pairs and commuting operators
                break
            end
        end
        for i in range(1, length(loc))
            for j in range(1, length(loc))
                if i % 2 == 1
                    adj = j == i + 1
                else
                    adj = j == i - 1
                end
                if (i < index || j < index) && adj 
                    @test comm(loc[i], loc[j]) == 0x01
                elseif i != j
                    @test comm(loc[i], loc[j]) == 0x0
                end
            end
        end   
    end
    #Test commutativise
    for n in [1, small_test_sizes...]
        t = zero(QuantumClifford.Tableau, rand(1:(2*n)), n)
        for i in eachindex(t) t[i] = random_pauli(n) end
        c, d = commutavise(t)
        for i in c
            for j in c
                @test comm(i, j) == 0x0
            end
        end
        for i in d
            for p in c 
                @test p[i] != (true, true)
            end
        end
        loc1= delete_columns(c, d)
        loc2 = logical_operator_canonicalize(t).tab
        for i in eachindex(delete_columns(c, d))
           for j in eachindex(loc1[i]) @test loc1[i][j] ==loc2[i][j] end
        end
    end
    #Test embed
    for n in [1,2,3,4,5]
        t = zero(QuantumClifford.Tableau, 2*n, n)
        for i in eachindex(t) t[i] = random_pauli(n) end
        e, d2, d1 = embed(t)
        s = Stabilizer(groupify(e))
        @test 2^(nqubits(s)) == length(s) #assumes commutativise works
        #find oringal tableau from embedded state, ignoring phases
        inverted = delete_columns(Stabilizer(normalizer(delete_columns(Stabilizer(e), d2).tab)), d1)
        original = Stabilizer(groupify(Stabilizer(t)))
        canonicalize!(inverted)
        canonicalize!(original)
        @test inverted.tab.xzs == original.tab.xzs
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
        paulis = QuantumClifford.pauligroup(n, phases=false)
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
    # #Test centralizer
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
        #Test contractor
        for n in [1, test_sizes...]
            s = random_stabilizer(n)
            subset = []
            for i in 1:nqubits(s) #create a random subset
                if rand(1:2) == 1 push!(subset, i) end
            end
            c = contractor(s, subset)
            count = 0
            for stabilizer in s 
                contractable = true
                for i in subset
                    if stabilizer[i] != (false, false) contractable = false end
                end
                if contractable count+=1 end
            end
            if length(c[1]) > 0
            @test count == length(c)
                for contracted in c
                    p = zero(PauliOperator, nqubits(s))
                    index = 0
                    for i in 1:nqubits(s)
                        if !(i in subset) 
                            index+=1 
                            p[i] = contracted[index] 
                        end
                    end
                    @test p in s || -1* p in s || 1im * p in s || -1im * p in s
                end
            else
                @test count ==0
            end
        end
    end
end