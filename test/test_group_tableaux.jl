using Random
using QuantumClifford

# Including sizes that would test off-by-one errors in the bit encoding.
test_sizes = [1,2,3,4,5,7,8,9,15,16,17] # Zero function(in groupify) slows down around 2^30(n=30),eventually breaks
# Zero function(in groupify) slows down around 2^30(n=30),eventually breaks
small_test_sizes = [1,2,3,4,5,7,8] # Pauligroup slows around n =9

@test_set "Group Tableaux" begin 
    #test groupify
    for n in [1,test_sizes...]
        s = random_stabilizer(n)
        s_test = copy(s)
        group = groupify(s)
        @test length(group) == 2^n
        unchanged = true
        for stabilizer in group
            apply!(s, stabilizer)
            if !( s== s_test)
                unchanged = false
            end 
            @test unchanged == true
        end 
    end
    #test get_generating_set
    for n in [1,small_test_sizes...]
        s = zero(Stabilizer, rand(1:(2*n)), n)
        for i in 1:length(s)
            s[i] = random_pauli(n)
        end
        group = groupify(s)
        gen_set = get_generating_set(group)
        new_group = groupify(gen_set)
        canonicalize!(group)
        canonicalize!(new_group)
        @test group == new_group
    end
    #test logical operator canonicalize
    for n in [1, test_sizes...]
        s = zero(Stabilizer, n, n)
        for i in 1:n
            s[i] = random_pauli(n)
        end
        loc, r = QuantumClifford.logical_operator_canonicalize(s)
        #@test groupify(s) == groupify(loc) current groupify only works for fully commutative generating sets
        index = 0
        for i in range(1, length(loc)-1)
            if comm(loc[i], loc[i+1]) == 0x0
                index  = i # index to split loc into non-commuting pairs and commuting operators
                break
            end
            i = i+1
        end
        for i in range(1, length(loc))
            for j in range(1, length(loc))
                if i % 2 == 1
                    adj = j == i + 1
                else
                    adj = j == i - 1
                end
                if i < index && adj 
                    @test comm(loc[i], loc[j]) == 0x0
                else
                    @test comm(loc[i], loc[j]) == 0x01
                end
    
            end
        end
    end
    #test normalize
    for n in [1, small_test_sizes...] # pauligroup is very slow at n=14
        s = random_stabilizer(n)
        normalized = normalize(s)
        stabilizers = pauligroup(n, true)
        for n_stabilizer in normalized
            for stabilizer in s
               @test comm(n_stabilizer, stabilizer) == 0x0
            end
        end
    end
    for stabilizer in stabilizers
        commutes = true
        for n_stabilizer in normalized
            if !(comm(n_stabilizer, stabilizer) == 0x0)
                commutes = false
            end
        end
        @test !commutes  || stabilizer in normalized
    end
    #test center
    for n in [1, test_sizes...]
        s = random_stabilizer(n)
        c = center(s)
        for c_stabilizer in c
            for stabilizer in s
               @test comm(c_stabilizer, stabilizer) == 0x0
            end
        end
        for stabilizer in s
            commutes = true
            for c_stabilizer in c
                @test comm(c_stabilizer, stabilizer) == 0x0
            end
            @test !commutes  || st in c
        end
            
    end
end
