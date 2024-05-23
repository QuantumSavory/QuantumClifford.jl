using Random
using QuantumClifford

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.


@test_set "Group Tableaux" begin 
    for n in [1,test_sizes...]
        s = random_stabilizer(n)
        group = groupify(s)
        unchanged = true
        for stabilizer in group
            apply!(s, stabilizer)
        end 
        @test unchanged == true
    end

    for n in [1,test_sizes...]
        s = random_stabilizer(n)
        group = groupify(s)
        gen_set = get_generating_set(group)
        @test s == gen_set
    end
end
