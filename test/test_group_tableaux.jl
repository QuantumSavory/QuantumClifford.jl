using Random
using QuantumClifford

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

#testing operations
#= for n in 1:10
    s = random_stabilizer(n)
    group = groupify(s)
    unchanged = true
    for stabilizer in group
        s1 = copy(s)
        apply!(s1, stabilizer)
        if (s == s1)
            print("a")
        end
    end 
    #@test unchanged == true
end
stab = S"XZ_
        +ZXZ
        +_ZX"

println(stab)
g, h_idx, ip_idx, z_idx = graphstate(stab);
plot = graphplot(g)
gates = graph_gatesequence(h_idx, ip_idx, z_idx)
println(gates)
println("")
gate = graph_gate(h_idx, ip_idx, z_idx, nqubits(stab))
println(gate)
println(typeof(gate))
for ga in vcat(gates...)
    print("a") 
    println(ga) 
    apply!(stab, ga) 
end
println(stab)
println("a")
#println(groupify(stab))
println("a")
print(sZ(1))
apply!(stab, sZ(1))
apply!(stab, sZ(2))
apply!(stab, sZ(3))
apply!(stab, stab[1])
println(stab) =#
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
