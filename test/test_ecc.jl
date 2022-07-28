
#SEPARATE IMPLMENTATION

#Syndrome circuit testing-------------------------------------
function test_naive_syndrome_circuit()
    @testset "Naive circuits" begin
        circuit = naive_syndrome_circuit(s,i) #TO DO: other syndrome circuits
        test_state = random_stabilizer(n)
        results_direct = project!(copy(test_state), s[i])
end
#--------------------------------------------------------------

#Enconding circuit testing-------------------------------------



#--------------------------------------------------------------