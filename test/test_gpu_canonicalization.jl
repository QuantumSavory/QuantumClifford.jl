@testitem "GPU Canonicalization" tags=[:cuda] begin
    using CUDA
    using QuantumClifford
    using Random
    using QuantumClifford: to_cpu, to_gpu, random_tableau
    
    if CUDA.functional()
        @testset "GPU canonicalize! correctness" begin
            for n in [10, 50, 100]
                cpu_stab = random_stabilizer(n)
                gpu_stab = to_gpu(cpu_stab)
                cpu_result = canonicalize!(copy(cpu_stab))
                gpu_result = canonicalize!(copy(gpu_stab))
                cpu_from_gpu = to_cpu(gpu_result)
                @test cpu_result == cpu_from_gpu
            end
        end
        
        @testset "GPU canonicalize! performance" begin
            n = 6000  # Large enough to see GPU benefit
            # Only linear independence is required to achieve full rank.
            # Occurs with unit probability for arbitrarily large qubit counts.
            cpu_stab = Stabilizer(random_tableau(n, n))
            gpu_stab = to_gpu(cpu_stab)
            canonicalize!(copy(gpu_stab))
            gpu_time = @elapsed canonicalize!(copy(gpu_stab))
            # Sanity check, note that for 6000 size Stabilizer cpu version takes 1900-2400 ms on average
            @test 1000 * gpu_time < 1900
        end
        
        @testset "GPU canonicalize! phases option" begin
            n = 50
            cpu_stab = random_stabilizer(n)
            gpu_stab = to_gpu(cpu_stab)
            gpu_result_no_phases = canonicalize!(copy(gpu_stab); phases=false)
            cpu_result_no_phases = canonicalize!(copy(cpu_stab); phases=false)
            
            @test to_cpu(gpu_result_no_phases) == cpu_result_no_phases
        end
    else
        @test_skip "CUDA not functional, skipping GPU canonicalization tests"
    end
end
