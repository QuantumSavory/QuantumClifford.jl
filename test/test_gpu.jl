using Test
using QuantumClifford
using CUDA

@testset "GPU" begin
    @test begin
        p = P"_IZXY"
        p_gpu = to_gpu(p)
        typeof(p_gpu.xz) <: CUDA.CuArray # todo this is a bad test because it depends on representation of data in QuantumClifford. Change later...
    end
    @test begin
        s = S"-XX
              +ZZ"
        s_gpu = to_gpu(s)
        typeof(tab(s_gpu).xzs) <: CUDA.CuArray # todo this is a bad test because it depends on representation of data in QuantumClifford. Change later...
    end

    @test begin
        s = to_gpu(S"XX XZ")
        op = SingleQubitOperator(sHadamard(1))
        apply!(s, op) # todo rewrite this with normal apply not _apply
        correct = S"ZX ZZ"
        to_cpu(s) == correct
    end

    @test begin
        s = random_stabilizer(10, 10)
        s_gpu = to_gpu(s)
        op = SingleQubitOperator(sHadamard(5))
        apply!(s_gpu, op) # todo rewrite this with normal apply not _apply
        apply!(s, op)
        to_cpu(s_gpu) == s
    end

    @test begin
        s = random_stabilizer(10, 10)
        s_gpu = to_gpu(s)
        op = sCNOT(2, 3)
        apply!(s_gpu, op) # todo rewrite this with normal apply not _apply 
        apply!(s, op)
        to_cpu(s_gpu) == s
    end

    @test begin
        # todo add this tests when we add GPU flag to pftrajectories function
        # circuite = [sHadamard(2), sHadamard(5), sCNOT(1, 2), sCNOT(2, 5), sMRZ(1), sMRZ(2), sMZ(4), sMZ(5)]
        # QuantumCliffordGPU.pftrajectories(circuite; trajectories=10_000)
        true
    end
end
