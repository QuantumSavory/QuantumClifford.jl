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
        # todo test MRZ and other random gates statistically
        circuit = [sHadamard(2), sHadamard(5), sCNOT(1, 2), sCNOT(2, 5), sMZ(1), sMZ(2), sMZ(4), sMZ(5)]
        ccircuit = if eltype(circuit) <: QuantumClifford.CompactifiedGate
            circuit
        else
            compactify_circuit(circuit)
        end

        CUDA.allowscalar(false) # todo how to add this to all tests.

        frames = QuantumClifford._create_pauliframe(ccircuit; trajectories=10)
        cpu_frames = to_cpu(frames)
        gpu_frames = to_gpu(frames)
        cpu_result = pftrajectories(cpu_frames, ccircuit)
        gpu_result = pftrajectories(gpu_frames, ccircuit)
        (cpu_result.frame == to_cpu(gpu_result.frame)) && (cpu_result.measurements == to_cpu(gpu_result.measurements))
    end
end
