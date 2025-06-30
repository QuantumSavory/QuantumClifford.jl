@testitem "GPU" tags=[:cuda] begin
    using QuantumClifford: to_cpu, to_gpu
    using CUDA

    function apply_single_qubit_and_compare(n, s, s_gpu)
        qbitidx = (((rand(Int) % n) + n) %n) + 1
        op = random_clifford1(qbitidx)
        apply!(s_gpu, op)
        apply!(s, op)
        if to_cpu(s_gpu) != s
            throw("result of cpu and gpu differ in single qubit operation. n=$n iteration=$iteration operation=$op")
        end
    end

    function apply_single_or_double_qubit_and_compare(n, s, s_gpu)
        op = if rand(Int) % 2 == 0
            qbitidx = (((rand(Int) % n) + n) %n) + 1
            random_clifford1(qbitidx)
        else
            # todo what other gates can we use? how to randomly generate two qubit operation?
            # can we use random_clifford(2) with this function?
            qbitidx1 = (((rand(Int) % n) + n) %n) + 1
            qbitidx2 = (((rand(Int) % n) + n) %n) + 1
            if qbitidx1 == qbitidx2
                qbitidx2 = (qbitidx1 % n) + 1
            end
            sCNOT(qbitidx1, qbitidx2)
        end
        apply!(s_gpu, op)
        apply!(s, op)
        if to_cpu(s_gpu) != s
            throw("result of cpu and gpu differ in single qubit operation. n=$n iteration=$iteration operation=$op")
        end
    end

    function approx(a, b, threshold)
        maximum(abs.(a - b)) <= threshold
    end


    @testset "GPU" begin
        CUDA.allowscalar(false) # make sure we are using GPU kernels and not iterating on indices

        @test begin
            p = random_pauli(3)
            p_gpu = to_gpu(p)
            typeof(p_gpu.xz) <: CUDA.CuArray # todo this is a bad test because it depends on representation of data in QuantumClifford. Change later...
        end
        @test begin
            s = random_stabilizer(3)
            s_gpu = to_gpu(s)
            typeof(tab(s_gpu).xzs) <: CUDA.CuArray # todo this is a bad test because it depends on representation of data in QuantumClifford. Change later...
        end

        @test begin
            for n in [2, 4, 8, 100, 500]
                s = random_stabilizer(n)
                s_gpu = to_gpu(s)
                for iteration in 1:100
                    apply_single_qubit_and_compare(n, s, s_gpu)
                end
            end
            true
        end

        @test begin
            for n in [2, 4, 8, 100, 500]
                s = random_stabilizer(n)
                s_gpu = to_gpu(s)
                for iteration in 1:100
                    apply_single_or_double_qubit_and_compare(n, s, s_gpu)
                end
            end
            true
        end

        @test begin
            # todo test MRZ and other random gates statistically
            circuit = [sHadamard(2), sHadamard(5), sCNOT(1, 2), sCNOT(2, 5), sMZ(1), sMZ(2), sMZ(4), sMZ(5)]
            ccircuit = if eltype(circuit) <: QuantumClifford.CompactifiedGate
                circuit
            else
                compactify_circuit(circuit)
            end

            for func in [identity, fastcolumn, fastrow]
                frames = QuantumClifford._create_pauliframe(ccircuit; trajectories=10)
                cpu_frames = func(to_cpu(frames))
                gpu_frames = func(to_gpu(frames))
                cpu_result = pftrajectories(cpu_frames, ccircuit)
                gpu_result = pftrajectories(gpu_frames, ccircuit)
                check = (cpu_result.frame == to_cpu(gpu_result.frame)) && (cpu_result.measurements == to_cpu(gpu_result.measurements))
                if !check
                    throw("pftrajectories produce wrong result after applying $func")
                end
            end
            true
        end

        @test begin
            # test fastrow
            for n in [2, 4, 8, 100, 500]
                s = fastrow(random_stabilizer(n))
                s_gpu = fastrow(to_gpu(s))
                for iteration in 1:100
                    # how to check if after each iteration we are still fastrow
                    apply_single_or_double_qubit_and_compare(n, s, s_gpu)
                end
            end
            true
        end

        @test begin
            # test fastcolumn
            for n in [2, 4, 8, 100, 500]
                s = fastcolumn(random_stabilizer(n))
                s_gpu = fastcolumn(to_gpu(s))
                for iteration in 1:100
                    # how to check if after each iteration we are still fastrow
                    apply_single_or_double_qubit_and_compare(n, s, s_gpu)
                end
            end
            true
        end

        @test begin
            # test applynoise
            N = 4
            trajectories = 10000
            f = PauliFrame(trajectories, N, N);
            f = to_gpu(f)
            p = 0.1
            measurements = pftrajectories(f, [
                                              sMZ(1, 1),
                                              sHadamard(2),
                                              sMZ(2, 2),
                                              NoiseOp(UnbiasedUncorrelatedNoise(p), [3, 4]),
                                              sMZ(3, 3),
                                              sHadamard(4),
                                              sMZ(4, 4)
                                             ]).measurements
            avg_result = to_cpu(sum(measurements, dims=1) / trajectories)
            error_threshold = 0.02
            approx(vec(avg_result), [0, .5, 2p/3, .5], error_threshold)
        end
    end
end
