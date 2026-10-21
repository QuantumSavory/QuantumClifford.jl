@testitem "ECC - quantum convolutional codes" tags=[:ecc] begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: decode, memory, frame_size, polynomial_form, shifted_commutes
    using Random

    # rate-1/3 code of Forney, Grassl & Guha (2007), Example 1 of Wilde, Krovi & Brun
    c = QuantumConvolutionalCode([[P"XXX", P"XZY"], [P"ZZZ", P"ZYX"]])
    pauli(g, n) = PauliOperator(0x0, g[1:n], g[n+1:end])
    w(P) = count(q -> P[q] != (false, false), 1:nqubits(P))

    @testset "code properties" begin
        @test frame_size(c) == 3 && memory(c) == 1 && rate(c) == 1 // 3
        @test polynomial_form(c) == ["XXX + D·XZY", "ZZZ + D·ZYX"]
        @test shifted_commutes(c)
        bad = QuantumConvolutionalCode([[P"XI", P"ZI"]])     # anticommutes with its own shift
        @test !shifted_commutes(bad)
        @test_throws ArgumentError TailBiting(bad, 4)
        @test_throws ArgumentError TailBiting(c, 1)
        @test_throws ArgumentError QuantumConvolutionalCode([[P"XX", P"X"]])
    end

    @testset "tail-biting code, N = $N" for N in (3, 4, 6, 8)
        tb = TailBiting(c, N)
        H = parity_checks(tb)
        @test code_n(tb) == 3N && code_k(tb) == N && code_s(tb) == 2N
        @test all(comm(H[i], H[j]) == 0 for i in 1:2N, j in 1:2N)
        @test size(faults_matrix(tb)) == (2N, 6N)
        # shift invariance: moving every qubit by one frame permutes the checks
        shift(P) = PauliOperator(0x0, circshift(xbit(P), 3), circshift(zbit(P), 3))
        @test sort(string.(shift.(collect(H)))) == sort(string.(collect(H)))
    end

    @testset "Viterbi decoder explains the syndrome with a lightest error" begin
        Random.seed!(42)
        for N in (4, 6)
            tb = TailBiting(c, N); d = ConvolutionalViterbiDecoder(tb); H = parity_checks(tb)
            for _ in 1:100
                # random errors of weight ≤ 1: decoded exactly up to stabilizers
                E = zero(PauliOperator, 3N)
                rand() < 0.8 && (E[rand(1:3N)] = rand(((true, false), (false, true), (true, true))))
                s = comm(E, H)
                G = pauli(decode(d, s), 3N)
                @test comm(G, H) == s
                @test w(G) <= w(E)
            end
        end
        # every single error, on every qubit, is corrected exactly (distance 3)
        tb = TailBiting(c, 5); d = ConvolutionalViterbiDecoder(tb); H = parity_checks(tb)
        for q in 1:15, P in ((true, false), (false, true), (true, true))
            E = zero(PauliOperator, 15); E[q] = P
            @test pauli(decode(d, comm(E, H)), 15) == E
        end
    end

    @testset "memory padded with an identity frame decodes the same" begin
        c2 = QuantumConvolutionalCode([[P"XXX", P"XZY", P"III"], [P"ZZZ", P"ZYX", P"III"]])
        @test memory(c2) == 2
        tb1, tb2 = TailBiting(c, 5), TailBiting(c2, 5)
        @test parity_checks(tb1) == parity_checks(tb2)
        d1, d2 = ConvolutionalViterbiDecoder(tb1), ConvolutionalViterbiDecoder(tb2)
        Random.seed!(1)
        for _ in 1:3
            s = rand(Bool, 10)
            @test w(pauli(decode(d1, s), 15)) == w(pauli(decode(d2, s), 15))
        end
    end

    @testset "evaluate_decoder: Viterbi beats the weight-1 lookup table" begin
        tb = TailBiting(c, 6)
        Random.seed!(7)
        pv = evaluate_decoder(ConvolutionalViterbiDecoder(tb), CommutationCheckECCSetup(0.02), 5000)
        pt = evaluate_decoder(TableDecoder(tb), CommutationCheckECCSetup(0.02), 5000)
        @test pv < pt
    end
end
