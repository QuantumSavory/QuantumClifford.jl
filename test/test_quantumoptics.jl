@testitem "Quantum optics" begin
    using QuantumOpticsBase
    using QuantumClifford: @S_str, random_stabilizer
    using LinearAlgebra
    #using QuantumCliffordQOpticsExt: _l0, _l1, _s₊, _s₋, _i₊, _i₋
    const qo = Base.get_extension(QuantumClifford, :QuantumCliffordQOpticsExt)
    const _l0 = qo._l0
    const _l1 = qo._l1
    const _s₊ = qo._s₊
    const _s₋ = qo._s₋
    const _i₊ = qo._i₊
    const _i₋ = qo._i₋

    @testset "conversion from Stabilizer to Ket" begin
        for n in 1:5
            stabs = [random_stabilizer(1) for _ in 1:n]
            stab = tensor(stabs...)
            translate = Dict(S"X"=>_s₊,S"-X"=>_s₋,S"Z"=>_l0,S"-Z"=>_l1,S"Y"=>_i₊,S"-Y"=>_i₋)
            kets = [translate[s] for s in stabs]
            ket = tensor(kets...)
            @test ket.data ≈ Ket(stab).data

            rstab = random_stabilizer(n)
            lstab = random_stabilizer(n)
            lket = Ket(rstab)
            rket = Ket(lstab)
            dotket = abs(lket'*rket)
            dotstab = abs(dot(lstab,rstab))
            @test (dotket<=1e-10 && dotstab==0) || dotket≈dotstab
        end
    end

    @testset "conversion from CliffordOperator to Operator" begin
        for n in 1:3
            for _c in 1:5
                cliff = random_clifford(n)
                U = Operator(cliff)
                for _c in 1:5
                    stab = random_stabilizer(n)
                    ψ₁ = Ket(stab)
                    ψ₂ = Ket(apply!(stab,cliff))
                    # test they are equal up to a phase
                    @test all(x->isnan(x)||abs(x)≈1 , (U*ψ₁).data ./ ψ₂.data)
                    @test abs(det(U.data))≈1
                end
            end
        end
    end

    @testset "conversion from GeneralizedStabilizer to Operator" begin
        for n in 1:5
            stab = random_stabilizer(n)
            @test dm(Ket(stab)) == Operator(GeneralizedStabilizer(stab))
        end
    end

    @testset "conversion from PauliOperator to Operator" begin
        for n in 1:5
            for _ in 1:10
                p = random_pauli(n)
                q = random_pauli(n)
                p̃ = Operator(p)
                q̃ = Operator(q)
                @test Operator(p*q) == p̃*q̃
                @test Operator(q*p) == q̃*p̃
            end
        end
    end


    tgate = sparse(identityoperator(SpinBasis(1//2)))
    tgate.data[2,2] = exp(im*pi/4)

    @testset "GeneralizedStabilizer/PauliChannel to QuantumOptics - explicit single-qubit Pauli channels" begin
        # manual checks
        @test Operator(pcT)≈tgate

        # single qubit checks
        for single_qubit_explicit_channel in [pcT]
            qo_gate = Operator(single_qubit_explicit_channel)
            for single_qubit_tableau in [S"X", S"Y", S"Z", S"-X", S"-Y", S"-Z"]
                sm = GeneralizedStabilizer(single_qubit_tableau)
                ψ = Ket(single_qubit_tableau)
                for rep in 1:8
                    apply!(sm, single_qubit_explicit_channel)
                    ψ = qo_gate*ψ
                    @test expect(Operator(sm), ψ) ≈ 1
                end
            end
        end

        # embedded checks
        for single_qubit_explicit_channel in [pcT]
            for n in 2:5
                i = rand(1:n)
                channel = embed(n,i,single_qubit_explicit_channel)
                qo_gate1 = Operator(single_qubit_explicit_channel)
                qo_gate = embed(basis(qo_gate1)^n, i, qo_gate1)
                stab = random_stabilizer(n)
                sm = GeneralizedStabilizer(stab)
                ψ = Ket(stab)
                for rep in 1:8
                    apply!(sm, channel)
                    ψ = qo_gate*ψ
                    @test expect(Operator(sm), ψ) ≈ 1
                end
            end
        end
    end
end
