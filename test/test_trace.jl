using Random
using QuantumClifford

using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

@testset "Partial traces" begin
    @testset "RREF canonicalization vs manual traceout" begin
        for N in test_sizes
            for n in [N,rand(N÷4:N÷2)]
                to_delete = randperm(N)[1:rand(N÷4:N÷2)]
                stab0 = random_stabilizer(n, N)
                id_paulis = zero(PauliOperator, N)
                # Trace out by doing projective measurements
                naive_stab = copy(stab0)
                for i in to_delete
                    naive_stab, anticom_index, result = project!(naive_stab, single_x(N,i))
                    if anticom_index!=0 && anticom_index<=length(naive_stab)
                        naive_stab[anticom_index] = id_paulis
                    end
                    naive_stab, anticom_index, result = project!(naive_stab, single_z(N,i))
                    if anticom_index!=0 && anticom_index<=length(naive_stab)
                        naive_stab[anticom_index] = id_paulis
                    end
                end
                canonicalize!(naive_stab)
                # Trace out by using the RREF canonical form
                stab = copy(stab0)
                stab, last_row = canonicalize_rref!(stab, to_delete)
                for i in last_row+1:n
                    stab[i] = id_paulis
                end
                canonicalize!(stab)
                # Confirm the results are the same
                @test stab == naive_stab
                @test mixed_stab_looks_good(stab[1:last_row])
                # Check the built-in traceout! functions for this
                s = traceout!(copy(stab0), to_delete)
                canonicalize!(s)
                @test stab == s
                # On MixedStabilizer instances
                s = traceout!(MixedStabilizer(copy(stab0), n), to_delete)
                canonicalize!(s)
                @test stab[1:last_row] == stabilizerview(s)
                @test mixed_stab_looks_good(s)
                # On MixedDestabilizer instances
                s = traceout!(MixedDestabilizer(copy(stab0)), to_delete)
                @test mixed_destab_looks_good(s)
                s = canonicalize!(stabilizerview(s))
                @test stab[1:last_row] == s
            end
        end
    end
end
@testset "Qubit resets" begin
    @test_throws DimensionMismatch reset_qubits!(S"XX YY",S"X",[1,2])
    @test_throws DimensionMismatch reset_qubits!(S"X",S"XX YY",[1,2])
    @test_throws DimensionMismatch reset_qubits!(MixedStabilizer(S"XX YY"),S"X",[1,2])
    @test_throws DimensionMismatch reset_qubits!(MixedStabilizer(S"X"),S"XX YY",[1,2])
    @test_throws DimensionMismatch reset_qubits!(MixedDestabilizer(S"XX YY"),S"X",[1,2])
    @test_throws DimensionMismatch reset_qubits!(MixedDestabilizer(S"X"),S"XX YY",[1,2])
    for N in test_sizes
        for R in [rand(N÷2:N*2÷3), N]
            R > 0 || continue
            s = random_stabilizer(R,N)
            nnew = rand(N÷4:N*2÷3)
            nnew > 0 || continue
            newstate = random_stabilizer(nnew)
            perm = randperm(N)[1:nqubits(newstate)]
            to_trace = setdiff(1:N,perm)
            resetop = Reset(newstate, perm)
            # Testing MixedDestabilizer
            md = MixedDestabilizer(s)
            mdr1 = reset_qubits!(copy(md), newstate,perm)
            @test mixed_destab_looks_good(mdr1)
            mdr2 = reset_qubits!(copy(mdr1),newstate,perm)
            mdro = apply!(copy(md), resetop)
            @test canonicalize!(copy(stabilizerview(mdr1)))==canonicalize!(copy(stabilizerview(mdr2)))==canonicalize!(copy(stabilizerview(mdro)))
            traceout!(mdr2,to_trace)
            mdr2v = stabilizerview(mdr2)
            @test canonicalize!(copy(mdr2v)[:,perm]) == canonicalize!(copy(newstate))
            # Testing MixedStabilizer
            ms = MixedStabilizer(s)
            msr1 = reset_qubits!(copy(ms), newstate,perm)
            @test mixed_stab_looks_good(msr1)
            msr2 = reset_qubits!(copy(msr1),newstate,perm)
            msro = apply!(copy(ms), resetop)
            @test msr1==msr2==msro
            traceout!(msr2,to_trace)
            msr2v = stabilizerview(msr2)
            @test canonicalize!(copy(msr2v)[:,perm]) == canonicalize!(copy(newstate))
            @test canonicalize!(msr2v) == canonicalize!(mdr2v)
            # Testing Stabilizer
            ss = R==N ? s : Stabilizer(tab(MixedStabilizer(s))) # Ensure the tableau is padded with Is
            ssr1 = reset_qubits!(copy(ss), newstate,perm)
            ssr2 = reset_qubits!(copy(ssr1),newstate,perm)
            ssro = apply!(copy(ss), resetop)
            @test canonicalize!(ssr1)==canonicalize!(ssr2)==canonicalize!(ssro)
            traceout!(ssr2,to_trace)
            ssr2v = stabilizerview(ssr2)
            c, x, z = canonicalize!(ssr2v, ranks=true)
            @test canonicalize!(copy(ssr2v)[:,perm])[1:z] == canonicalize!(copy(newstate))
            @test canonicalize!(msr2v) == c[1:z]
            # Compare different datastractures
            @test canonicalize!(copy(stabilizerview(mdr1)))==canonicalize!(copy(stabilizerview(msr1)))==canonicalize!(ssr1[1:mdr1.rank])
        end
    end
end
