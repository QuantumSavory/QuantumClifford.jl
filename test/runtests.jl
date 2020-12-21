using QuantumClifford, Test, Random, Documenter
using AbstractAlgebra
using Quantikz: circuit2string
using QuantumClifford: stab_looks_good, mixed_stab_looks_good, destab_looks_good, mixed_destab_looks_good
using QuantumClifford: CNOTcol, SWAPcol, Hadamardcol, Phasecol, CliffordIdcol
using QuantumClifford.Experimental.NoisyCircuits
using Base.Threads: @threads, nthreads

macro mythreads(arg)
    return arg # TODO bypass this as multithreaded tests do not work yet
    if nthreads()==1
        return arg
    else
        return esc(:( @threads $arg ))
    end
end

test_sizes = [10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

function doset(descr)
    if length(ARGS) == 0
        return true
    end
    for a in ARGS
        if occursin(lowercase(a), lowercase(descr))
            return true
        end
    end
    return false
end

function tests()

Random.seed!(42)

if doset("Doctests") && VERSION >= v"1.6"
@testset "Doctests" begin
    DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)
    doctest(QuantumClifford)
end
end

if doset("Pauli Operators")
@testset "Pauli Operators" begin
    @testset "Parsing, constructors, and properties" begin
        @test P"-iXYZ" == PauliOperator(0x3, 3, vcat(BitArray([1,1,0]).chunks, BitArray([0,1,1]).chunks))
        @test P"-iXYZ" == PauliOperator(0x3, Bool[1,1,0], Bool[0,1,1])
        @test xbit(P"-iXYZ") == Bool[1,1,0]
        @test zbit(P"-iXYZ") == Bool[0,1,1]
        @test P"-iXYZ".xz == UInt64[0x03, 0x06]
        @test P"-iXYZ".phase[] == 0x03 # TODO why is this failing?
        @test P"-iXYZ".nqubits == 3
        @test size(P"-iXYZ") == (3,)
    end
    @testset "Indexing" begin
        @test eachindex(P"IXYZ") == 1:4
        @test P"IXYZ"[3] == (true, true)
        p = P"IXYZ"
        @test p[[3,2]] == P"YX"
        p[4] = (true,false)
        @test p == P"IXYX"
    end
    @testset "Elementary operations" begin
        @test P"X"*P"Z" == P"-iY"
        @test comm(P"XX",P"YY") == 0x0
        @test comm(P"XZ",P"YZ") == 0x1
        @test prodphase(P"XX",P"YY") == 0x2
        @test prodphase(P"ZZZ",P"XXX") == 0x3
    end
    @testset "Commutation implies real phase" begin
        @mythreads for i in 1:10
            @mythreads for n in test_sizes
                p1,p2 = random_pauli(n; nophase=true), random_pauli(n; nophase=true)
                com = comm(p1,p2)==0x0
                p = prodphase(p1,p2)
                rea = p==0x0 || p==0x2
                @test (com && rea) || (!com && !rea)
            end
        end
    end
    @testset "Single qubit Paulis and their action" begin
        @mythreads for i in 1:10
            @mythreads for n in test_sizes
                ix, iy, iz = rand(1:n), rand(1:n), rand(1:n)
                px = single_x(n,ix)
                py = single_y(n,iy)
                pz = single_z(n,iz)
                s1 = random_stabilizer(n)
                s2 = copy(s1)
                apply!(s1,px)
                apply_single_x!(s2,ix)
                @test s1==s2
                apply!(s1,py)
                apply_single_y!(s2,iy)
                @test s1==s2
                apply!(s1,pz)
                apply_single_z!(s2,iz)
                @test s1==s2
            end
        end
    end
end
end

if doset("Pure and Mixed state initialization")
@testset "Pure and Mixed state initialization" begin
    @testset "Destabilizer initialization" begin
        @mythreads for n in test_sizes
            @test destab_looks_good(Destabilizer(random_stabilizer(n)))
        end
    end
    @testset "Mixed destabilizer initialization" begin
        @mythreads for n in test_sizes[2:end]
            @test mixed_destab_looks_good(MixedDestabilizer(random_stabilizer(rand(n÷2+1:n-4),n)))
        end
    end
end
end

if doset("Stabilizer canonicalization")
@testset "Stabilizer canonicalization" begin
    @testset "Default canonicalization" begin
        s = S"- XZZZZ_____
              - _YZY___YX_
              - __XXZ__YX_
              + Z_Z_Y__YXZ
              + _____Z____
              + __________
              + __________
              + ______YYX_
              + __________
              + __________"
        canonicalize!(s)
        t = S"- XZZZZ_____
              - _YZY___YX_
              - __XXZ__YX_
              + Z_Z_Y__YXZ
              + ______YYX_
              + _____Z____
              + __________
              + __________
              + __________
              + __________"
        @test s == t
        for n in test_sizes
            @test stab_looks_good(random_stabilizer(n))
        end
    end
    @testset "Gottesman canonicalization" begin
        @mythreads for n in test_sizes
            @mythreads for nrows in [n, rand(n÷3:n*2÷3)]
                rs = random_stabilizer(nrows,n)
                c = canonicalize!(copy(rs))
                g, _, _, perm1, perm2 = canonicalize_gott!(copy(rs))
                c1 = canonicalize!(colpermute!(colpermute!(copy(rs),perm1),perm2))
                cg = canonicalize!(copy(g))
                @test cg == c1
                @test stab_looks_good(g)
            end
        end
    end
    @testset "Canonicalization of complex tableaus" begin
        @mythreads for n in test_sizes
            @mythreads for nrows in [n, rand(n÷3:n*2÷3)]
                rs = random_stabilizer(nrows,n)
                c  = canonicalize_rref!(copy(rs),1:n)[1]
                dc = canonicalize_rref!(Destabilizer(copy(rs)),1:n)[1]
                mc = canonicalize_rref!(MixedDestabilizer(copy(rs)),1:n)[1]
                @test stabilizerview(mc) == stabilizerview(dc) == c
                @test stab_looks_good(c)
                @test destab_looks_good(dc)
                @test mixed_destab_looks_good(mc)
            end
        end
    end
end
end

if doset("Projective measurements")
@testset "Projective measurements" begin
    @testset "Stabilizer representation" begin
        s = S"XXX
              ZZI
              IZZ"
        ps, anticom, res = project!(copy(s), P"ZII")
        ps = canonicalize!(ps)
        @test anticom==1 && isnothing(res) && ps == S"ZII
                                                      IZI
                                                      IIZ"
        @test stab_looks_good(ps)

        ps, anticom, res = project!(copy(s), P"-XXX")
        @test anticom==0 && res[]==0x2 && ps == canonicalize!(copy(s))
        @test stab_looks_good(ps)

        ps, anticom, res = project!(copy(s), P"-XXX"; keep_result=false)
        @test anticom==0 && isnothing(res) && ps == s
        @test stab_looks_good(ps)

        @mythreads for n in test_sizes
            s = random_stabilizer(n)
            m = random_pauli(n;nophase=true)
            ps, anticom, res = project!(copy(s),m)
            @test anticom==0x0 || ps[anticom]==m
            @test stab_looks_good(ps)
            m = single_z(n,1)
            ps, anticom, res = project!(copy(s),m)
            @test anticom==0x0 || ps[anticom]==m
            @test stab_looks_good(ps)
            m = single_x(n,1)
            ps, anticom, res = project!(copy(s),m)
            @test anticom==0x0 || ps[anticom]==m
            @test stab_looks_good(ps)
        end
    end
    @testset "Destabilizer representation" begin
        @mythreads for n in test_sizes
            s = canonicalize!(random_stabilizer(n))
            m = random_pauli(n;nophase=true)
            ps, anticom, res = project!(copy(s),m)
            dps, danticom, dres = project!(Destabilizer(copy(s)),m)
            @test destab_looks_good(dps)
            @test anticom==danticom && res==dres && canonicalize!(ps)==canonicalize!(stabilizerview(dps))
            m = single_z(n,1)
            ps, anticom, res = project!(copy(s),m)
            dps, danticom, dres = project!(Destabilizer(copy(s)),m)
            @test destab_looks_good(dps)
            @test anticom==danticom && res==dres && canonicalize!(ps)==canonicalize!(stabilizerview(dps))
            m = single_x(n,1)
            ps, anticom, res = project!(copy(s),m)
            dps, danticom, dres = project!(Destabilizer(copy(s)),m)
            @test destab_looks_good(dps)
            @test anticom==danticom && res==dres && canonicalize!(ps)==canonicalize!(stabilizerview(dps))
        end
    end
    @testset "Anticommutation indices and NA results" begin
        s = S" XXX
              -ZZI"
        ds = Destabilizer(copy(s))
        ms = MixedStabilizer(copy(s))
        mds = MixedDestabilizer(copy(s))

        p = P"IZZ"
        ps, a, r = project!(copy(s),p)
        @test stab_looks_good(ps)
        @test a==0 && isnothing(r)
        @test_throws BadDataStructure pds, a, r = project!(copy(ds),p)
        pms, a, r = project!(copy(ms),p)
        @test stab_looks_good(pms)
        @test pms.rank==3
        @test a==0 && isnothing(r)
        pmds, a, r = project!(copy(mds),p)
        @test mixed_destab_looks_good(pmds)
        @test pmds.rank==3
        @test a==0 && isnothing(r)

        p = P"ZZI"
        ps, a, r = project!(copy(s),p)
        @test stab_looks_good(ps)
        @test a==0 && r==0x2
        @test_throws BadDataStructure pds, a, r = project!(copy(ds),p)
        pms, a, r = project!(copy(ms),p)
        @test stab_looks_good(pms)
        @test pms.rank==2
        @test a==0 && r==0x2
        pmds, a, r = project!(copy(mds),p)
        @test mixed_destab_looks_good(pmds)
        @test pmds.rank==2
        @test a==0 && r==0x2
        @test canonicalize!(ps)==canonicalize!(stabilizerview(pms))==canonicalize!(stabilizerview(pmds))

        p = P"XZZ"
        ps, a, r = project!(copy(s),p)
        @test stab_looks_good(ps)
        @test a==2 && isnothing(r)
        pds, a, r = project!(copy(ds),p)
        @test destab_looks_good(pds)
        @test a==2 && isnothing(r)
        pms, a, r = project!(copy(ms),p)
        @test stab_looks_good(pms)
        @test pms.rank==2
        @test a==2 && isnothing(r)
        pmds, a, r = project!(copy(mds),p)
        @test mixed_destab_looks_good(pmds)
        @test pmds.rank==2
        @test a==2 && isnothing(r)
        @test canonicalize!(ps)==canonicalize!(stabilizerview(pms))==canonicalize!(stabilizerview(pds))==canonicalize!(stabilizerview(pmds))
    end
    @testset "Mixed Destabilizer projection on logical operator" begin
        stab = one(MixedDestabilizer, 2,4)
        projzl = single_z(4,1)
        projzr = single_z(4,4)
        projxl = single_x(4,1)
        projxr = single_x(4,4)
        s, a, r = project!(copy(stab), projzl)
        @test mixed_destab_looks_good(s)
        @test a==0 && r==0x0       && stabilizerview(s)==S"Z___
                                                           _Z__"
        s, a, r = project!(copy(stab), projxl)
        @test mixed_destab_looks_good(s)
        @test a==1 && isnothing(r) && stabilizerview(s)==S"X___
                                                           _Z__"
        s, a, r = project!(copy(stab), projzr)
        @test mixed_destab_looks_good(s)
        @test a==0 && isnothing(r) && stabilizerview(s)==S"Z___
                                                           _Z__
                                                           ___Z"
        s, a, r = project!(copy(stab), projxr)
        @test mixed_destab_looks_good(s)
        @test a==0 && isnothing(r) && stabilizerview(s)==S"Z___
                                                           _Z__
                                                           ___X"
    end
end
end

if doset("Partial traces")
@testset "Partial traces" begin
    @testset "RREF canonicalization vs manual traceout" begin
        @mythreads for N in test_sizes
            @mythreads for rep in 1:5
                to_delete = randperm(N)[1:rand(N÷4:N÷2)]
                stab0 = random_stabilizer(N)
                id_paulis = zero(PauliOperator, N)
                # Trace out by doing projective measurements
                naive_stab = copy(stab0)
                for i in to_delete
                    naive_stab, anticom_index, result = project!(naive_stab, single_x(N,i))
                    if anticom_index!=0
                        naive_stab[anticom_index] = id_paulis
                    end
                    naive1_stab, anticom_index, result = project!(naive_stab, single_z(N,i))
                    if anticom_index!=0
                        naive_stab[anticom_index] = id_paulis
                    end
                end
                for i in 1:N
                    for j in to_delete
                        naive_stab[i,j] = (false,false)
                    end
                end
                canonicalize!(naive_stab)
                # Trace out by using the RREF canonical form
                stab = copy(stab0)
                stab, last_row = canonicalize_rref!(stab, to_delete)
                for i in last_row+1:N
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
                s = traceout!(MixedStabilizer(copy(stab0), N), to_delete)
                canonicalize!(s)
                @test stab[1:last_row] == stabilizerview(s)
                @test mixed_stab_looks_good(s)
                # On MixedDestabilizer instances
                s = traceout!(MixedDestabilizer(copy(stab0)), to_delete)
                #@test mixed_destab_looks_good(s) #TODO why is this test failing
                s = canonicalize!(stabilizerview(s))
                @test stab[1:last_row] == s
            end
        end
    end
end
end

if doset("GF(2) representations")
@testset "GF(2) representations" begin
    @testset "Equivalence of GF(2) Gaussian elimination and Stabilizer canonicalization" begin
        @mythreads for n in test_sizes
            @mythreads for rep in 1:5
                s = random_stabilizer(n)[randperm(n)[1:rand(n÷2+1:n)]]
                cs = canonicalize!(copy(s));
                H = stab_to_gf2(cs);
                cH = gf2_gausselim!(stab_to_gf2(s));
                @test H==cH
            end
        end
    end
    @testset "GF(2) H and G matrices" begin
        @mythreads for n in test_sizes
            @mythreads for rep in 1:5
                H = random_invertible_gf2(n)[randperm(n)[1:rand(n÷2+1:n)],:]
                H = gf2_gausselim!(H)
                G = gf2_H_to_G(H)
                @test sum(G*H' .%2)==0;
            end
        end
    end
end
end

if doset("Clifford Operators")
@testset "Clifford Operators" begin
    @testset "Permutations of qubits" begin
        # TODO (see the column version)
    end
    @testset "Tensor products" begin
        # TODO (see the column version)
    end
    @testset "Clifford acting on Stabilizer" begin
        @mythreads for size in test_sizes
            s = random_stabilizer(size)
            gates = vcat([CNOT, Hadamard, Phase], repeat([CliffordId],size-4))
            gates_perm = randperm(size-1)
            gates = gates[gates_perm]
            big_gate = reduce(⊗,gates)

            s1 = apply!(copy(s),big_gate)
            @test stab_looks_good(s1)

            igates_perm = perm_inverse(gates_perm)
            s2 = copy(s)
            canonicalize!(s2)
            s2 = apply!(s2, CNOT, [igates_perm[1],igates_perm[1]+1])
            canonicalize!(s2)
            s2 = apply!(s2, Hadamard, [igates_perm[2]+(igates_perm[1]<igates_perm[2])])
            canonicalize!(s2)
            s2 = apply!(s2, Phase, [igates_perm[3]+(igates_perm[1]<igates_perm[3])])

            @test canonicalize!(s1) == canonicalize!(s2)
        end
    end
end
end

if doset("Clifford Operators (column representation)")
@testset "Clifford Operators (column representation)" begin
    @testset "Permutations of qubits" begin
        function naive_permute(c::CliffordColumnForm,p::AbstractArray{T,1} where T) # TODO this is extremely slow stupid implementation
            ops = QuantumClifford.getallpaulis(c)
            CliffordColumnForm([ops[i][p] for i in 1:2*c.nqubits][vcat(p,p.+c.nqubits)])
        end
        for c in [CNOTcol, CliffordIdcol⊗Hadamardcol, CNOTcol⊗CNOTcol, tensor_pow(CNOTcol,6), tensor_pow(CNOTcol,7), tensor_pow(CNOTcol,6)⊗Phasecol, tensor_pow(CNOTcol,7)⊗Phasecol]
            for rep in 1:5
                p = randperm(nqubits(c))
               @test permute(c,p) == naive_permute(c,p)
            end
        end
        for i in 1:5
            p = randperm(125)
            c = rand([CliffordIdcol, Hadamardcol, Phasecol], 125)
            @test ⊗(c[p]...) == naive_permute(⊗(c...), p) == permute(⊗(c...), p)
        end
    end
    @testset "Tensor products" begin
        function naive_mul(l::CliffordColumnForm, r::CliffordColumnForm) # TODO this is extremely slow stupid implementation
            opsl = QuantumClifford.getallpaulis(l)
            opsr = QuantumClifford.getallpaulis(r)
            onel = zero(opsl[1])
            oner = zero(opsr[1])
            opsl = [l⊗oner for l in opsl]
            opsr = [onel⊗r for r in opsr]
            CliffordColumnForm(vcat(opsl[1:end÷2],opsr[1:end÷2],opsl[end÷2+1:end],opsr[end÷2+1:end]))
        end
        function naive_tensor_pow(op::CliffordColumnForm,power::Integer,mem::Dict{Integer,CliffordColumnForm})
            if power==1
                return op
            elseif haskey(mem,power)
                return mem[power]
            end
            half,rest = divrem(power,2)
            phalf = get!(mem,half) do
                naive_tensor_pow(op,half,mem)
            end
            res = naive_mul(phalf,phalf)
            if rest!=0
                prest = get!(mem,rest) do
                    naive_tensor_pow(op,half,mem)
                end
                res = naive_mul(res,prest)
            end
            res
        end
        function naive_tensor_pow(op::CliffordColumnForm,power::Integer)
            naive_tensor_pow(op,power,Dict{Integer,CliffordColumnForm}())
        end
        for s in test_sizes
            for g in [CNOTcol,Hadamardcol,CNOTcol⊗Phasecol]
                @test naive_tensor_pow(g,s)==tensor_pow(g,s)
            end
        end
        for g1 in [CNOTcol,Hadamardcol,CNOTcol⊗Phasecol,CliffordIdcol]
            for g2 in [CNOTcol,Hadamardcol,CNOTcol⊗Phasecol,CliffordIdcol]
                @test naive_mul(g1,g2)==g1⊗g2
            end
        end
        c1 = tensor_pow(CNOTcol,32)
        @test naive_mul(c1,c1) == c1⊗c1
        c1 = naive_tensor_pow(Hadamardcol,33)
        @test naive_mul(c1,c1) == c1⊗c1
        c2 = naive_tensor_pow(Hadamardcol,32)
        @test naive_mul(c1,c2) == c1⊗c2
        @test naive_mul(c2,c1) == c2⊗c1
    end
    @testset "Clifford acting on Stabilizer" begin
        @mythreads for size in test_sizes
            s = random_stabilizer(size)
            gates = vcat([CNOTcol, Hadamardcol, Phasecol], repeat([CliffordIdcol],size-4))
            gates_perm = randperm(size-1)
            gates = gates[gates_perm]
            big_gate = reduce(⊗,gates)

            s1 = apply!(copy(s),big_gate; phases=false)
            @test stab_looks_good(s1)

            igates_perm = perm_inverse(gates_perm)
            s2 = copy(s)
            s2 = apply!(s2, CNOTcol, [igates_perm[1],igates_perm[1]+1]; phases=false)
            s2 = apply!(s2, Hadamardcol, [igates_perm[2]+(igates_perm[1]<igates_perm[2])]; phases=false)
            s2 = apply!(s2, Phasecol, [igates_perm[3]+(igates_perm[1]<igates_perm[3])]; phases=false)

            @test s1 == s2
        end
    end
end
end

if doset("Bug fixes - regression tests")
@testset "Bug fixes - regression tests" begin
    @testset "Redundant row permutations in `project!(::MixedDestabilizer)`" begin
        # Fixed in 41ed1d3c
        destab =  S"+ ZX_Y_YXZ
                    + XY_Y____
                    + _Z_XXY__
                    + _ZYXXY__
                    + X__Y_ZXZ
                    + X__YXZXZ
                    + ___YXXZZ
                    + _______Z"
        stab =    S"+ X_______
                    + _X_Y____
                    + __ZY____
                    + __Z_____
                    + ___YZY__
                    + X__YZYZZ
                    + X____YZZ
                    + ______YX"
        t = MixedDestabilizer(vcat(destab,stab), 8)
        @test mixed_destab_looks_good(t)
        c = copy(stabilizerview(t)[[1,3,5,7]])
        traceout!(t,[1,4,3,6])
        @test mixed_destab_looks_good(t)
        project!(t,c[1])
        @test mixed_destab_looks_good(t)
        project!(t,c[2])
        @test mixed_destab_looks_good(t) # This used to fail because anticomlog==rank+1 leading to a repeated row permutation
    end
end
end

end

function noisycircuits_tests()
if doset("Noisy Circuits")
@testset "Noisy Circuits" begin

if doset("Monte Carlo sims")
@testset "Monte Carlo sims" begin
    @testset "Purification examples" begin
        g1 = SparseGate(CNOT, [1,3])
        g2 = SparseGate(CNOT, [2,4])
        m = BellMeasurement([X,X],[3,4])
        good_bell_state = S"XX
                            ZZ"
        canonicalize_rref!(good_bell_state)
        v = VerifyOp(good_bell_state,[1,2])
        n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.01))
        with_purification = mctrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v], trajectories=500)
        @test with_purification[:detected_failure] > 5
        @test with_purification[:undetected_failure] > 10
        @test with_purification[:true_success] > 430
        without_purification = mctrajectories(good_bell_state⊗good_bell_state, [n,v], trajectories=500)
        @test without_purification[:detected_failure] == 0
        @test without_purification[:undetected_failure] > 10
        @test without_purification[:true_success] > 450
        nonoise = mctrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v], trajectories=10)
        @test nonoise[:detected_failure] == 0
        @test nonoise[:undetected_failure] == 0
        @test nonoise[:true_success] == 10
    end
end
end

if doset("Perturbative expansion sims")
@testset "Perturbative expansion sims" begin
    @testset "Purification examples comparison to MC" begin
        compare(a,b, symbol) = abs(a[symbol]/500-b[symbol]) / (a[symbol]/500+b[symbol]+1e-5) < 0.3
        g1 = SparseGate(CNOT, [1,3])
        g2 = SparseGate(CNOT, [2,4])
        m = BellMeasurement([X,X],[3,4])
        good_bell_state = S"XX
                            ZZ"
        canonicalize_rref!(good_bell_state)
        v = VerifyOp(good_bell_state,[1,2])
        n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.01))
        mc = mctrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v], trajectories=500)
        pe = petrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v])
        @test compare(mc,pe,:detected_failure)
        @test compare(mc,pe,:undetected_failure)
        @test compare(mc,pe,:true_success)
        mc = mctrajectories(good_bell_state⊗good_bell_state, [n,v], trajectories=500)
        pe = petrajectories(good_bell_state⊗good_bell_state, [n,v])
        @test compare(mc,pe,:detected_failure)
        @test compare(mc,pe,:undetected_failure)
        @test compare(mc,pe,:true_success)
        mc = mctrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v], trajectories=500)
        pe = petrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v])
        @test compare(mc,pe,:detected_failure)
        @test compare(mc,pe,:undetected_failure)
        @test compare(mc,pe,:true_success)
    end
    
    @testset "Symbolic" begin
        R, (e,) = PolynomialRing(RealField, ["e"])
        unity = R(1);
        
        good_bell_state = S"XX
                            ZZ"
        canonicalize_rref!(good_bell_state)
        initial_state = good_bell_state⊗good_bell_state

        g1 = SparseGate(CNOT, [1,3]) # CNOT between qubit 1 and qubit 3 (both with Alice)
        g2 = SparseGate(CNOT, [2,4]) # CNOT between qubit 2 and qubit 4 (both with Bob)
        m = BellMeasurement([X,X],[3,4]) # Bell measurement on qubit 3 and 4
        v = VerifyOp(good_bell_state,[1,2]) # Verify that qubit 1 and 2 indeed form a good Bell pair
        epsilon = e # The error rate
        n = NoiseOpAll(UnbiasedUncorrelatedNoise(epsilon))

        # This circuit performs a depolarization at rate `epsilon` to all qubits,
        # then bilater CNOT operations
        # then a Bell measurement
        # followed by checking whether the final result indeed corresponds to the correct Bell pair.
        circuit = [n,g1,g2,m,v]

        pe_symbolic = petrajectories(initial_state, circuit, branch_weight=unity) # perturbative expansion
        @test pe_symbolic[:undetected_failure] == -162.0*e^4 + 162.0*e^3 + -54.0*e^2 + 6.0*e
        @test pe_symbolic[:detected_failure]   == -108.0*e^4 + 108.0*e^3 + -36.0*e^2 + 4.0*e
        @test pe_symbolic[:true_success]       == 27.0*e^4 + -54.0*e^3 + 36.0*e^2 + -10.0*e + 1.0
    end
end
end

if doset("Quantikz diagrams")
@testset "Quantikz diagrams" begin
    noise = UnbiasedUncorrelatedNoise(0.1)
    @test circuit2string([
            SparseGate(CNOT, [1,4]),
            SparseGate(CNOT, [3,2]),
            SparseGate(CPHASE, [1,2]),
            SparseGate(SWAP, [2,4]),
            SparseGate(CNOT*CNOT, [1,3]),
            NoiseOp(noise,[1,3]),
            NoiseOpAll(noise),
            NoisyGate(SparseGate(CNOT*CNOT, [2,4]),noise),
            
            ]) == "\\begin{quantikz}[transparent]\n\\qw & \\ctrl{0} & \\qw & \\ctrl{0} & \\qw & \\gate[3,label style={yshift=0.2cm}]{} & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\qw & \\qw\\\\\n\\qw & \\qw & \\targ{}\\vqw{0} & \\ctrl{-1} & \\swap{0} & \\linethrough & \\qw & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\gate[3,label style={yshift=0.2cm}]{} & \\qw\\\\\n\\qw & \\qw & \\ctrl{-1} & \\qw & \\qw &  & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\linethrough & \\qw\\\\\n\\qw & \\targ{}\\vqw{-3} & \\qw & \\qw & \\swap{-2} & \\qw & \\qw & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} &  & \\qw\n\\end{quantikz}"
end
end

end
end

end

tests()
noisycircuits_tests()