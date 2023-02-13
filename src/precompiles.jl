function _precompile_()
    ds = random_destabilizer(3);
    s = random_stabilizer(3);
    canonicalize!(s);
    canonicalize_rref!(s);
    canonicalize_gott!(s);
    canonicalize!(s,phases=false)
    canonicalize_gott!(s,phases=false)
    canonicalize_rref!(s,phases=false)
    c = random_clifford(3)
    apply!(s, c);
    apply!(s, tCNOT, [1,2]);
    apply!(s, sCNOT(1,2));
    project!(s, P"XXX");
    s = S"XX
          ZZ"
    md = MixedDestabilizer(s)
    p = P"XX"
    c = C"X_
          _X
          Z_
          _Z"
    p*md
    c*md
    apply!(md,p,phases=false)
    apply!(md,c,phases=false)
    random_clifford(3)
    random_pauli(3)
    for op in [ sHadamard
                sId1
                sInvPhase
                sPhase
                sX
                sY
                sZ
            ]#subtypes(QuantumClifford.AbstractSingleQubitOperator)
        op(2)*s
    end
    for op in [ sCNOT
                sSWAP
                sCPHASE
            ]#subtypes(QuantumClifford.AbstractTwoQubitOperator)
        op(2,1)*s
    end
    project!(s, p)
    project!(md,p)
end

using SnoopPrecompile

# precompilation causes allocation performance bugs for <v1.8 https://github.com/JuliaLang/julia/issues/35972
VERSION > v"1.8" && @precompile_setup begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        _precompile_()
    end
end
