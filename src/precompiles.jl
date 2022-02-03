function _my_precompile_()
    s = random_stabilizer(3);
    canonicalize!(s);
    canonicalize_rref!(s);
    canonicalize_gott!(s);
    canonicalize!(s,phases=false)
    canonicalize_gott!(s,phases=false)
    canonicalize_rref!(s,phases=false)
    c = random_clifford(3)
    apply!(s, c);
    apply!(s, CNOT, [1,2]);
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
            ]#subtypes(QuantumClifford.AbstractTwoQubitOperator)
        op(2,1)*s
    end
    project!(s, p)
    project!(md,p)
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    _my_precompile_()
end

_precompile_()