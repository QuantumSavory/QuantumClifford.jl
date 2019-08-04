# Pauli Operators

They are stored in memory as a phase (a single byte where `0x0,0x1,0x2,0x3` corresponds to $1,i,-1,-i$) and two bit-arrays, for X and for Z components.

You can create them with a `P` string.

```jldoctest
julia> P"-iXZ"
```

Or by specifying phase and X/Z components:

```jldoctest
julia> PauliOperator(0x0,Bool[0,1,0],Bool[0,0,1])
```

Both underscore and I can be used for identity.

```jldoctest
julia> P"I_XYZ"
```

Multiplication with scalars or other Pauli operators works as expected,
as well as tensor products of Pauli operators.

```jldoctest
julia> -1im*P"X"
julia> P"X" * P"Z"
julia> P"X" ⊗ P"Z"
```

One can check for commutativity with `comm`.

```jldoctest
julia> comm(P"X",P"Z")
julia> comm(P"XX",P"ZZ")
```

And check the phase of a product with `prodphase`.

```jldoctest
julia> prodphase(P"X", P"Z")
julia> prodphase(P"X", P"iZ")
julia> prodphase(P"X",P"Y")
```

Indexing operations are available.

```jldoctest
julia> p = P"IXZY";
julia> p[1], p[2], p[3], p[4]
julia> p = P"III";
julia> p[2] = (true, true);
julia> p
```
Including various permutation methods and fancy indexing:

```jldoctest
julia> permute(P"XYZ", [3,1,2])
julia> P"IXYZ"[[2,3]]
```

The operator is represented in memory by bit arrays (much denser than using byte
arrays).

```jldoctest
julia> p = P"-IXYZ";
julia> p.nqbits, p.phase, p.xz
```

The convenience properties `xbit` and `zbit` give you Bool (GF2) vectors.
TODO: this should be a separate function.

```jldoctest
julia> P"XYZI".xbit
```

# Stabilizers

They are stored in memory as a phase list and a bit-matrix for X and Z
components. They can be created by an `S` string, or with a number of different
constructors.

```jldoctest
julia> S"-XX
         +ZZ"
julia> Stabilizer([P"-XX",P"+ZZ"])
julia> Stabilizer([0x2, 0x0],
    Bool[1 1;
         0 0],
    Bool[0 0;
         1 1])
julia> Stabilizer([0x2, 0x0],
    Bool[1 1;
         0 0;
         0 0;
         1 1])
```


Direct sums can be performed,

```jldoctest
julia> S"-XX" ⊕ S"ZZ"
```

Indexing operations are available, including fancy indexing. 2D indexing,
into the Pauli operators is not implemented (TODO).

```jldoctest
julia> s = S"-XXX
             +ZZX
             +III";
julia> s[2]
julia> s[1:2]
julia> s[[1,3]]
julia> s[[true, false, true]]
julia> s[1] = P"+YYX";
julia> s
```

Consistency at creation is not verified so nonsensical stabilizers can be
created, both in terms of content and shape.

```jldoctest
julia> S"iX
         +Z"
```

Similarly to the Pauli operators, a bit array representation is used.

```jldoctest stab
# Representation in memory
julia> s = S"-XXX
             +ZZI
             -IZZ";
julia> s.phases, s.nqbits, s.xzs
```

And there are convenience functions that can extract the corresponding binary
check matrix.

```jldoctest stab
julia> stab_to_gf2(s)
```

# Canonicalization of Stabilizers

Canonicalization (akin to Gaussian elimination over F(2,2)) is implemented.

```jldoctest
julia> s = S"-XXX
             +ZZX
             +III";
julia> canonicalize!(s)
```

If phases are inconsequential, the operations can be faster by not tracking and updating them.

```jldoctest
julia> s = S"-XXX
             +ZZX
             +III";
julia> canonicalize!(s; phases=false)
```

These operations are in place (as customarily signified by "!").

```jldoctest
julia> s = S"-XXX
             +ZZX
             +III";
julia> canonicalize!(s; phases=false);
julia> s
```

# Projective measurements

To observe the effect of different projections, we will start with a GHZ state.

```jldoctest proj
julia> s = S"-XXX
             +ZZI
             -IZZ"
```

The `project!` function returns the new stabilizer, the index where the anticommutation was detected,
and the result of the projection (`nothing` being an undetermined result). For instance
here we project on an operator that does not commute with all stabilizer generators.

```jldoctest proj
julia> project!(copy(s), P"ZII")
```

Or we can project on a commuting operator, hence no anticommuting terms (the index is zero),
and the result is perfectly determined (-1, or in our convention to represent the phase, 0x2).

```jldoctest proj
julia> project!(copy(s), P"-ZZI")
```

When the projection is consistent with the stabilizer (i.e. the mesurement result is not `nothing`),
this would trigger an expensive canonicalization procedure in order to calculate the measurement
result (unless we are using more advanced data structure to represent the state).
If all you want to know is whether the projection is consistent with the stabilizer, but
you do not care about the measurement result, you can skip the canonicalization and calculation of
the result.

```jldoctest proj
julia> project!(copy(s), P"-ZZI", keep_result=false)
```

Lastly, in either case, you can skip the calculation of the phases as well, if they are unimportant.

```jldoctest proj
julia> project!(copy(s), P"ZZI", phases=false)
```

# Generating a Pauli operator with Stabilizer generators

i.e. checking for independence. This particular function requires the stabilizer to be already canonicalized.

```jldoctest gen
julia> s = S"-XXX
             +ZZI
             -IZZ";
julia> s = canonicalize!(s)
```

It modifies the Pauli operator in place, reducing it to identity if possible.
The leftover phase is present to indicate if the phase itself could not have been canceled.
The list of indices specifies which rows of the stabilizer were used to generated the
desired Pauli operator.


```jldoctest gen
julia> generate!(P"XYY", s)
```

Phases can be neglected, for higher performance.

```jldoctest gen
julia> generate!(P"XYY", s, phases=false)
```

If the Pauli operator can not be generated by the stabilizer, `nothing` value is returned.

```jldoctest gen
julia> generate!(P"ZZZ", s)
julia> generate!(P"XZX", s)
julia> generate!(P"YYY", s)
```

# Clifford Operators

A number of predefined Clifford operators are available.

```jldoctest
julia> Hadamard
julia> Phase
julia> CNOT
julia> CliffordId
```

Chaining and tensor products are possible (but slow (TODO)). Same for qubit permutations.

```jldoctest
julia> Hadamard ⊗ Phase
julia> Hadamard * Phase
```


```jldoctest
julia> permute(CNOT, [2,1])
```

You can create custom Clifford operators with C-strings or with a list of Pauli
operators. TODO: creating them by using a Stabilizer or by using boolean matrices.


```jldoctest
julia> C"-ZZ
         +_Z
         -X_
         +XX"
julia> CliffordOperator([P"-ZZ", P"_Z", P"-X_", P"XX"])
```

Naturally, the operators can be applied to stabilizer states. This includes
high performance in-place operations (and the phase can be neglected with
`phases=false` for faster computation).

```jldoctest
julia> CNOT * S"X_"
julia> s = S"X_";
julia> apply!(s,CNOT)
```

TODO: Small Clifford operators can be applied to large stabilizers by specifying the qubit indices.

Pauli operators act as Clifford operators too (but they are rather boring, as they only change signs).

```jldoctest
julia> P"XII" * S"ZXX"
```

# Destabilizers

Slightly abusing the name: What we call "destabilizers" here is a stabilizer and
its destabilizing operators saved together. They are initialized from a
stabilizer.


```jldoctest destab
julia> s=S"-XXX
           -ZZI
           +IZZ";
julia> d = calculate_destabilizer(s)
```

With convenience properties to extract only the stabilizer and destabilizer
pieces:

```jldoctest destab
julia> d.stabilizer
julia> d.destabilizer
```

Importantly commuting projections are much faster when tracking the destabilizer
as canonicalization is not necessary.

```jldoctest destab
julia> project!(d,P"ZZI")
```

Non-commuting projections are just as fast as when using only stabilizers.

```jldoctest destab
julia> project!(d,P"ZZZ")
```

Clifford operations can be applied the same way they are applied to stabilizers.

```jldoctest destab
julia> apply!(d,CNOT⊗Hadamard)
```

# TODO: Mixed Stabilizer States

Currently we deal manually with mixed states, as they are not implemented inside the library.
