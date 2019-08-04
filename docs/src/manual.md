# Manual

```@meta
DocTestSetup = quote
    using SimpleClifford
end
```

# Pauli Operators

They are stored in memory as a phase (a single byte where `0x0,0x1,0x2,0x3` corresponds to $1,i,-1,-i$) and two bit-arrays, for X and for Z components.

You can create them with a `P` string.

```jldoctest
julia> P"-iXZ"
-iXZ
```

Or by specifying phase and X/Z components:

```jldoctest
julia> PauliOperator(0x0,Bool[0,1,0],Bool[0,0,1])
+ _XZ
```

Both underscore and I can be used for identity.

```jldoctest
julia> P"I_XYZ"
+ __XYZ
```

Multiplication with scalars or other Pauli operators works as expected,
as well as tensor products of Pauli operators.

```jldoctest
julia> -1im*P"X"
-iX

julia> P"X" * P"Z"
-iY

julia> P"X" ⊗ P"Z"
+ XZ
```

One can check for commutativity with `comm`.

```jldoctest
julia> comm(P"X",P"Z")
0x01

julia> comm(P"XX",P"ZZ")
0x00
```

And check the phase of a product with `prodphase`.

```jldoctest
julia> prodphase(P"X", P"Z")
0x03

julia> prodphase(P"X", P"iZ")
0x00

julia> prodphase(P"X",P"Y")
0x01
```

Indexing operations are available.

```jldoctest
julia> p = P"IXYZ";

julia> p[1], p[2], p[3], p[4]
((false, false), (true, false), (true, true), (false, true))

julia> p = P"III";

julia> p[2] = (true, true);

julia> p
+ _Y_
```
Including various permutation methods and fancy indexing:

```jldoctest
julia> permute(P"XYZ", [3,1,2])
+ ZXY

julia> P"IXYZ"[[2,3]]
+ XY
```

The operator is represented in memory by bit arrays (much denser than using byte
arrays).

```jldoctest
julia> p = P"-IXYZ";

julia> p.nqbits, p.phase, p.xz
(4, 0x02, UInt64[0x0000000000000006, 0x000000000000000c])
```

The convenience properties `xbit` and `zbit` give you Bool (GF2) vectors.
TODO: this should be a separate function.

```jldoctest
julia> P"XYZI".xbit
4-element BitArray{1}:
 1
 1
 0
 0
```

# Stabilizers

They are stored in memory as a phase list and a bit-matrix for X and Z
components. They can be created by an `S` string, or with a number of different
constructors.

```jldoctest
julia> S"-XX
         +ZZ"
- XX
+ ZZ

julia> Stabilizer([P"-XX",P"+ZZ"])
- XX
+ ZZ

julia> Stabilizer([0x2, 0x0],
                  Bool[1 1;
                       0 0],
                  Bool[0 0;
                       1 1])
- XX
+ ZZ
```


Direct sums can be performed,

```jldoctest
julia> S"-XX" ⊕ S"ZZ"
- XX__
+ __ZZ
```

Indexing operations are available, including fancy indexing. 2D indexing,
into the Pauli operators is not implemented (TODO).

```jldoctest
julia> #TODO
```

Consistency at creation is not verified so nonsensical stabilizers can be
created, both in terms of content and shape.

```jldoctest
julia> S"iX
         +Z"
+iX
+ Z
```

Similarly to the Pauli operators, a bit array representation is used.

```jldoctest stab
julia> s = S"-XXX
             +ZZI
             -IZZ"
- XXX
+ ZZ_
- _ZZ

julia> s.phases, s.nqbits, s.xzs
(UInt8[0x02, 0x00, 0x02], 3, UInt64[0x0000000000000007 0x0000000000000000; 0x0000000000000000 0x0000000000000003; 0x0000000000000000 0x0000000000000006])
```

And there are convenience functions that can extract the corresponding binary
check matrix.

```jldoctest stab
julia> stab_to_gf2(s)
3×6 BitArray{2}:
 1  1  1  0  0  0
 0  0  0  1  1  0
 0  0  0  0  1  1
```

# Canonicalization of Stabilizers

Canonicalization (akin to Gaussian elimination over F(2,2)) is implemented.

```jldoctest
julia> s = S"-XXX
             +ZZX
             +III";

julia> canonicalize!(s)
+ YY_
+ ZZX
+ ___
```

If phases are inconsequential, the operations can be faster by not tracking and updating them.

```jldoctest
julia> s = S"-XXX
             +ZZX
             +III";

julia> canonicalize!(s; phases=false)
- YY_
+ ZZX
+ ___
```

These operations are in place (as customarily signified by "!").

```jldoctest
julia> s = S"-XXX
             +ZZX
             +III";

julia> canonicalize!(s; phases=false);

julia> s
- YY_
+ ZZX
+ ___
```

# Projective measurements

To observe the effect of different projections, we will start with a GHZ state.

```jldoctest proj
julia> s = S"-XXX
             +ZZI
             -IZZ";
```

The `project!` function returns the new stabilizer, the index where the anticommutation was detected,
and the result of the projection (`nothing` being an undetermined result). For instance
here we project on an operator that does not commute with all stabilizer generators.

```jldoctest proj
julia> project!(copy(s), P"ZII")
(+ Z__
+ ZZ_
- _ZZ, 1, nothing)
```

Or we can project on a commuting operator, hence no anticommuting terms (the index is zero),
and the result is perfectly determined (-1, or in our convention to represent the phase, 0x2).

```jldoctest proj
julia> project!(copy(s), P"-ZZI")
(- XXX
- Z_Z
- _ZZ, 0, 0x02)
```

When the projection is consistent with the stabilizer (i.e. the mesurement result is not `nothing`),
this would trigger an expensive canonicalization procedure in order to calculate the measurement
result (unless we are using more advanced data structure to represent the state).
If all you want to know is whether the projection is consistent with the stabilizer, but
you do not care about the measurement result, you can skip the canonicalization and calculation of
the result.

```jldoctest proj
julia> project!(copy(s), P"-ZZI", keep_result=false)
(- XXX
+ ZZ_
- _ZZ, 0, nothing)
```

Lastly, in either case, you can skip the calculation of the phases as well, if they are unimportant.

```jldoctest proj
julia> project!(copy(s), P"ZZI", phases=false)
(- XXX
+ Z_Z
- _ZZ, 0, 0x00)
```

# Generating a Pauli operator with Stabilizer generators

i.e. checking for independence. This particular function requires the stabilizer to be already canonicalized.

```jldoctest gen
julia> s = S"-XXX
             +ZZI
             -IZZ";

julia> s = canonicalize!(s)
- XXX
- Z_Z
- _ZZ
```

It modifies the Pauli operator in place, reducing it to identity if possible.
The leftover phase is present to indicate if the phase itself could not have been canceled.
The list of indices specifies which rows of the stabilizer were used to generated the
desired Pauli operator.


```jldoctest gen
julia> generate!(P"XYY", s)
(- ___, [1, 3])
```

Phases can be neglected, for higher performance.

```jldoctest gen
julia> generate!(P"XYY", s, phases=false)
(+ ___, [1, 3])
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
X ⟼ + Z
Z ⟼ + X

julia> Phase
X ⟼ + Y
Z ⟼ + Z

julia> CNOT
X_ ⟼ + XX
_X ⟼ + _X
Z_ ⟼ + Z_
_Z ⟼ + ZZ

julia> CliffordId
X ⟼ + X
Z ⟼ + Z
```

Chaining and tensor products are possible (but slow (TODO)). Same for qubit permutations.

```jldoctest
julia> Hadamard ⊗ Phase
X_ ⟼ + Z_
_X ⟼ + _Y
Z_ ⟼ + X_
_Z ⟼ + _Z

julia> Hadamard * Phase
X ⟼ - Y
Z ⟼ + X
```


```jldoctest
julia> permute(CNOT, [2,1])
X_ ⟼ + X_
_X ⟼ + XX
Z_ ⟼ + ZZ
_Z ⟼ + _Z
```

You can create custom Clifford operators with C-strings or with a list of Pauli
operators. TODO: creating them by using a Stabilizer or by using boolean matrices.


```jldoctest
julia> C"-ZZ
         +_Z
         -X_
         +XX"
X_ ⟼ - ZZ
_X ⟼ + _Z
Z_ ⟼ - X_
_Z ⟼ + XX

julia> CliffordOperator([P"-ZZ", P"_Z", P"-X_", P"XX"])
X_ ⟼ - ZZ
_X ⟼ + _Z
Z_ ⟼ - X_
_Z ⟼ + XX
```

Naturally, the operators can be applied to stabilizer states. This includes
high performance in-place operations (and the phase can be neglected with
`phases=false` for faster computation).

```jldoctest
julia> CNOT * S"X_"
+ XX

julia> s = S"X_";

julia> apply!(s,CNOT)
+ XX
```

TODO: Small Clifford operators can be applied to large stabilizers by specifying the qubit indices.

Pauli operators act as Clifford operators too (but they are rather boring, as they only change signs).

```jldoctest
julia> P"XII" * S"ZXX"
- ZXX
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
+ Z__
+ X__
+ _X_
━━━━━
- XXX
- Z_Z
+ _ZZ
```

With convenience properties to extract only the stabilizer and destabilizer
pieces:

```jldoctest destab
julia> d.stabilizer
- XXX
- Z_Z
+ _ZZ

julia> d.destabilizer
+ Z__
+ X__
+ _X_
```

Importantly commuting projections are much faster when tracking the destabilizer
as canonicalization is not necessary.

```jldoctest destab
julia> project!(d,P"ZZI")
(+ Z__
+ X__
+ _X_
━━━━━
- XXX
- Z_Z
+ _ZZ, 0, 0x02)
```

Non-commuting projections are just as fast as when using only stabilizers.

```jldoctest destab
julia> project!(d,P"ZZZ")
(- XXX
+ _XX
+ X_X
━━━━━
+ ZZZ
- Z_Z
+ _ZZ, 1, nothing)
```

Clifford operations can be applied the same way they are applied to stabilizers.

```jldoctest destab
julia> apply!(d,CNOT⊗Hadamard)
- X_Z
+ _XZ
+ XXZ
━━━━━
+ _ZX
- Z_X
+ ZZX
```

# TODO: Mixed Stabilizer States

Currently we deal manually with mixed states, as they are not implemented inside the library.
