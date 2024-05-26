# [Stabilizer Tableau Algebra Manual](@id Stabilizer-Tableau-Algebra-Manual)

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

The library consists of two main parts: Tools for working with the algebra of Stabilizer tableaux and tools specifically for efficient Circuit Simulation. This chapter discusses the former "lower level" Stabilizer tableau algebra tools.

# Pauli Operators

The [`PauliOperator`](@ref) object represents multi-qubit Pauli operator
(``±\{1,i\}\{I,Z,X,Y\}^{\otimes n}``). It is stored in memory as a phase (a
single byte where `0x0,0x1,0x2,0x3` corresponds to $1,i,-1,-i$) and two
bit-arrays, for X and for Z components.

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

Multiplication with scalars or other Pauli operators works as expected, as well
as tensor products of Pauli operators.

```jldoctest
julia> -1im*P"X"
-iX

julia> P"X" * P"Z"
-iY

julia> P"X" ⊗ P"Z"
+ XZ
```

One can check for commutativity with [`comm`](@ref).

```jldoctest
julia> comm(P"X",P"Z")
0x01

julia> comm(P"XX",P"ZZ")
0x00
```

And check the phase of a product with [`prodphase`](@ref).

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

Including fancy indexing:

```jldoctest
julia> P"IXYZ"[[2,3]]
+ XY

julia> P"IXYZ"[[false,true,true,false]]
+ XY
```

The operator is represented in memory by bit arrays (much denser than using byte
arrays).

```jldoctest
julia> p = P"-IXYZ";

julia> p.nqubits, p.xz
(4, UInt64[0x0000000000000006, 0x000000000000000c])
```

Views that give just the X or Z components of the `xz` bitarray are available through [`xview`](@ref) and [`zview`](@ref).

```jldoctest
julia> xview(P"XYZI")
1-element view(::Vector{UInt64}, 1:1) with eltype UInt64:
 0x0000000000000003
```

The convenience methods [`xbit`](@ref) and [`zbit`](@ref) give you Bool (GF2) vectors.

```jldoctest
julia> xbit(P"XYZI")
4-element Vector{Bool}:
 1
 1
 0
 0
```

# [Stabilizers](@id Stabilizers)

A [`Stabilizer`](@ref) object is a tableau of Pauli operators. When the tableau is
meant to represent a (pure or mixed) stabilizer state, all of these operators
should commute (but that is not enforced, rather `Stabilizer` is a generic
tableau data structure). It is stored in memory as a phase list and a bit-matrix
for X and Z components. It can be instantiated by an `S` string, or with a
number of different constructors.

!!! tip "Stabilizers and Destabilizers"
    In many cases you probably would prefer to use the [`MixedDestabilizer`](@ref)
    data structure, as it caries a lot of useful additional information, like tracking
    rank and destabilizer operators. `Stabilizer` has mostly a pedagogical value, and it
    is also used for slightly faster simulation of a particular subset of Clifford
    operations.
    See also the [data structures discussion page](@ref Choosing-Appropriate-Data-Structure).

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

julia> Stabilizer(Bool[1 1 0 0;
                       0 0 1 1])
+ XX
+ ZZ
```

Direct sums can be performed,

```jldoctest
julia> S"-XX" ⊗ S"ZZ"
- XX__
+ __ZZ
```

Indexing operations are available, including fancy indexing. Be careful about how phase information gets transferred during sub-indexing.

```jldoctest
julia> s = S"-XYZ
             -ZIX
             +XIZ";

julia> s[1]
- XYZ

julia> s[1,2]
(true, true)

julia> s[[3,1]]
+ X_Z
- XYZ

julia> s[[3,1],[2]]
+ _
- Y

julia> s[3][3]
(false, true)

julia> getindex(s, 1)
- XYZ

julia> getindex(s, 3, 1)
(true, false)

julia> setindex!(s, P"Z", 1)
+ Z__
- Z_X
+ X_Z

julia> setindex!(s, P"ZYX", 1)
+ ZYX
- Z_X
+ X_Z

julia> setindex!(s, (true, true), 1, 1)
+ YYX
- Z_X
+ X_Z

julia> firstindex(s)
1

julia> lastindex(s)
3

julia> lastindex(s, 2)
3
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

julia> phases(s), tab(s).xzs
(UInt8[0x02, 0x00, 0x02], UInt64[0x0000000000000007 0x0000000000000000 0x0000000000000000; 0x0000000000000000 0x0000000000000003 0x0000000000000006])
```

And there are convenience functions that can extract the corresponding binary
check matrix.

```jldoctest stab
julia> stab_to_gf2(s)
3×6 Matrix{Bool}:
 1  1  1  0  0  0
 0  0  0  1  1  0
 0  0  0  0  1  1
```

Stabilizer is copied and assigned.
```jldoctest stabilizer
julia> s₁ = copy(s)
- XXX
+ ZZ_
- _ZZ
```

# [Canonicalization of Stabilizers](@id Canonicalization-of-Stabilizers)

Canonicalization (akin to Gaussian elimination over F(2,2)) is implemented in
the [`canonicalize!`](@ref) function.
Besides the default canonicalization prescription,
alternative ones are available as described in the
[canonicalization page](@ref Canonicalization-operations).

```jldoctest
julia> s = S"-XXX
             +ZZX
             +III";

julia> canonicalize!(s)
+ YY_
+ ZZX
+ ___
```

If phases are inconsequential, the operations can be faster by not tracking and
updating them.

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

# [Projective Measurements](@id Projective-Measurements)

The [`project!`](@ref) function is used to perform generic projective measurements.

!!! tip "Single qubit projections"
    If you know your Pauli measurement operator acts on a single qubit, there are
    much faster projection functions available, discussed in the next section.
    Namely [`projectX!`](@ref), [`projectY!`](@ref), and [`projectZ!`](@ref).

To observe the effect of different projections, we will start with a GHZ state.

```jldoctest proj
julia> s = S"-XXX
             +ZZI
             -IZZ";
```

The [`project!`](@ref) function returns the new stabilizer, the index where the
anticommutation was detected, and the result of the projection (`nothing` being
an undetermined result). For instance here we project on an operator that does
not commute with all stabilizer generators.

```jldoctest proj
julia> project!(copy(s), P"ZII")[1]
+ Z__
+ ZZ_
- _ZZ
```

Importantly, when there is an undetermined result, we return `nothing` **and
leave the phase of the new stabilizer the same as the phase of the projection
operator**. If you want to perform a Monte Carlo simulation, you need to
randomize the phase of the stabilizer at the anticommuting index yourself. For
instance, one can do:

```jldoctest proj
julia> newstate, anticomindex, result = project!(copy(s), P"XII")
       if isnothing(result)
           phases(newstate)[anticomindex] = rand([0x0,0x2])
       end
       result, anticomindex
(nothing, 2)
```

Of course, this is a rather cumbersome way to run a simulation, so we also provide
[`projectrand!`](@ref) which does the necessary randomization automatically,
for cases where you do not need the fine grained control of `project!`.

We can project on a commuting operator, hence no anticommuting terms (the
index is zero), and the result is perfectly determined (-1, or in our convention
to represent the phase, 0x2).

```jldoctest proj
julia> project!(copy(s), P"-ZZI")
(Stabilizer 3×3, 0, 0x02)
```

When the projection is consistent with the stabilizer (i.e. the measurement
result is not `nothing`), this would trigger an expensive canonicalization
procedure in order to calculate the measurement result (unless we are using more
advanced data structures to represent the state, which are discussed later). If
all you want to know is whether the projection is consistent with the
stabilizer, but you do not care about the measurement result, you can skip the
canonicalization and calculation of the result.

```jldoctest proj
julia> project!(copy(s), P"-ZZI", keep_result=false)
(Stabilizer 3×3, 0, nothing)
```

Lastly, in either case, you can skip the calculation of the phases as well, if
they are unimportant.

```jldoctest proj
julia> project!(copy(s), P"ZZI", phases=false)
(Stabilizer 3×3, 0, 0x00)
```

## Sparse single-qubit measurements

In many circumstances only a single-qubit operator is being measured. In that case one should use the [`projectX!`](@ref), [`projectY!`](@ref), and [`projectZ!`](@ref) functions as they are much faster thanks to tracking only a single qubit.
They have versions that randomize the phase as necessary as well:  [`projectXrand!`](@ref), [`projectYrand!`](@ref), and [`projectZrand!`](@ref).

## Gate-like interface

If you do not need all this boilerplate, and especially if you want to perform the randomization automatically, you can use the gate-like "symbolic" objects [`sMX`](@ref), [`sMY`](@ref), and [`sMZ`](@ref), that perform the measurement and the necessary randomization of phase. If the measurement result is to be stored, you can use the [`Register`](@ref) structure that stores both stabilizer tableaux and bit values.

```
julia> state = Register(ghz(3), [false,false])
Register{Vector{UInt8}, Matrix{UInt64}}(Rank 3 stabilizer
+ Z__
+ _X_
+ __X
═════
+ XXX
+ ZZ_
+ Z_Z
═════
, Bool[0, 0])

julia> apply!(state, sMX(3,2)) # which qubit is measured and in which bit it is stored
Register{Vector{UInt8}, Matrix{UInt64}}(Rank 3 stabilizer
+ Z__
+ _X_
+ Z_Z
═════
+ XXX
+ ZZ_
- __X
═════
, Bool[0, 1])

julia> bitview(state)
2-element Vector{Bool}:
 0
 1
```

Or you can use the [`projectXrand!`](@ref), [`projectYrand!`](@ref), and [`projectZrand!`](@ref) if you prefer a function-call interface.


# [Partial Traces](@id Partial-Traces)

Partial trace (using [`traceout!`](@ref)) over even a single qubit might cause many of them to decohere due to entanglement.

```jldoctest
julia> ghz = S"XXX
               ZZ_
               _ZZ";

julia> traceout!(ghz, [1])
+ _ZZ
+ ___
+ ___
```

This is somewhat more elegant when the datastructure being used explicitly supports mixed states.

```jldoctest
julia> ghz = MixedStabilizer(S"XXX
                               ZZ_
                               _ZZ");

julia> traceout!(ghz, [1])
+ _ZZ
```


# [Generating a Pauli Operator with Stabilizer Generators](@id Generating-a-Pauli-Operator-with-Stabilizer-Generators)

The [`generate!`](@ref) function attempts to generate a Pauli operator by
multiplying together the operators belonging to a given stabilizer (or reports
their independence). This particular function requires the stabilizer to be
already canonicalized.

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
The leftover phase is present to indicate if the phase itself could not have
been canceled. The list of indices specifies which rows of the stabilizer were
used to generated the desired Pauli operator.

```jldoctest gen
julia> generate!(P"XYY", s)
(- ___, [1, 3])
```

Phases can be neglected, for higher performance.

```jldoctest gen
julia> generate!(P"XYY", s, phases=false)
(+ ___, [1, 3])
```

If the Pauli operator can not be generated by the stabilizer, `nothing` value is
returned.

```jldoctest gen
julia> generate!(P"ZZZ", s)

julia> generate!(P"XZX", s)

julia> generate!(P"YYY", s)
```

# [Clifford Operators](@id Clifford-Operators)

The [`CliffordOperator`](@ref) structure represents a linear mapping between
stabilizers (which should also preserve commutation relationships, but that is
not checked at instantiation).
These are n-qubit dense tableaux, representing an operation on n-qubit states. For single- or two-qubit gates, it is much more efficient to use small sparse [symbolic clifford operators](@ref Symbolic-Clifford-Operators).
A number of predefined Clifford operators are available,
their name prefixed with `t` to mark them as dense tableaux.

```jldoctest
julia> tHadamard
X₁ ⟼ + Z
Z₁ ⟼ + X

julia> tPhase
X₁ ⟼ + Y
Z₁ ⟼ + Z

julia> tCNOT
X₁ ⟼ + XX
X₂ ⟼ + _X
Z₁ ⟼ + Z_
Z₂ ⟼ + ZZ

julia> tId1
X₁ ⟼ + X
Z₁ ⟼ + Z
```

Chaining and tensor products are possible. Same for qubit permutations.

```jldoctest
julia> tHadamard ⊗ tPhase
X₁ ⟼ + Z_
X₂ ⟼ + _Y
Z₁ ⟼ + X_
Z₂ ⟼ + _Z

julia> tHadamard * tPhase
X₁ ⟼ - Y
Z₁ ⟼ + X

julia> permute(tCNOT, [2,1])
X₁ ⟼ + X_
X₂ ⟼ + XX
Z₁ ⟼ + ZZ
Z₂ ⟼ + _Z
```

You can create custom Clifford operators with C-strings or with a list of Pauli
operators.

```jldoctest
julia> C"-ZZ
         +_Z
         -X_
         +XX"
X₁ ⟼ - ZZ
X₂ ⟼ + _Z
Z₁ ⟼ - X_
Z₂ ⟼ + XX

julia> CliffordOperator([P"-ZZ", P"_Z", P"-X_", P"XX"])
X₁ ⟼ - ZZ
X₂ ⟼ + _Z
Z₁ ⟼ - X_
Z₂ ⟼ + XX
```

Naturally, the operators can be applied to stabilizer states. This includes high
performance in-place operations (and the phase can be neglected with
`phases=false` for faster computation).

```jldoctest
julia> tCNOT * S"X_"
+ XX

julia> s = S"X_";

julia> apply!(s,tCNOT)
+ XX
```

Sparse applications where a small Clifford operator is applied only on a particular subset of a larger stabilizer is also possible, but in such circumstances it is useful to consider using [symbolic operators](@ref Symbolic-Clifford-Operators) too.

```jldoctest
julia> s = S"Z_YX";

julia> apply!(s, tCNOT, [4,2]) # Apply the CNOT on qubits 4 and 2
+ ZXYX
```

Pauli operators act as Clifford operators too (but they are rather boring, as
they only change signs).

```jldoctest
julia> P"XII" * S"ZXX"
- ZXX
```

Internally, the `CliffordOperator` structure simply stores the tableau representation of the operation.

The `apply!` function is efficiently multithreaded for `CliffordOperators`. To start multithreaded Julia, use `julia -t<N>`
where `<N>` specifies the number of threads.

# [Symbolic Clifford Operators](@id Symbolic-Clifford-Operators)

Much faster implementations for a number of common Clifford operators are available. They are stored as special
named structs, instead of as a full tableau. These are the subtypes of `AbstractSingleQubitOperator` and
`AbstractTwoQubitOperator`. Currently these are:

```@example subtypes
using QuantumClifford # hide
using InteractiveUtils # hide
subtypes(QuantumClifford.AbstractSingleQubitOperator)
```

```@example subtypes
subtypes(QuantumClifford.AbstractTwoQubitOperator)
```

Generally, they have the prefix `s` for symbolic/small/sparse.
They are used slightly differently, as one needs to specify the qubits on which they act while instantiating them:

```jldoctest
julia> sHadamard(2)
sHadamard on qubit 2
X₁ ⟼ + Z
Z₁ ⟼ + X

julia> sHadamard(2)*S"XXX"
+ XZX

julia> sCNOT(2,3)*S"XYY"
- XXZ
```

The `apply!` function is efficiently multithreaded for these symbolic operators as well. To start multithreaded Julia, use `julia -t<N>`
where `<N>` specifies the number of threads.

Symbolic projectors on single qubits also exist: [`sMX`](@ref), [`sMY`](@ref), [`sMZ`](@ref). When used with the [`Register`](@ref) state representation, they can store the measurement results in the corresponding classical register.

# Destabilizers

Slightly abusing the name: What we call "destabilizers" here is a stabilizer and
its destabilizing operators saved together. They are implemented with the
[`Destabilizer`](@ref) object and are initialized from a stabilizer.

```jldoctest destab
julia> s=S"-XXX
           -ZZI
           +IZZ";

julia> d = Destabilizer(s)
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ _X_
+ __X
𝒮𝓉𝒶𝒷━
- XXX
- ZZ_
- Z_Z
```

They have convenience methods to extract only the stabilizer and destabilizer
pieces:

```jldoctest destab
julia> stabilizerview(d)
- XXX
- ZZ_
- Z_Z

julia> destabilizerview(d)
+ Z__
+ _X_
+ __X
```

Importantly commuting projections are much faster when tracking the destabilizer
as canonicalization is not necessary (an ``\mathcal{O}(n^2)`` complexity because it avoids
the expensive ``\mathcal{O}(n^3)`` canonicalization operation).

```jldoctest destab
julia> project!(d,P"ZZI")
(Destablizer 3×3, 0, 0x02)
```

Non-commuting projections are just as fast as when using only stabilizers.

```jldoctest destab
julia> project!(d,P"ZZZ")
(Destablizer 3×3, 1, nothing)
```

Clifford operations can be applied the same way they are applied to stabilizers.

```jldoctest destab
julia> apply!(d,tCNOT⊗tHadamard)
𝒟ℯ𝓈𝓉𝒶𝒷
- X_Z
+ XXZ
+ X__
𝒮𝓉𝒶𝒷━
+ _ZX
- _Z_
- Z_X
```

# Mixed States

Both the `Stabilizer` and `Destabilizer` structures have more general forms
that enable work with mixed stabilizer states.
They are the [`MixedStabilizer`](@ref) and [`MixedDestabilizer`](@ref) structures,
described in [Mixed States](@ref Mixed-Stabilizer-States).
More information that can be seen in the [data structures page](@ref Choosing-Appropriate-Data-Structure),
which expands upon the algorithms available for each structure.

# Random States and Circuits

[`random_clifford`](@ref), [`random_stabilizer`](@ref), and [`enumerate_cliffords`](@ref) can be used for the generation of random states.
