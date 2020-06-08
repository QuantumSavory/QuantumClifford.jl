# [Choosing Appropriate Data Structure](@id Choosing-Appropriate-Data-Structure)

There are four different data structures used to represent stabilizer states. If
you will never need projective measurements you probably would want to use
[`Stabilizer`](@ref). If you require projective measurements, but only on pure
states, [`Destabilizer`](@ref) should be the appropriate data structure. If
mixed stabilizer states are involved, [`MixedStabilizer`](@ref) would be
necessary.

[`Stabilizer`](@ref) is simply a list of Pauli operators in a tableau form. As a
data structure it does not enforce the requirements for a pure stabilizer state
(the rows of the tableau do not necessarily commute, nor are they forced to be
Hermitian; the tableau might be underdetermined, redundant, or contradictory).
It is up to the user to ensure that the initial values in the tableau are
meaningful and consistent.

[`canonicalize!`](@ref), [`project!`](@ref), and [`generate!`](@ref) can accept
an under determined (mixed state) `Stabilizer` instance and operate correctly.
`canonicalize!` can also accept a redundant `Stabilizer` (i.e. not all rows are
independent), leaving as many identity rows at the bottom of the canonicalized
tableau as the number of redundant stabilizers in the initial tableau.

`canonicalize!` takes ``\mathcal{O}(n^3)`` steps. `generate!` expects a
canonicalized input and then takes ``\mathcal{O}(n^2)`` steps. `project!` takes
``\mathcal{O}(n^3)`` for projecting on commuting operators due to the need to
call `canonicalize!` and `generate!`. If the projections is on an anticommuting
operator (or if `keep_result=false`) then it takes ``\mathcal{O}(n^2)`` steps.

[`MixedStabilizer`](@ref) provides explicit tracking of the rank of the mixed
state and works properly when the projection is on a commuting operator not in
the stabilizer (see table below for details). Otherwise it has the same
performance as `Stabilizer`.

The canonicalization can be made unnecessary if we track the destabilizer
generators. There are two data structures capable of that.

[`Destabilizer`](@ref) stores both the destabilizer and stabilizer states.
`project!` called on it never requires a stabilizer canonicalization, hence it
runs in ``\mathcal{O}(n^2)``. However, `project!` will raise an exception if you
try to project on a commuting state that is not in the stabilizer as that would
be an expensive operation (``\mathcal{O}(n^3)``).

[`MixedDestabilizer`](@ref) tracks both the destabilizer operators and the
logical operators in addition to the stabilizer generators. It does not require
canonicalization for measurements and its `project!` operations always takes
``\mathcal{O}(n^3)``.

| projection | `Stabilizer` | `MixedStabilizer` | `Destabilizer` | `MixedDestabilizer` |
|---|---|---|---|---|
| on anticommuting operator | correct result in ``\mathcal{O}(n^2)`` steps | same as `Stabilizer` | same as `Stabilizer` | same as `Stabilizer` |
| on commuting operator in the stabilizer | ``\mathcal{O}(n^3)``, or ``\mathcal{O}(n^2)`` if `keep_result=false` | ``\mathcal{O}(n^3)`` | ``\mathcal{O}(n^2)`` if the state is pure, throws exception otherwise | ``\mathcal{O}(n^2)`` |
| on commuting operator out of the stabilizer[^1] | ``\mathcal{O}(n^3)``, but the user needs to manually include the new operator to the stabilizer | ``\mathcal{O}(n^3)`` | not applicable if the state is pure, throws an error otherwise | ``\mathcal{O}(n^3)`` |

[^1]:

    This can occur only if the state being projected is mixed. Both `Stabilizer`
    and `Destabilizer` can be used for mixed states by simply providing fewer
    stabilizer generators than qubits at initialization. This can be useful for
    low-level code that tries to avoid the extra memory cost of using
    `MixedStabilizer` and `MixedDestabilizer` but should be avoided otherwise.
    `project!` works correctly or raises an explicit warning on all 4 data
    structures.
