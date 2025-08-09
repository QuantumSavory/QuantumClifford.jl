const THROW_BOUNDS =
"Unable to perform the requested operation due to encountering a mismatch \
between the provided index(indices) and the pertinent size(s) of the other \
argument(s)."

const THROW_SIZE =
"Unable to perform the requested operation due to encountering a mismatch \
between the pertinent size(s) of the provided arguments."

const THROW_NQUBITS =
"Unable to perform the requested operation due to encountering a mismatch \
between the number of qubits in the provided arguments."

const THROW_REPR_NEEDS_COMMUTATIVE_ALGEBRA = 
"The group algebra must be commutative when using a single `repr` function,\
which is not the case here. Please specify separate `A_repr` and `B_repr` \
instead of a single `repr`. The default choice of\
`A_repr=right_repr_matrix, B_repr=left_repr_matrix` is frequently sufficient."

const THROW_MUST_BE_DIVISIBLE_BY_3 =
"The input parameters must be divisible by 3"

const THROW_REQUIRED_POSITIVE_ARG = 
"The input parameter should be positive"

const THROW_INVALID_XZ_COMPONENTS = 
"`xzcomponents` should be `:split` or `:together`"

const THROW_INVALID_CODE_DIMENSION =
"Dimension of the input code must be at least 2"    

const THROW_INVALID_X_METACHECKS =
"`X`-metachecks (`Mx`) require `D ≥ 4`"    

const THROW_INVALID_Z_METACHECKS = 
"`Z`-metachecks (`Mz`) require `D ≥ 3`"    


const THROW_INVALID_CUDA_ARG = 
"First argument to @run_cuda should be a function call"

const THROW_MODEL_MEMORY_LIMIT =
"Model exceeded memory limits"

const THROW_MODEL_TIME_LIMIT =
"Model exceeded time limit"

const THROW_INVALID_PARAMETERS_GOLAY =
"Invalid parameters: `n` must be either 24 or 23 to obtain a valid code."

const THROW_INVALID_PARAMETERS_HAMMING =
"Invalid parameters: `r` must be ≥ 2 to obtain a valid code."

const THROW_INVALID_PARAMETERS_REEDMULLER =
"Invalid parameters: r must be non-negative and r ≤ m in order to valid code."

const THROW_COLOR_CODES_ODD =
"Only odd distance triangular color codes are allowed.\nRefer to https://arxiv.org/abs/1108.5738"

const THROW_COLOR_CODES_MIN_DIST =
"Smallest allowed distance is 3.\nRefer to https://arxiv.org/abs/1108.5738"

const THROW_DELFOSSE_823_MIN_BLOCKS =
"The number of blocks must be at least 1 to construct a valid code."

const THROW_DELFOSSE_MIN_BLOCKS =
"The number of blocks must be at least 2 to construct a valid code."

const THROW_DELFOSSE_SELF_ORTHOGONAL =
"The `Reed-Muller` parity check matrix must be 'self-orthogonal' to construct a self-dual
CSS `DelfosseReichardt` code. Use `search_self_orthogonal_rm_codes` to search for good parameters for `Reed-Muller` codes
that provide `self-orthogonal` seeds."

const THROW_DELFOSSE_REP_MULTIPLE =
"The number of blocks must be a multiple of 2."

const THROW_INVALID_QUANTUM_REED_MULLER =
"Invalid parameters: m must be bigger than 2 in order to have a valid code."

const THROW_INVALID_PARAMETERS_TILLICH_ZEMOR =
"Conditions for the existence of `M` in `H = [C | M]` are not satisfied."

const THROW_INVALID_PARAMETERS_CLIFFORD = 
"Input tableau should be of size 2n×n (top half is the X mappings and the bottom half are the Z mappings)."

const THROW_INVALID_ACTION_QUBITS =
"The tableau and the provided operator need to act on the same number of qubits.
Consider specifying an array of indices as a third argument to the `apply!` function to avoid this error."

const THROW_INVALID_PARAMETERS_BCH_ARG = 
"""Conditions for valid BCH(m, t) arguments:
- m ≥ 3
- t < 2ᵐ⁻¹
- m*t ≤ 2ᵐ - 1
"""

const THROW_MISSING_HECKE =
"You've invoked a function which depends on the package `Hecke` but you have not installed or imported it yet.
Immediately after you import `Hecke`, the function will be available."

const THROW_MISSING_OSCAR =
"You've invoked a function which depends on the package `Oscar` but you have not installed or imported it yet.
Immediately after you import `Oscar`, the function will be available."

const THROW_MISSING_LDPCDecoders =
"You've invoked a function which depends on the package `LDPCDecoders` but you have not installed or imported it yet.
Immediately after you import `LDPCDecoders`, the function will be available."

const THROW_MISSING_PyQDecoders =
"You've invoked a function which depends on the package `PyQDecoders` but you have not installed or imported it yet.
Immediately after you import `PyQDecoders`, the function will be available."

const THROW_INVALID_PARAMETERS_NOISE =
"The input noise should be between 0 and 1."

const THROW_INVALID_LOGICAL_OPERATOR =
"`logical_operator_type` must be :X or :Z"

const THROW_CHECKS_MISSING =
"Codes of the type used do not have separate X and Z parity checks, either because they are not a CSS code
and thus inherently do not have separate checks, or because its separate checks are not yet implemented in this library."

const THROW_INCONSISTENT_TABLEAU =
"the tableau is inconsistent (check if it is clip-canonicalized and Hermitian)"

const THROW_PHASE_FUNCTIONALITY_UNAVAILABLE_GROUPIFY =
"in groupify phases=true functionality not yet implemented"

const THROW_MIXED_FUNCTIONALITY_UNAVAILABLE_DOTPROD =
"Only pure (not mixed) states are supported when calculating inner product."

const THROW_BAD_DATA_STRUCTURE =
"You've used a type(`Stabilizer` or `Destabilizer`) does not permit automatic tracking of the rank.
Use `length`, `trusted_rank`, the `MixedDestabilizer` type, or track the rank manually."

const THROW_NOT_STABILIZER_STATE =
"""The argument you have provided for good_state is not a logical state within the codespace.
Expected a pure qubit stabilizer state (i.e. independent stabilizer generators on be equal to number of qubits)"""

const THROW_MIXED_FUNCTIONALITY_UNAVAILABLE_STABILIZER =
"""Attempting to convert a `Stabilizer`-like object to `GeneralizedStabilizer` object failed, 
because the initial state does not represent a pure state. 
Currently only pure states can be used to initialize a `GeneralizedStabilizer` mixture of stabilizer states."""

const THROW_PAULI_BOUNDS =
"""
You are attempting to construct a `PauliChannel` but have provided Pauli operators
that are not all of the same size (same number of qubits).
Please ensure that all of the Pauli operators being provided of of the same size.
"""

const THROW_EMBEDDING_REQUIRED =
"The number of qubits in the GeneralizedStabilizer state does not match the number of qubits in the Pauli channel.
 Use `embed` to pad the PauliChannel so it acts on the correct number of qubits."

const THROW_INEFFICIENT_DESTABILIZER =
"`Destabilizer` can not efficiently (faster than n^3) detect whether you are projecting on a stabilized or a logical operator.
Switch to one of the `Mixed*` data structures."

const THROW_ONLY_FULL_RANK_DESTABILIZER_SUPPORTED =
"Only full-rank `Destabilizer` can be converted to `MixedDestabilizer` with specific rank. Try not specifying `r`."

const THROW_PHASES =
"Only {±1,±i} are permitted as phases."

const THROW_ROW_MISMATCH_TABLEAUX =
"All input Tableaux/Stabilizers musT have the same number of rows."

const THROW_NO_ZERO_QUBIT = 
"Qubit indices have to be larger than zero, 
but you are attempting to create a gate acting on a qubit with a non-positive index. Ensure indexing always starts from 1."

const THROW_INVALID_CONVERSION_TO_SINGLEQUBIT =
"You are trying to convert a multiqubit `CliffordOperator` into a symbolic `SingleQubitOperator`."

const THROW_INVALID_QUBIT_OPERATOR_SIZE_ALLOCATION =
"Set a larger `n`, otherwise the `SingleQubitOperator` or `TwoQubitOperator`(whichever is being used)  can not fit in the allocated `CliffordOperator`.
 Use `compact=true` if you want to discard the target index."

const THROW_INVALID_TWO_QUBIT_GATE =
"Failed to create a two qubit gate because the two qubits it acts on have the same index. The qubit indices have to be different."

const THROW_COMPACT_QUBIT_OPERATOR_SHAPE_ERROR =
"Set `n=k` when compacting (`compact=true`) a `k`-qubit operator, so that it results in a `k×k` `CliffordOperator`."

const THROW_INVALID_BASIS_REP =
"`basis` should be one of :X, :Y, or :Z"

