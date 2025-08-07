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