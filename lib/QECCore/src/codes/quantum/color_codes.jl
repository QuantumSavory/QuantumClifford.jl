abstract type ColorCode <: AbstractCSSCode end

"""Planar color codes that encode a single logical qubit."""
abstract type TriangularCode <: ColorCode end

"""
    $TYPEDEF

Triangular code following the `4.8.8` tiling. Constructor take a distance `d` as input.

### Example

Here is `[[17,1, 5]]` color code following the `4.8.8` tiling:

```jldoctest
julia> import HiGHS; import JuMP; # hide

julia> using QuantumClifford.ECC: Triangular666, DistanceMIPAlgorithm; # hide

julia> c = Triangular488(5);

julia> code = Stabilizer(c)
+ XXXX_____________
+ X_X_XX___________
+ __XX_XX__XX__XX__
+ ____XX__XX_______
+ ______XX__XX_____
+ _______X___X___XX
+ ________XX__XX___
+ __________XX__XX_
+ ZZZZ_____________
+ Z_Z_ZZ___________
+ __ZZ_ZZ__ZZ__ZZ__
+ ____ZZ__ZZ_______
+ ______ZZ__ZZ_____
+ _______Z___Z___ZZ
+ ________ZZ__ZZ___
+ __________ZZ__ZZ_

julia> distance(c, DistanceMIPAlgorithm(solver=HiGHS))
5
```

More information can be seen in [landahl2011color](@cite)

### Fields
    $TYPEDFIELDS
"""
struct Triangular488 <: TriangularCode
    """The distance of the code."""
    d::Int

    function Triangular488(d)
        if d%2!=1
            throw(ArgumentError(THROW_COLOR_CODES_ODD))
        elseif d<3
            throw(ArgumentError(THROW_COLOR_CODES_MIN_DIST))
        end
        return new(d)
    end
end

"""
    $TYPEDEF

Triangular code following the `6.6.6` tiling. Constructor take a distance `d` as input.

### Example

Here is `[[19,1, 5]]` color code following the `6.6.6` tiling:

```jldoctest
julia> import HiGHS; import JuMP; # hide

julia> using QuantumClifford.ECC: Triangular666, DistanceMIPAlgorithm; # hide

julia> c = Triangular666(5);

julia> code = Stabilizer(c)
+ XXXX_______________
+ _X_X_XX____________
+ __XXXX_XX__________
+ ____X__X__XX_______
+ _________X___X___XX
+ _____XX_XX__XX_____
+ _______XX__XX__XX__
+ __________XX__XX___
+ ____________XX__XX_
+ ZZZZ_______________
+ _Z_Z_ZZ____________
+ __ZZZZ_ZZ__________
+ ____Z__Z__ZZ_______
+ _________Z___Z___ZZ
+ _____ZZ_ZZ__ZZ_____
+ _______ZZ__ZZ__ZZ__
+ __________ZZ__ZZ___
+ ____________ZZ__ZZ_

julia> distance(c, DistanceMIPAlgorithm(solver=HiGHS))
5
```

More information can be seen in [landahl2011color](@cite)

### Fields
    $TYPEDFIELDS
"""
struct Triangular666 <: TriangularCode
    """The distance of the code."""
    d::Int
    function Triangular666(d)
        if d%2!=1
            throw(ArgumentError(THROW_COLOR_CODES_ODD))
        elseif d<3
            throw(ArgumentError(THROW_COLOR_CODES_MIN_DIST))
        end
        return new(d)
    end
end

Triangular488() = Triangular488(3) # smallest d
Triangular666() = Triangular666(3) # smallest d

parity_matrix_x(c::TriangularCode) = _colorcode_get_check_matrix(c)
parity_matrix_z(c::TriangularCode) = _colorcode_get_check_matrix(c)

function _colorcode_get_check_matrix(c::Triangular666)
    n = code_n(c)
    num_checks = (n-1)÷2
    num_layers = (c.d-1)÷2
    checks = zeros(Bool, num_checks, n)

    i = 1
    checks_written = 0
    for layer in 1:num_layers
        # extend half 6-gons from last iteration to full hexagons
        num_6_gons = layer-1
        checks_written -= num_6_gons
        for j in 1:(num_6_gons)
            init_pos = i + 2*(j-1)
            checks[checks_written+1,init_pos+2*(layer-1)+1] = 1
            checks[checks_written+1,init_pos+2*(layer-1)+2] = 1
            checks_written+=1
        end

        # red trapezoid
        checks[checks_written+1,i] = 1
        checks[checks_written+1,i+(layer-1)*2+1] = 1
        checks[checks_written+1,i+2+4*(layer-1)] = 1
        checks[checks_written+1,i+2+4*(layer-1)+1] = 1
        checks_written += 1

        # blue trapezoid
        checks[checks_written+1,i+1+4*(layer-1)] = 1
        checks[checks_written+1,i+1+4*(layer-1)+layer*2] = 1
        checks[checks_written+1,i+5+8*(layer-1)] = 1
        checks[checks_written+1,i+5+8*(layer-1)+1] = 1
        checks_written += 1

        # red hexagons
        for j in 1:(layer-1)
            checks[checks_written+1,i+(j-1)*2+1] = 1
            checks[checks_written+1,i+(j-1)*2+2] = 1
            checks[checks_written+1,i+(j-1)*2+2+2*(layer-1)] = 1
            checks[checks_written+1,i+(j-1)*2+2+2*(layer-1)+1] = 1
            checks[checks_written+1,i+(j-1)*2+2+2*(layer-1)+layer*2] = 1
            checks[checks_written+1,i+(j-1)*2+2+2*(layer-1)+layer*2+1] = 1

            checks_written += 1
        end

        # blue hexagons
        for j in 1:(layer-1)
            init_pos = i+(j-1)*2+(layer-1)*2+1
            checks[checks_written+1, init_pos] = 1
            checks[checks_written+1, init_pos+1] = 1
            checks[checks_written+1, init_pos+2*layer] = 1
            checks[checks_written+1, init_pos+2*layer+1] = 1
            checks[checks_written+1, init_pos+4*layer] = 1
            checks[checks_written+1, init_pos+4*layer+1] = 1

            checks_written += 1
        end

        # green half 6gons
        for j in 0:(layer-1)
            init_pos = i+2*j+2+4*(layer-1)
            checks[checks_written+1,init_pos] = 1
            checks[checks_written+1,init_pos+1] = 1
            checks[checks_written+1,init_pos+2*layer] = 1
            checks[checks_written+1,init_pos+2*layer+1] = 1
            checks_written += 1
        end

        i += 4+6*(layer-1)
    end

    return checks
end

"""Returns the binary matrix defining the X stabilizers for the Triangular488 code. The Z stabilizers are the same.
Based on Fig. 2 of [landahl2011color](@cite)"""
function _colorcode_get_check_matrix(c::Triangular488)
    n = code_n(c)
    num_checks = (n-1)÷2
    num_layers = (c.d-1)÷2
    checks = zeros(Bool, num_checks, n)

    i = 1
    checks_written = 0
    for layer in 1:num_layers
        # Convert half 8-gons from previous layer into full 8-gons
        num_8_gons = layer-1
        checks_written -= num_8_gons
        for j in 1:(num_8_gons)
            checks[checks_written+1,i+2*layer+1+(j-1)*2] = 1
            checks[checks_written+1,i+2*layer+j*2] = 1

            offset = 0
            if layer == num_layers && num_layers%2==0
                offset = 1
            end

            checks[checks_written+1,(2*layer+1)+i+2*layer+1+(j-1)*2 - offset] = 1
            checks[checks_written+1,(2*layer+1)+i+2*layer+j*2 - offset] = 1

            checks_written += 1
        end

        # red 4-gons
        for j in 0:(layer-1)
            checks[checks_written+1,i+j*2] = 1
            checks[checks_written+1,i+1+j*2] = 1
            checks[checks_written+1,i+2*layer+j*2] = 1
            checks[checks_written+1,i+2*layer+1+j*2] = 1
            checks_written += 1
        end

        if layer%2 == 1
            # green half 8-gons on left side
            checks[checks_written+1,i] = 1
            checks[checks_written+1,i+2*layer] = 1
            checks[checks_written+1,i+4*layer] = 1
            checks[checks_written+1,i+4*layer+1] = 1
            checks_written += 1
        else
            # blue half 8-gon on right side
            checks[checks_written+1,i+1+(layer-1)*2] = 1
            checks[checks_written+1,i+2*layer+1+(layer-1)*2] = 1

            # when d=5,9,13,... the final row qubits is indexed slightly differently.
            offset = 0
            if layer == num_layers
                offset = 1
            end

            checks[checks_written+1,(i+4*layer)+layer*2-offset] = 1
            checks[checks_written+1,(i+4*layer)+1+layer*2-offset] = 1
            checks_written += 1
        end

        # blue/green half 8-gons on the bottom
        for j in 0:(layer-1)
            checks[checks_written+1,i+2*layer+j*2] = 1
            checks[checks_written+1,i+2*layer+1+j*2] = 1

            offset = 0
            if layer == num_layers && num_layers%2==0
                offset = 1
            end

            checks[checks_written+1,(2*layer+1)+i+2*layer+j*2 - offset] = 1
            checks[checks_written+1,(2*layer+1)+i+2*layer+1+j*2 - offset] = 1

            checks_written += 1
        end

        i += 4*layer
    end
    return checks
end

"""Returns which qubits each stabilizer touches when given a parity check matrix. Useful for debugging/testing."""
function _colorcode_get_qubit_indices(matrix::Matrix{Bool})
    for j in 1:size(matrix)[1] println(findall(>(0),matrix[j,:])) end
    return
end

# From https://arxiv.org/abs/1108.5738 Fig. 2's caption:
code_n(c::Triangular488) = (c.d^2+2c.d-1)÷2
code_n(c::Triangular666) = (3*c.d^2+1)÷4
code_k(c::TriangularCode) = 1
