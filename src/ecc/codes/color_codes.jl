abstract type ColorCode <: AbstractECC end

"""Planar color codes that encode a single logical qubit."""
abstract type TriangularCode <: ColorCode end

"""Triangular code following the 4.8.8 tiling. Constructor take a distance d as input."""
struct Triangular4_8_8 <: TriangularCode
    d::Int
    function Triangular4_8_8(d)
        if d%2!=1
            throw(DomainError("only odd distance trianglular color codes are allowed.\nRefer to https://arxiv.org/abs/1108.5738"))
        end
        return new(d)
    end
end

"""Triangular code following the 6.6.6 tiling. Constructor take a distance d as input."""
struct Triangular6_6_6 <: TriangularCode
    d::Int
    function Triangular6_6_6(d)
        if d%2!=1
            throw(DomainError("only odd distance trianglular color codes are allowed.\nRefer to https://arxiv.org/abs/1108.5738"))
        end
        return new(d)
    end
end

Triangular4_8_8() = Triangular4_8_8(3) # smallest d
Triangular6_6_6() = Triangular6_6_6(3) # smallest d

function parity_checks(c::TriangularCode)
    matrix = get_check_matrix(c)

    num_checks, qubits = size(matrix)
    return Stabilizer(vcat(matrix,zeros(Bool,num_checks,qubits)), vcat(zeros(Bool,num_checks,qubits), matrix))
end

function get_check_matrix(c::Triangular6_6_6) 
    # TODO make code look a bit nicer by copying the style of the "init_pos" code blocks
    n = code_n(c)
    num_checks = (n-1)/2 |> Int 
    num_layers = (c.d-1)/2 |> Int
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

"""Returns the binary matrix defining the X stabilizers for the Triangular4_8_8 code. The Z stabilizers are the same.
Based on Fig. 2 of https://arxiv.org/abs/1108.5738"""
function get_check_matrix(c::Triangular4_8_8)
    n = code_n(c)
    num_checks = (n-1)/2 |> Int 
    num_layers = (c.d-1)/2 |> Int
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
function get_qubit_indices(matrix::Matrix{Bool})
    for j in 1:size(matrix)[1] println(findall(>(0),matrix[j,:])) end
    return
end

function iscss(::ColorCode)
    return true
end

"""From https://arxiv.org/abs/1108.5738 Fig. 2's caption:"""
code_n(c::Triangular4_8_8) = 0.5*c.d^2+c.d-0.5 |> Int
code_n(c::Triangular6_6_6) = 0.75*c.d^2+.25 |> Int
code_k(c::TriangularCode) = 1