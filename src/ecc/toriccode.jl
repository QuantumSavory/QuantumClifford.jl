import .ECC

using LinearAlgebra
using .ECC

struct Toric <: AbstractECC 
    n
end 

code_n(c::Toric) = c.n

#Parity checks ----------------------------------

function checks(c::Toric) 
    #not all n can make square lattices
    available_n = [4,7,12,17,24,31] #these values only go up to 9 sided grids - working on a generation function
    
    for i in available_n 
        if c.n == i
            for j in range(0,1,1000)
                if i <= j^2
                    grid = zeros(Int8,j, j)
                    z_locations = zeros(j, j)
                    x_locations = zeros(j, j)
                    z_n = 0
                    x_n = 0

                    for i in grid
                        #Z gates
                        if i == b && b < (c.n-2n) && b_row < rows
                            z_locations[i,z_n] = 1
                            z_locations[i+n,z_n] = 1
                            z_locations[i+n+1,z_n] = 1
                            z_locations[i+n+1+n,z_n] = 1
                            #TODO not count edges
                            #TODO not count unused qubits
                        end

                        #X gates
                        if isdefined(grid[i+n],0) && isdefined(grid[i+(2*c.n)+1],0) && isdefined(grid[i+(3*c.n)+1],0)
                            x_locations[i,x_n] = 1
                            x_locations[i+n,x_n] = 1
                            x_locations[i+(2*n)+1,x_n] = 1
                            x_locations[i+(3*n)+1,x_n] = 1                   
                        elseif isdefined(grid[i+n],0) 
                            x_locations[i,n] = 1
                            x_locations[i+n,x_n] = 1
                        elseif isdefined(grid[i+(2*n)+1],0) 
                            x_locations[i,x_n] = 1
                            x_locations[i+(2*n)+1,x_n] = 1
                        elseif isdefined(grid[i+(3*n)+1],0) 
                            x_locations[i,x_n] = 1
                            x_locations[i+(3*n)+1,x_n] = 1
                        end

                        z_n += 1 #next z
                        x_n += 1 #next x

                    end #for 

                    #making X & Z into bool
                    z_locations = !=(0).(Z)
                    x_locations = !=(0).(X)

                    return Stabilizer(x_locations,z_locations)

                end #if
            end #for
        end #if
    end #for

end #function

parity_checks(c::Toric) = checks(c)
#-----------------------------------------------------

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------

distance(c::Toric) = sqrt(c.n)