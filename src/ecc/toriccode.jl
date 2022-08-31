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
    available_n = [4,7,12,17,24,31] #these values only go up to 9 sided grids - working on a generation function #not count unused qubits
    
    for i in available_n 
        if c.n == i
            for j in range(0,1000)
                if i <= j^2
                    grid = zeros(Int8,j, j)
                    z_locations = zeros(j, j)
                    x_locations = zeros(j, j)
                    z_n = 1
                    x_n = 1

                    for location in range(1,j*j)
                        #Z gates
                        if (location + 1)%(j) !=0 && (location - 1)%(j) !=0 && z_n < (j) && (location%j) !=0 #i is not the last layer,is not at the end of a row, not outside of bounds
                            z_locations[location,z_n] = 1
                            z_locations[location+1,z_n] = 1
                            z_locations[location,z_n+1] = 1
                            z_locations[location+1,z_n+1] = 1
                            println("Z has run")
                        end

                        #X gates -not running as expected
                        if isdefined(location,x_n+1) && isdefined(location+1,x_n+1) && isdefined(location+1,x_n+2)
                            x_locations[location,x_n] = 1
                            x_locations[location,x_n+1] = 1
                            x_locations[location+1,x_n+1] = 1
                            x_locations[location+1,x_n+2] = 1                   
                            println("X has run")
                        elseif isdefined(location,x_n+1) 
                            x_locations[location,x_n] = 1
                            x_locations[location,x_n+1] = 1
                        elseif isdefined(location+1,x_n+1) 
                            x_locations[location,x_n] = 1
                            x_locations[location+1,x_n+1] = 1
                        elseif isdefined(location+1,x_n+2) 
                            x_locations[location,x_n] = 1
                            x_locations[location+1,x_n+2] = 1
                        end

                        z_n += 1 #next z
                        x_n += 1 #next x

                    end #for 

                    #making X & Z into bool
                    Z = !=(0).(z_locations)
                    X = !=(0).(x_locations)

                    return Stabilizer(X,Z)

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
