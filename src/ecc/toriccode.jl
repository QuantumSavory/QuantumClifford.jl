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
    available_n = [4,12,24] #these values only go up to 9 sided grids - working on a generation function #not count unused qubits
    for i in available_n
        if c.n == i
            if i == 4
                j = 1
            else
                j = i/4 -1
            end

            j = convert(Int8, j)
            z_locations = zeros(c.n, c.n)
            x_locations = zeros(c.n, c.n)
            z_n = 1
            x_n = 1

            #Row classification
            oddrows = []
            evenrows = []
            for i in range(1,c.n+1)
                if i % ((2j)+(j/2)) == 0 && i != 0
                    push!(evenrows,i-4)
                    push!(evenrows,i-3)
                    push!(oddrows,i-2)    
                    push!(oddrows,i-1)  
                    push!(oddrows,i) 
                end               
            end
            
            for location in range(1,c.n)
                #Z gates
                #within the rox bounds
                #last location within length bounds
                #want to make sure w don't take midle row values  
                if z_n <= (c.n) && (location+2j+1) <= (c.n) && location in evenrows
                    z_locations[z_n,location] = 1
                    z_locations[z_n,location+j] = 1
                    z_locations[z_n,location+j+1] = 1
                    z_locations[z_n,location+2j+1] = 1
                #=
                elseif z_locations[location] !=1 && z_n < (c.n) && (location) < (c.n) 
                    z_locations[location,z_n] = 0
                =#
                end

                #X gates -not running as expected
                if (location+1+2j) <= (c.n) && x_n <= (c.n)
                    for e in evenrows
                        for o in oddrows
                            for l in evenrows
                                for p in oddrows
                                    
                                    if (location == e && location + j == o) || (location == o && location + j == e) # top and left
                                        
                                        if ((location == e && location + j +1 != o) || (location == o && location + j +1 != e)) && ((location == e && location + j != o) || (location == o && location + j != e))#making sure we are not at a border
                                            
                                            if (location == o && location+j == e && location+1+j == l && x_n,location+2j+1 == p) || (location == e && location+j == o && location+1+j == p && x_n,location+2j+1 == l)
                                                # 4 sides
                                                x_locations[x_n,location] = 1
                                                x_locations[x_n,location+j] = 1
                                                x_locations[x_n,location+1+j] = 1
                                                x_locations[x_n,location+2j+1] = 1
                                            
                                            elseif (location == o && location+j == e && x_n,location+2j+1 == p) || (location == e && location+j == o && x_n,location+2j+1 == l) #top left
                                                x_locations[x_n,location] = 1
                                                x_locations[x_n,location+j] = 1
                                                x_locations[x_n,location+2j+1] = 1
                                            
                                            elseif (location == o && location+j+1 == e && x_n,location+2j+1 == p) || (location == e && location+j+1 == o && x_n,location+2j+1 == l) #top right
                                                x_locations[x_n,location] = 1
                                                x_locations[x_n,location+1+j] = 1
                                                x_locations[x_n,location+2j+1] = 1
                                            end
                                        end 

                                    elseif (location+j) <= (c.n) && x_n <= (c.n) #right bottom
                                        if (location == e && location + j == o) || (location == o && location + j == e)
                                            x_locations[x_n,location] = 1
                                            x_locations[x_n,locationj] = 1
                                        end

                                    elseif (location+j+j/2) <= (c.n) && x_n <= (c.n) #left bottom
                                        if (location == e && location + j +j/2 == o) || (location == o && location +j +j/2 == e)
                                            x_locations[x_n,location] = 1
                                            x_locations[x_n,location+j+j/2] = 1
                                        end
                                    end
                                end
                            end
                        end
                    end
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

end #function

parity_checks(c::Toric) = checks(c)
#-----------------------------------------------------

parity_matrix(c::Toric) = stab_to_gf2(parity_checks(c))

#Encoding circuit ----------------------------------

encoding_circuit(c::Toric) = [] #TODO
#-----------------------------------------------------

distance(c::Toric) = sqrt(c.n)

code_k(c::Toric) = 1