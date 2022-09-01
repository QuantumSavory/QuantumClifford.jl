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
            rowcounter = 1
            endevenrow = j-1
            counterqineven = 2

            for location in range(1,c.n*c.n)
                if (rowcounter % j == 0) && ((counterqineven - j/2) % j != 0) 
                    if endevenrow !=0 && endevenrow % (j+1) == 0
                        push!(evenrows,location)
                        counterqineven += 1
                        endevenrow += 1
                        rowcounter += 1
                    else
                        push!(evenrows,location)
                        endevenrow += 1
                    end
                elseif location < j*rowcounter && location <= c.n
                    push!(oddrows,location)
                elseif location == j*rowcounter && location <= c.n #issue 1st 
                    push!(oddrows,location)
                    rowcounter += 1
                    counterqineven = 0
                end #if
            end #for
            @show evenrows
            @show oddrows
            @show rowcounter
            @show counterqineven
            @show endevenrow

            for location in range(1,c.n*c.n)
                #Z gates
                #within the rox bounds
                #last location within length bounds
                #want to make sure w don't take midle row values  
                if z_n <= (c.n) && (location+2j+1) <= (c.n) && location in oddrows
                    @show location
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
                    x_locations[x_n,location] = 1
                    x_locations[x_n,location+1] = 1
                    x_locations[x_n,location+1+j] = 1
                    x_locations[x_n,location+2j+1] = 1                   
                elseif (location+1) <= (c.n) && x_n <= (c.n)
                    x_locations[x_n,location] = 1
                    x_locations[x_n,location+1] = 1
                elseif (location+1+j) <= (c.n) && x_n <= (c.n)
                    x_locations[x_n,location] = 1
                    x_locations[x_n,location+1+j] = 1
                #=
                elseif isdefined(location,x_n) && isdefined(location+1+j,x_n) 
                    x_locations[location,x_n] = 1
                    x_locations[location+1+j,x_n] = 1
                =#
                end #pretty sure there are more combos

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