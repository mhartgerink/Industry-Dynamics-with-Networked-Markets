#=First we define an updater for each individual firm. It constructs the sequence 
of numbers that corresponds to the size of the set under the possible combinations
=#
function n_calculator(final_array, initi_array)
    for i in 1:(length(final_array)-1)
        final_array[i+1, 1] = initi_array[i+1,1]+final_array[i,1]
    end
    return final_array
end
#= This function calculates the sum of the sequence of the last firm added, and 
    and therefore is equivalent to the value defined in the theorem.
    =#
function final_calculator(n_f, n_m, ω)
    action_size = 2^n_m*(ω+1)
    final_array = ones(action_size)
    initi_array = ones(action_size)
    for i in 1:n_f-1
        final_array = n_calculator(final_array, initi_array)
        initi_array = final_array
    end
    size = sum(initi_array)
    return size
end

#This matrix is for the values listed in the thesis
A = [[2 2 10]
[2 3 10]
[2 4 10]
[2 7 10]
[2 4 5]
[2 4 15]
[2 4 25]
[2 4 30]
[3 2 10]
[3 3 10]
[3 4 10]
[3 7 10]
[3 3 5]
[3 3 15]
[3 3 25]
[3 3 30]
[4 2 10]
[4 3 10]
[4 4 10]
[4 7 10]
[4 3 5]
[4 3 15]
[4 3 25]
[4 3 30]
[5 2 5]
[5 2 10]
[5 3 5]
[5 4 3]
[6 2 5]
[6 3 4]
]
#The function below constructs the table in the thesis
function cardcalc(Matrix)
    Reduction = zeros(size(Matrix, 1))
    Reducedcard = zeros(size(Matrix, 1))
    Nred = zeros(size(Matrix, 1))
    for i in axes(Matrix, 1)
        n_f = Matrix[i, 1]
        n_m = Matrix[i,2]
        maxω = Matrix[i,3]
        Reducedcard[i] = ((1+maxω)*2^n_m)^(n_f)
        Nred[i] = final_calculator(n_f, n_m, maxω)
        Reduction[i] = Reducedcard[i]/Nred[i]
    end
    
    Matrix = hcat(Matrix, Reducedcard, Nred, Reduction)
    tab = Table(Firms = Matrix[:,1], Markets = Matrix[:,2], Maxω = Matrix[:,3], Nonreduced = Matrix[:, 4], Reduced = Matrix[:,5], Reduction = Matrix[:,6])
    CSV.write("Cardinality.csv",  tab; delim = ",")
    return Matrix            
end
cardcalc(A)
