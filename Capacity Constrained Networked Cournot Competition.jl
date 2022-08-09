using Random, Distributions, LinearAlgebra, InvertedIndices, ProgressMeter, NonNegLeastSquares, JLD, SparseArrays
Random.seed!(123)

n_m, maxω = 3,10
α, β, c = rand(Uniform(0.3,1.5), n_m), rand(Uniform(8,10), n_m), rand(Uniform(0.2,0.5))

function importmatrix(n_m, maxω)
    LC = load("LC_matrix_[$n_m, $maxω].jld", "Individual")
    RSS = load("RSS_matrix_[$n_m, $maxω].jld", "States")
    GE = load("Graph_Equiv_[$n_m, $maxω].jld", "Equiv")
    CE = load("Capacity_Equiv_[$n_m, $maxω].jld", "Equiv")
    CD = load("Capgraph_difference_[$n_m, $maxω].jld", "Diff")
    MI = load("Max_investment_[$maxω].jld", "i_level")
    II = load("II_matrix_[$maxω].jld", "Initial")
    CS = load("capacity_states_[$n_m, $maxω].jld", "capacity")
    IL = load("Link_strategy_[$n_m].jld", "pr")
    IGS = load("I_Graph_state_[$n_m, $maxω].jld", "graph")
    OGS = load("O_Graph_state_[$n_m, $maxω].jld", "graph")
    return LC, RSS, GE, CE, CD, MI, II, CS, IL, IGS, OGS
end

LC, RSS, GE, CE, CD, MI, II, CS, IL, IGS, OGS = importmatrix(n_m,maxω)[[1,2,3,4,5,6,7,8,9, 10, 11]]



function U_Profit(RSS::Matrix{Int}, LC::Matrix{Int}, α::Vector{Float64}, β::Vector{Float64}, c::Float64, error, n_m, maxω)
    Profit_vector = zeros(Float64, size(RSS, 2))

    Approx = zeros(Float64, size(RSS, 2), size(RSS, 1))
    for i in axes(RSS, 2)
        
        f_1, f_2 = Vector(LC[RSS[1, i],:]), Vector(LC[RSS[2, i],:])
        e_1, e_2 = f_1[2:end], f_2[2:end]
        c_1, c_2 = e_1.*α, e_2.*α
        cap_1, cap_2 = convert(Float64, LC[RSS[1, i], 1]), convert(Float64, LC[RSS[2,i], 1])

        profit_1 = 0
        if cap_1 ==  0 || sum(e_1) == 0 
            profit_1= 0
        elseif cap_2 == 0 || sum(e_2) == 0
            MR = zeros(Float64, n_m, n_m)
            @inbounds for a in 1:n_m
                MR[a, a] = 2*c_1[a]
            end

            @inbounds for b in 1:n_m
                a = n_m+1-b
                if e_1[a] == 0
                    MR = MR[1:end.!= a, 1:end.!= a]
                else
                    nothing
                end
            end
            
            MR = vcat(MR,(ones(Int, sum(e_1)).*(-c))')
            MR = hcat(MR, ones(Int, size(MR,1)))
            β_new = vcat(e_1.*β)
            s_b = size(β_new, 1)
            #reduce the vector
            @inbounds for a in axes(β_new,1)
                b = s_b+1-a
                if β_new[b] == 0
                    β_new = β_new[1:end.!=b]
                else 
                    nothing
                end
            end
            β_new = vcat(β_new, [0])
            q = nonneg_lsq(MR,β_new)
            
            #Calculate the values for each firm
            
            quantity_1 = zeros(n_m)            
            quantity_1[findall(x->x==1, e_1)] = q[1:sum(e_1)]
            quantity_2 = zeros(n_m)

        else
                     
            #First calculate the unconstrained equilibrium
            MR = zeros(Float64, 2*n_m, 2*n_m)
            #Construct the MR-Matrix
            @inbounds for a in 1:2*n_m
                if a <= n_m
                    MR[a, a] = 2*c_1[a]
                    MR[a, n_m+a] = c_2[a]
                else
                    MR[a, a] = 2*c_2[a-n_m]
                    MR[a, a-n_m] = c_1[a-n_m]
                end
            end
            #Reducing the columns and rows to ensure non-singularity of the matrix
            @inbounds for b in 1:n_m
                a = n_m+1-b
                if e_2[a] == 0
                    MR = MR[1:end.!= a+n_m, 1:end.!= a+n_m]
                else
                    nothing
                end
            end
            @inbounds for b in 1:n_m
                a = n_m+1-b
                if e_1[a] == 0
                    MR = MR[1:end.!= a, 1:end.!= a]
                else
                    nothing
                end
            end

            # We now add the costs and the MR=MC conditions#
            MR = vcat(hcat(MR, vcat(hcat(ones(Int, sum(e_1)), zeros(Int, sum(e_1))), hcat(zeros(Int, sum(e_2)), ones(Int, sum(e_2))))), hcat(vcat(hcat(ones(Int, sum(e_1))', zeros(Int, sum(e_2))'), hcat(zeros(Int, sum(e_1))', ones(Int, sum(e_2))')).*(-c),I(2)*1))
            #construct the β
            
            β_new = vcat(e_1.*β, e_2.*β)
            s_b = size(β_new, 1)
            #reduce the vector
            @inbounds for a in axes(β_new,1)
                b = s_b+1-a
                if β_new[b] == 0
                    β_new = β_new[1:end.!=b]
                else 
                    nothing
                end
            end
            β_new = vcat(β_new, [0,0])
            
            #Calculate the values for each firm
            q = nonneg_lsq(MR,β_new)
            quantity_1 = zeros(n_m)            
            quantity_1[findall(x->x==1, e_1)] = q[1:sum(e_1)]
            quantity_2 = zeros(n_m)
            quantity_2[findall(x->x==1, e_2)] = q[sum(e_1)+1:sum(e_1)+sum(e_2)]

            
        end

        if cap_1 == 0 || sum(e_1)==0
            nothing
        else
            if sum(quantity_1) > cap_1+error || sum(quantity_2) > cap_2+error 
                
                #Calculate the correct quantities for each constraint
                if cap_2 == 0 || sum(e_2) == 0
                    MR[end, end] = 0
                    MR = hcat(MR, vcat(zeros(size(MR, 1)-1), 1))
                    MR = vcat(MR, hcat(ones(Int, sum(e_1))', zeros(Int, 2)'))
                    β_new = vcat(β_new, cap_1)
                    q = nonneg_lsq(MR, β_new)
                    quantity_1 = zeros(n_m)            
                    quantity_1[findall(x->x==1, e_1)] = q[1:sum(e_1)]
                    if sum(quantity_1)>cap_1+error
                        print(i, "  No q2, error")
                    end
                end

                #CHECKING CONSTRAINTS
                if sum(quantity_1) > cap_1+error && sum(quantity_2) <= cap_2+error 

                    quantity_1, quantity_2 = profit_constraint(MR, e_1, e_2, β_new, cap_1, 0, error)[[1,2]]

                elseif sum(quantity_1) <= cap_1+error && sum(quantity_2) > cap_2+error 

                    quantity_1, quantity_2 = profit_constraint(MR, e_1, e_2, β_new, 0, cap_2, error)[[1,2]]

                elseif sum(quantity_1) > cap_1+error && sum(quantity_2) > cap_2+error 

                    quantity_1, quantity_2 = profit_constraint(MR, e_1, e_2, β_new, cap_1, cap_2, error)[[1,2]]
                    

                end
            end
             
            
            quantity = quantity_1+quantity_2
            price = zeros(Float64, n_m)
            @inbounds for a in 1:n_m
                price[a] = β_new[a]-α[a]*quantity[a]
            end
            profit_1 = 0
            @inbounds for a in 1:n_m
                profit_1 = profit_1+price[a]*quantity_1[a]
            end
            if sum(quantity_1) > cap_1+error && sum(quantity_2) > cap_2+error 
                Approx[i, :] =[sum(quantity_1)-cap_1, sum(quantity_2)]       
            elseif sum(quantity_1) > cap_1+error && sum(quantity_2) <= cap_2+error 
                Approx[i, 1] = sum(quantity_1)-cap_1   
            elseif sum(quantity_1) <= cap_1+error && sum(quantity_2) > cap_2+error 
                Approx[i, 2] = sum(quantity_2)            
            end
            
        end
        if i == 130
            println(quantity_1, cap_1, quantity_2, cap_2)
        end
        if profit_1 > 0
            Profit_vector[i] = profit_1
        else
            Profit_vector[i] = 0
        end

    end
    error = size(findall(x->x!=0, Approx), 1)/(size(Approx,1)*size(Approx,2))
    sd_array = []
    for i in axes(Approx, 1), j in axes(Approx,2)
        if Approx[i,j] == 0
            nothing
        else
            sd_array = vcat(sd_array, Approx[i,j])
        end
    end
    av = sum(sd_array)/size(sd_array,1)
    sd = std(sd_array)
    #save("Profit_[3,2]_se.jld", "prof", Profit_vector)
    #save("Error_[3,2]_se.jld", "err", Approx)

    return nothing
end

function profit_constraint(A, b_1, b_2, β, cap_1, cap_2, error)
    
    e_1 = copy(b_1)
    e_2 = copy(b_2)
    MR = copy(A)
    if cap_1 > 0 && cap_2 > 0

        MR = MR[1:sum(e_1)+sum(e_2)+2, 1:sum(e_1)+sum(e_2)+2]
        β_new = β[1:sum(e_1)+sum(e_2)+2]
        #Make sure that MC and MR are allowed to differ
        MR[sum(e_1)+sum(e_2)+1, sum(e_1)+sum(e_2)+1] = 0
        MR[sum(e_1)+sum(e_2)+2, sum(e_1)+sum(e_2)+2] = 0
        #Impose constraints
        con_1 = hcat(ones(Int, sum(e_1))', zeros(Int, size(MR,2)-sum(e_1))')
        con_2 = hcat(zeros(Int, sum(e_1))', ones(Int, sum(e_2))', zeros(Int, 2)')

        vcon_1 = vcat(zeros(Int, sum(e_1+e_2)), 1, zeros(Int, 3))
        vcon_2 = vcat(zeros(Int, sum(e_1+e_2)+1), 1, zeros(Int, 2))
        #Add the constraints
        MR = vcat(MR, con_1, con_2)
        MR = hcat(MR, vcon_1, vcon_2)
        #Add the constraints to the β_new
        β_new = vcat(β_new, cap_1, cap_2)
        
        q = nonneg_lsq(MR,β_new)
        quantity_1 = zeros(n_m)            
        quantity_1[findall(x->x==1, e_1)] = q[1:sum(e_1)]
        quantity_2 = zeros(n_m)
        quantity_2[findall(x->x==1, e_2)] = q[sum(e_1)+1:sum(e_1)+sum(e_2)]
        
        


        if sum(q[1:sum(e_1)]) > cap_1 + error || sum(q[sum(e_1)+1:sum(e_1)+sum(e_2)]) > cap_2 + error
            MR = MR[1:sum(e_1)+sum(e_2)+2, 1:sum(e_1)+sum(e_2)+2]
            β_new = β[1:sum(e_1)+sum(e_2)+2]
            #Make sure that MC and MR are allowed to differ
            
            MR[sum(e_1)+sum(e_2)+2, sum(e_1)+sum(e_2)+2] = 0
            #Impose constraints
            
            con_2 = hcat(zeros(Int, sum(e_1))', ones(Int, sum(e_2))', zeros(Int, 2)')

            vcon_2 = vcat(zeros(Int, sum(e_1+e_2)+1), 1, zeros(Int, 1))
            #Add the constraints
            MR = vcat(MR, con_2)
            MR = hcat(MR, vcon_2)
            #Add the constraints to the β_new
            β_new = vcat(β_new, cap_2)
            
            q = nonneg_lsq(MR,β_new)
            quantity_1 = zeros(n_m)            
            quantity_1[findall(x->x==1, e_1)] = q[1:sum(e_1)]
            quantity_2 = zeros(n_m)
            quantity_2[findall(x->x==1, e_2)] = q[sum(e_1)+1:sum(e_1)+sum(e_2)]


        end     

    elseif cap_1 > 0  && cap_2 == 0
        MR = MR[1:sum(e_1)+sum(e_2)+2, 1:sum(e_1)+sum(e_2)+2]
        β_new = β[1:sum(e_1)+sum(e_2)+2]
        #Make sure that MC and MR are allowed to differ
        MR[sum(e_1)+sum(e_2)+1, sum(e_1)+sum(e_2)+1] = 0
        
        #Impose constraints
        con_1 = hcat(ones(Int, sum(e_1))', zeros(Int, size(MR,2)-sum(e_1))')
        
        vcon_1 = vcat(zeros(Int, sum(e_1+e_2)), 1, zeros(Int, 2))

        #Add the constraints
        MR = vcat(MR, con_1)
        MR = hcat(MR, vcon_1)
        #Add the constraints to the β_new
        β_new = vcat(β_new, cap_1)
        
        q = nonneg_lsq(MR,β_new)
        quantity_1 = zeros(n_m)            
        quantity_1[findall(x->x==1, e_1)] = q[1:sum(e_1)]
        quantity_2 = zeros(n_m)
        quantity_2[findall(x->x==1, e_2)] = q[sum(e_1)+1:sum(e_1)+sum(e_2)]
        if sum(q[1:sum(e_1)]) > cap_1 + error || sum(q[sum(e_1)+1:sum(e_1)+sum(e_2)]) > cap_2 + error
            MR = MR[1:sum(e_1)+sum(e_2)+2, 1:sum(e_1)+sum(e_2)+2]
            β_new = β[1:sum(e_1)+sum(e_2)+2]
            #Make sure that MC and MR are allowed to differ
            
            MR[sum(e_1)+sum(e_2)+2, sum(e_1)+sum(e_2)+2] = 0
            #Impose constraints
            
            con_2 = hcat(zeros(Int, sum(e_1))', ones(Int, sum(e_2))', zeros(Int, 2)')

            vcon_2 = vcat(zeros(Int, sum(e_1+e_2)+1), 1, zeros(Int, 1))
            #Add the constraints
            MR = vcat(MR, con_2)
            MR = hcat(MR, vcon_2)
            #Add the constraints to the β_new
            β_new = vcat(β_new, cap_2)
            
            q = nonneg_lsq(MR,β_new)
            quantity_1 = zeros(n_m)            
            quantity_1[findall(x->x==1, e_1)] = q[1:sum(e_1)]
            quantity_2 = zeros(n_m)
            quantity_2[findall(x->x==1, e_2)] = q[sum(e_1)+1:sum(e_1)+sum(e_2)]


        end     

    elseif cap_1 ==0 && cap_2 > 0
        
        MR = MR[1:sum(e_1)+sum(e_2)+2, 1:sum(e_1)+sum(e_2)+2]
        β_new = β[1:sum(e_1)+sum(e_2)+2]
        #Make sure that MC and MR are allowed to differ
        
        MR[sum(e_1)+sum(e_2)+2, sum(e_1)+sum(e_2)+2] = 0
        #Impose constraints
        
        con_2 = hcat(zeros(Int, sum(e_1))', ones(Int, sum(e_2))', zeros(Int, 2)')

        vcon_2 = vcat(zeros(Int, sum(e_1+e_2)+1), 1, zeros(Int, 1))
        #Add the constraints
        MR = vcat(MR, con_2)
        MR = hcat(MR, vcon_2)
        #Add the constraints to the β_new
        β_new = vcat(β_new, cap_2)
        
        q = nonneg_lsq(MR,β_new)
        quantity_1 = zeros(n_m)            
        quantity_1[findall(x->x==1, e_1)] = q[1:sum(e_1)]
        quantity_2 = zeros(n_m)
        quantity_2[findall(x->x==1, e_2)] = q[sum(e_1)+1:sum(e_1)+sum(e_2)]
          
        
    end

    
    return quantity_1, quantity_2
end


@time begin
U_Profit(RSS, LC, α, β, c, 0.005, n_m, maxω)
end


#error = load("Error_[3,2]_se.jld", "err")



