using JLD, Random, ProgressMeter, Distributions, SIMD, SparseArrays, QuadGK

Random.seed!(123)
# State Space

n_m = 3
maxω = 10

function importmatrix(n_m, maxω)
    LC = load("LC_matrix_[$n_m, $maxω].jld", "Individual")
    RSS = load("RSS_matrix_[$n_m, $maxω].jld", "States")
    GE = load("Graph_Equiv_[$n_m, $maxω].jld", "Equiv")
    CE = load("Capacity_Equiv_[$n_m, $maxω].jld", "Equiv")
    CD = load("Capgraph_difference_[$n_m, $maxω].jld", "Diff")
    #MI = load("Max_investment_[$maxω]_last.jld", "i_level")
    #II = load("II_matrix_[$maxω]_last.jld", "Initial")
    #CS = load("capacity_states_[$n_m, $maxω].jld", "capacity")
    #IL = load("Link_strategy_[$n_m]_last.jld", "pr")
    #IGS = load("I_Graph_state_[$n_m, $maxω].jld", "graph")
    #OGS = load("O_Graph_state_[$n_m, $maxω].jld", "graph")
    #MS = load("Merger_states_[$n_m, $maxω].jld", "states")
    #GEn = load("Entry_Diff_[$n_m]_last.jld","diff")
    #GEx = load("Exit_Diff_[$n_m]_last.jld", "diff")
    #CCNCC = load("Profit_[$n_m,$maxω]_se.jld", "prof")
    return LC, RSS, GE, CE, CD, MI, II, CS, IL, IGS, OGS, MS, GEn, GEx, CCNCC
end

#SET UP *** VARIABLES ***


#INVESTMENT STAGE
β, cost_investment, exponent, minη, maxη = 0.925, 1.5, 2.0, -2, 2

#*** GENERAL MATRICES ***
LC, RSS, GE, CE, CD, MI, II, CS, IL, IGS, OGS, MS, GEn, GEx, CCNCC = importmatrix(n_m,maxω)[[1,2,3,4,5,6,7,8,9,10,11,12,13,14, 15]]

#INTIAL MATRICES
V = CCNCC./(1-β)

#=*** GENERAL FUNCTIONS/DISTRIBUTIONS ***
    convolution of two binomial distributions
=#

function Δω(Δω::Int, x::Int, ω::Int)
    p::Float64 = 0
    @inbounds for i in Δω:min(x, (ω+Δω))::Int
        p = p+pdf(Binomial(x, 0.8), i)*pdf(Binomial(ω, 0.3), i - Δω)::Float64
    end
    return p
end


#INVESTMENT STAGE

function I_policy(c::Float64, exponent::Float64, β::Float64, II::Matrix{Float64}, CD::Matrix{Int}, GE::Matrix{Int}, RSS::Matrix{Int}, V::Vector{Float64}, CS::Vector{Any}, MI::Vector{Any})
    MATRIX = zeros(Float64, size(RSS, 2), size(II, 2))
    Value = zeros(Float64, size(RSS, 2))
    
    for i in axes(RSS, 2)
        EV_matrix = zeros(Float64, MI[i]+1) 
        graph_equivalents = GE[Tuple(findfirst(x->x==i, GE))[1],:]
        ω_differences = CD[(i-1)*2+1:i*2,:]

        if RSS[1,i] == RSS[2, i]
            #If symmetric, firms use the same investment strategies
            opponent_strategy = @view(II[i,:])
            if MI[i] == 0
                MATRIX[i, 1] == 1.0
                EV = 0.0
                for j in axes(graph_equivalents,1), k in axes(opponent_strategy,1)
                    EV = EV + V[graph_equivalents[j]]*Δω(ω_differences[1,j], 0, CS[i])*(Δω(ω_differences[2,j], k-1, CS[i])*opponent_strategy[k])
                end
                EV_matrix[1] = β*EV #cost of investment equals 0
            else
                @inbounds for a in 1:MI[i]+1::Int
                EV = EV_matrix[a]::Float64
                @inbounds for j in axes(graph_equivalents,1), k in axes(opponent_strategy, 1)
                    EV = EV + V[graph_equivalents[j]]*Δω(ω_differences[1,j], a-1, CS[i])*(Δω(ω_differences[2,j], k-1, CS[i])*opponent_strategy[k])
                end
                EV_matrix[a] = -c*(a-1)^exponent+β*EV::Float64
                end
            end
        else
            #If asymmetric, firms use different investment strategies
            if RSS[1, i] == RSS[2, i+1]
                #Its counterpart is at i+1, by construction of RSS
                opponent_strategy = @view(II[i+1,:])
                if MI[i] == 0
                    MATRIX[i, 1] == 1.0
                    EV = 0.0
                    for j in axes(graph_equivalents,1), k in axes(opponent_strategy,1)
                        EV = EV + V[graph_equivalents[j]]*Δω(ω_differences[1,j], 0, CS[i])*(Δω(ω_differences[2,j], k-1, CS[i])*opponent_strategy[k])
                    end
                    EV_matrix[1] = β*EV #cost of investment equals 0
                else
                    @inbounds for a in 1:MI[i]+1::Int
                        EV = EV_matrix[a]::Float64
                        @inbounds for j in axes(graph_equivalents,1), k in axes(opponent_strategy, 1)
                            EV = EV + V[graph_equivalents[j]]*Δω(ω_differences[1,j], a-1, CS[i])*(Δω(ω_differences[2,j], k-1, CS[i])*opponent_strategy[k])
                        end
                        EV_matrix[a] = -c*(a-1)^exponent+β*EV::Float64
                    end
                end
            else
                #Its counterpart is at i-1, by construction of RSS
                opponent_strategy = @view(II[i-1,:])
                if MI[i] == 0
                    MATRIX[i, 1] == 1.0
                    EV = 0.0
                    for j in axes(graph_equivalents,1), k in axes(opponent_strategy,1)
                        EV = EV + V[graph_equivalents[j]]*Δω(ω_differences[1,j], 0, CS[i])*(Δω(ω_differences[2,j], k-1, CS[i])*opponent_strategy[k])
                    end
                    EV_matrix[1] = β*EV #cost of investment equals 0
                else
                    
                    @inbounds for a in 1:MI[i]+1::Int
                        
                        EV = EV_matrix[a]
                        @inbounds for j in axes(graph_equivalents,1), k in axes(opponent_strategy, 1)
                            EV = EV + V[graph_equivalents[j]]*Δω(ω_differences[1,j], a-1, CS[i])*(Δω(ω_differences[2,j], k-1, CS[i])*opponent_strategy[k])
                            #p[a] = p[a] + Δω(ω_differences[1,j], a-1, CS[i])*(Δω(ω_differences[2,j], k-1, CS[i])*opponent_strategy[k])
                            #println(p)
                        end
                        
                        EV_matrix[a] = -c*(a-1)^exponent+β*EV::Float64
                    end
                    
                end
            end
        end
        if MI[i] == 0
            Value[i] = sum(EV_matrix)
            MATRIX[i, 1] = 1.0
        else
            MAX = zeros(Float64, size(EV_matrix, 1))
            @simd for j in axes(EV_matrix,1)
                MAX[j] = maximum(EV_matrix[1:end.!=j])
            end
            @inbounds for j in axes(EV_matrix,1), k in axes(EV_matrix,1)
                MATRIX[i, j] = MATRIX[i, j]+(1-cdf(Uniform(-5,5), (MAX[j]-EV_matrix[j])))*1/(size(EV_matrix,1)^2)
                #The probability of obtaining the shock is squared, as this loop also ranges over all values of k, which entails that the probability of obtaining the shock is higher as well
                MATRIX[i, k] = MATRIX[i, k] + (1/size(EV_matrix,1))*cdf(Uniform(-5,5), EV_matrix[k]-EV_matrix[j])*(EV_matrix[k] in MAX[j])
            end
            Value[i] = sum(EV_matrix.*MATRIX[i, 1:size(EV_matrix,1)])/(MI[i]+1)
        end 
    end
    
    return Value, MATRIX
end

#MERGER STAGE

function merger(LC, RSS, maxω, n_m)
    merger_states = zeros(Int, size(RSS, 2))
    b = findall(x->x==size(LC,1), RSS[2,:])
    for i in axes(RSS, 2)
        
        new_state = LC[RSS[1, i], :]+LC[RSS[2,i], :]
        if new_state[1] > maxω
            new_state[1] = maxω
        else
            nothing
        end

        for a in 2:size(new_state,1)::Int
            if new_state[a] == 2
                new_state[a] = 1
            else
                nothing
            end
        end
        new_LC_index = 0
        for a in axes(LC,1)
            if LC[a,:] == new_state
                new_LC_index = a
            else
                nothing
            end
        end

        for a in b
            if new_LC_index == RSS[1, a]
                merger_states[i] = a
                
            else
                nothing
            end
        end
        
    end
    save("Merger_states_[$n_m, $maxω].jld", "states", merger_states)
    return merger_states

end


function M_policy(RSS::Matrix{Int}, V::Vector{Float64}, MS::Vector{Int}, minsyn, maxsyn)
    EV = zeros(Float64, size(RSS, 2))
    Policy = zeros(Float64, size(RSS, 2))
    @inbounds for i in axes(RSS, 2)
        reservation_value = V[i]
        merger_value = V[MS[i]]
        #symmetric
        if RSS[1,i] == RSS[2, i]
            opp_reservation_value = V[i]
            Pr_merge = 1-cdf(Uniform(minsyn, maxsyn), opp_reservation_value+reservation_value-merger_value)
            EV[i] = 1/2*reservation_value+1/2*(Pr_merge*(merger_value-opp_reservation_value)+(1-Pr_merge)*(reservation_value))
            Policy[i] = 1-cdf(Uniform(0,maxsyn), (reservation_value-merger_value+opp_reservation_value))
        else
            if RSS[1, i] == RSS[2, i+1]
                opp_reservation_value = V[i]
                Pr_merge = 1-cdf(Uniform(minsyn, maxsyn),opp_reservation_value+reservation_value-merger_value)
                EV[i] = 1/2*reservation_value+1/2*(Pr_merge*(merger_value-opp_reservation_value)+(1-Pr_merge)*(reservation_value))
                Policy[i] = 1-cdf(Uniform(minsyn,maxsyn), (reservation_value-merger_value+opp_reservation_value))
            else
                opp_reservation_value = V[i]
                Pr_merge = 1-cdf(Uniform(minsyn, maxsyn),opp_reservation_value+reservation_value-merger_value)
                EV[i] = 1/2*reservation_value+1/2*(Pr_merge*(merger_value-opp_reservation_value)+(1-Pr_merge)*(reservation_value))
                Policy[i] = 1-cdf(Uniform(minsyn,maxsyn), (reservation_value-merger_value+opp_reservation_value))
            end
        end
    end

    return EV, Policy
end

#ENTRY/EXIT STAGE

#Irwin-Hall distributions
function CDF_convolution_IW(x::Float64)
    f(a) = pdf(Uniform(0,1), a)*(1-cdf(Uniform(0,1), a-x))
    b = quadgk(a->  f(a), 0, 1)[1]
    return b
end

function EE_policy(RSS::Matrix{Int}, CE::Matrix{Int}, V::Vector{Float64}, IL::Matrix{Float64}, n_m, OGS, IGS, GEn, GEx, ϕ, ψ)
    value_vector = zeros(Float64, size(RSS, 2))
    policy_matrix = zeros(Float64, size(RSS, 2), 2^n_m)
    
    for i in axes(RSS, 2) 
        EV_matrix = zeros(Float64, 2^n_m)
        capacity_equivalents = CE[Tuple(findfirst(x->x==i, CE))[1],:]
        #CREATE THE MATRIX FOR THE EXPECTED VALUE FUNCTIONS FOR EACH CHOICE
        if RSS[1,i] == RSS[2, i]
            opponent_strategy = IL[i,:]
            for a in axes(capacity_equivalents, 1)
                α = capacity_equivalents[a]
                EV_matrix[IGS[α]]= EV_matrix[IGS[α]]+V[α]*opponent_strategy[OGS[α]]

            end
        else
        #If asymmetric, firms use different strategies
            if RSS[1, i] == RSS[2, i+1]
                opponent_strategy = IL[i+1,:]::Vector{Float64}
                for a in axes(capacity_equivalents, 1)
                    α = capacity_equivalents[a]
                    EV_matrix[IGS[α]]= EV_matrix[IGS[α]]+V[α]*opponent_strategy[OGS[α]]

                end
            else
            #Its counterpart is at i-1, by construction of RSS
                opponent_strategy = IL[i-1,:]::Vector{Float64}
                for a in axes(capacity_equivalents, 1)
                    α = capacity_equivalents[a]
                    EV_matrix[IGS[α]]= EV_matrix[IGS[α]]+V[α]*opponent_strategy[OGS[α]]

                end 
            end
        end
        # The expected value of each choice is generated
        # For each choice, we now subtract cost/benefit obtained from entering or leaving the market.
        for a in axes(EV_matrix,1)
            EV_matrix[a] = EV_matrix[a]-GEn[IGS[i], a]*ϕ+GEx[IGS[i], a]*ψ
        end 

        #PROBABILITY
        # Each choice is subject to a random shock, which reflects uncertainty about the future
        
        policy = ones(Float64, size(EV_matrix))
        for a in axes(policy,1), j in axes(policy,1)
            if a == j 
                nothing
            else

                policy[a] = policy[a]*CDF_convolution_IW(EV_matrix[a]-EV_matrix[j])

            end
        end
        policy_matrix[i, :] = policy
        
        value_vector[i]= sum(EV_matrix.*policy)
    end  
    return value_vector, policy_matrix
end

function iterator(iterations, c::Float64, exponent::Float64, β::Float64, II::Matrix{Float64}, CD::Matrix{Int}, GE::Matrix{Int}, RSS::Matrix{Int}, V::Vector{Float64}, CS, MI, CE, IL, n_m, OGS, IGS, GEn, GEx, ϕ, ψ, CCNCC, MS, minsyn, maxsyn)
    
    diff = zeros(Float64, size(RSS, 2))'
    av_investment_change = zeros(Float64, iterations)
    
    i_policy_diff = zeros(Float64, size(II, 1), size(II,2))
    
    av_ee_change = zeros(Float64, iterations)
    
    ee_policy_diff = zeros(Float64, size(IL, 1), size(IL,2))
    
    m_policy_change = zeros(Float64, iterations)
    
    a = 0
    M_P = zeros(Float64, size(RSS, 2))
    @showprogress for _=1:iterations
        
        a = a+1
        #MERGER STAGE
        pre_merger = II
        V, new_II = I_policy(c, exponent, β, II, CD, GE, RSS, V, CS, MI)[[1, 2]]
        #To test for convergence
        II = (II*(a-1)+new_II)/a
        #Average and standard deviation     
        i_policy_diff = pre_merger-new_II
        av_investment_change[a] = sum(broadcast(abs, i_policy_diff))/(size(i_policy_diff,1)*size(i_policy_diff,2))
        #ENTRY EXITY STAGE
        pre_IL = IL
        V, new_IL = EE_policy(RSS, CE, V, IL, n_m, OGS, IGS, GEn, GEx, ϕ, ψ)[[1,2]]
        IL = (IL*(a-1)+new_IL)/a
        ee_policy_diff = pre_IL - new_IL
        av_ee_change[a] = sum(broadcast(abs, ee_policy_diff))/(size(ee_policy_diff,1)*size(ee_policy_diff,2))
        #PRODUCTION STAGE
        V = V + CCNCC 
        diff = vcat(diff, V')
        #MERGER STAGE
        I_M_P = M_P
        V, M_P = M_policy(RSS, V, MS, minsyn, maxsyn)[[1,2]]
        change_M_P = broadcast(abs, M_P - I_M_P)
        m_policy_change[a] = sum(change_M_P)/(size(change_M_P,1)*size(change_M_P,2))
        
    end
    save("Value_[$n_m, $maxω, $iterations]_3_2.jld", "Values", V)
    save("Investment_policy_[$n_m, $maxω,$iterations]_3_2.jld", "Policies", II)
    save("Entry_Policy_[$n_m, $maxω,$iterations]_3_2.jld", "Policies", IL)
    save("Progress_[$n_m, $maxω,$iterations]_3_2.jld", "Values", diff[2:end, :])
    save("Merger_Policy_[$n_m, $maxω, $iterations]_3_2.jld", "Policies", M_P)
    save("Entry_Policy_conv_av_[$n_m, $maxω]_3_2.jld", "change", av_ee_change)
    save("Investment_Policy_conv_av_[$n_m, $maxω]_3_2.jld", "change", av_investment_change)
    save("Merger_Policy_conv_av_[$n_m, $maxω]_3_2.jld", "change", m_policy_change)
    return diff[end-7:end,:]
end

#function iterator(it, RSS, V, MS)
#    for _=1:it
#        A = M_policy(RSS, V, MS, -2, 2)
#        V = A
#    end
#    return V
#end

#iterator(1000, cost_investment, exponent, β, II, CD, GE, RSS, V, CS, MI, CE, IL, n_m, OGS, IGS, GEn, GEx, 15,13, CCNCC, MS, -1, 1)


#findall(x->x == 66, RSS[2,:])


CCNCC[130]
