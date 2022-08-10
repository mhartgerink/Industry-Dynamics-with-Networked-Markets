#= STATE SPACE MATRICES
***
This code generates the matrices for a specification of the model, for 2 firms, that are called in the 
code.
***
=#
using JLD, ProgressMeter, Random, SparseArrays
Random.seed!(123)

n_m, maxω = 3,2

function importmatrix(n_m, maxω)
    LC = load("LC_matrix_[$n_m, $maxω].jld", "Individual")
    AG = load("State_graphs_[$n_m, $maxω].jld", "stat_graph")
    RSS = load("RSS_matrix_[$n_m, $maxω].jld", "States")
    GE = load("Graph_Equiv_[$n_m, $maxω].jld", "Equiv")
    CE = load("Capacity_Equiv_[$n_m, $maxω].jld", "Equiv")
    CD = load("Capgraph_difference_[$n_m, $maxω].jld", "Diff")
    MI = load("Max_investment_[$maxω]_last.jld", "i_level")
    II = load("II_matrix_[$maxω]_last.jld", "Initial")
    CS = load("capacity_states_[$n_m, $maxω].jld", "capacity")
    IL = load("Link_strategy_[$n_m]_last.jld", "pr")
    IGS = load("I_Graph_state_[$n_m, $maxω].jld", "graph")
    OGS = load("O_Graph_state_[$n_m, $maxω].jld", "graph")
    return LC, AG, RSS, GE, CE, CD, MI, II, CS, IL, IGS, OGS
end


#LC, AG, RSS, GE, CE, CD, MI, II, CS, IL, IGS, OGS = importmatrix(n_m,maxω)[[1,2,3,4,5,6,7,8,9, 10, 11,12]]

#= Matrix for possible combinations of individual link and capacity choices.
    Each possible link vector is combined with each possible capacity choice and is ordered as:
    1: Lowest row-index for highest capacity choice
    2: Within each capacity choice, edges are lexicographically ordered.
    *Example*
    [2 1 1 1;
    2 1 1 0;
    2 1 0 1;
    2 1 0 0;
    2 0 0 0;
    etc. ]
    INPUT: number of markets, maximal capacity level
    OUTPUT: SparseMatrixCSC{Int64, Int64}
    =#

function LC_choices(n_m::Int, maxω::Int)
    #First generate all possible link vectors
    e = [1, 0]::Vector{Int}

    for _=1:(n_m-1)::Int
        e = vcat(sparse(hcat(ones(Int, size(e,1)), e)), sparse(hcat(zeros(Int, size(e, 1)), e)))
        # new matrix for 'i+1' markets
    end
    # Add the capacity choice
    ω_0 = zeros(Int, size(e,1))
    LC_matrix = hcat(ω_0, e)
    for i in 1:maxω::Int
        LC_matrix = convert(Matrix, vcat(hcat(fill(i, size(e,1)) ,e), LC_matrix))
    end
    save("LC_matrix_[$n_m, $maxω].jld", "Individual", LC_matrix)
    save("State_graphs_[$n_m, $maxω].jld", "stat_graph", e)
    return LC_matrix
end

LC = LC_choices(n_m, maxω)

#= Reduced state space for 2 firms
    ***
    We define the reduced state space by a set of column vectors of size n_f, where each element 
    corresponds to the index of the choice of the firm in the LC_matrix. The reduced state space therefore
    consists of all possible lexicographical orderings of size n_f with the domain (1+max\omega)*n_m
    =#

function Reduced_State_Space(LC, n_m::Int, maxω::Int)
    RSS_matrix = [0,0]
    for i in 1:size(LC, 1)::Int
        state = [i,i]
        for j in i+1:size(LC,1)::Int
            state = hcat(state, [i,j], [j,i])
        end
        RSS_matrix = hcat(RSS_matrix, state)
    end
    RSS_matrix = RSS_matrix[:, 2:end]
    save("RSS_matrix_[$n_m, $maxω].jld", "States", RSS_matrix)
    return RSS_matrix
end

RSS = Reduced_State_Space(LC, n_m, maxω)

#= GRAPH EQUIVALENT STATES
***
This matrix contains all graph equivalent states.
=#

function Graph_Equivalent_States(RSS, n_m::Int, maxω::Int)
    
    #All possible unique graphs#
    graphs = [0 0]
    for i in 1:2^n_m::Int, j in 1::Int:2^n_m::Int
        graphs = vcat(graphs, [i j])
    end
    graphs = graphs[1:end.!=1, :]
    #Equivalent link vectors
    equivalent_links = zeros(Int, 2^n_m, 1+maxω)
    for i in 1:2^n_m
        equivalent_links[i, 1:end]=collect(i:2^n_m:((1+maxω)*2^n_m))'
    end
    #Define the Reduced state space in terms of link vectors
    for i in axes(RSS, 2)
        RSS[:, i] = [Tuple(findfirst(x->x==RSS[1,i], equivalent_links))[1], Tuple(findfirst(x->x==RSS[2,i], equivalent_links))[1]]
    end
    
    graph_equivalents = zeros(Int, size(graphs,1), convert(Int, ceil(size(RSS, 2)/size(graphs,1))+(1+maxω)*2^n_m))
    for i in axes(graphs, 1)
        eq = []
        for j in axes(RSS, 2)
            if graphs[i,:] == RSS[:, j]
               eq = vcat(eq, j)
               
            else
                nothing
            end 
        end
        graph_equivalents[i, 1:size(eq', 2)]= eq'
    end

    graph_equivalents = sort!(graph_equivalents, dims = 2)
    i=1
    while sum(graph_equivalents[:, i]) ==0
        graph_equivalents = graph_equivalents[:, i+1:end]
    end
    save("Graph_Equiv_[$n_m, $maxω].jld", "Equiv", graph_equivalents)
    return graph_equivalents
end

GE = Graph_Equivalent_States(RSS, n_m, maxω)

#CAPACITY EQUIVALENTS

function ω_equivalents(n_m::Int, maxω::Int)
    #Import RSS_matrix
    RSS = load("RSS_matrix_[$n_m, $maxω].jld", "States")
    #All possible unique graphs#
    ω = [0 0]
    for i in 1:(1+maxω)::Int, j in 1::Int:(1+maxω)::Int
        ω = vcat(ω, [i j])
    end
    ω = ω[1:end.!=1, :]
    
    
    #equivalent capacity 
    ω_equiv = zeros(Int, 1+maxω, 2^n_m)
    for i in 1:(1+maxω)
        ω_equiv[i, 1:end]=collect((i-1)*2^n_m+1:1:i*2^n_m)'
    end 

    for i in axes(RSS, 2)
        RSS[:, i] = [Tuple(findfirst(x->x==RSS[1,i], ω_equiv))[1], Tuple(findfirst(x->x==RSS[2,i], ω_equiv))[1]]
    end

    ω_equivalents = zeros(Int, size(ω,1), convert(Int, size(RSS, 2)/size(ω,1)))
    for i in axes(ω, 1)
        eq = []
        for j in axes(RSS, 2)
            if ω[i,:] == RSS[:, j]
               eq = vcat(eq, j)  
            else
                nothing
            end 
        end
        ω_equivalents[i, 1:size(eq', 2)]= eq'
    end

    save("Capacity_Equiv_[$n_m, $maxω].jld", "Equiv", ω_equivalents)
    return ω_equivalents
end

CE = ω_equivalents(n_m, maxω)
#STATE SPACE READER 

function Space_Reader(LC, RSS::Matrix{Int}, index::Int)
    state = RSS[:, index]::Vector{Int}
    actual_state = Matrix(vcat(LC[state[1], :]', LC[state[2], :]'))
    return actual_state
end

#= INVESTMENT STAGE
    ***
We now construct the matrices for the investment stage=#

function max_investment(RSS::Matrix{Int}, LC::Matrix{Int}, maxω)
    max_investment = []
    for i in axes(RSS,2)
        max_investment = vcat(max_investment, maxω - Space_Reader(LC, RSS, i)[1][1])
        
    end
    save("Max_investment_[$maxω]_last.jld", "i_level", max_investment)
    return max_investment
end

MI = max_investment(RSS, LC, maxω)

function ω_difference(LC::Matrix{Int}, RSS::Matrix{Int}, GE::Matrix{Int})
    ω_difference = zeros(Int, size(RSS,2)*2, size(GE,2))::Matrix

    for i in axes(RSS, 2)
        state_ω = Space_Reader(LC, RSS, i)[:, 1]
        ge = GE[Tuple(findall(x->x==i, GE)[1])[1],:]
        for j in axes(ge,1)
            ι = ge[j]
            ω_difference[(i-1)*2+1:i*2,j] = (Space_Reader(LC, RSS, ι)[:, 1]-state_ω)::Vector{Int}
        end
    end
    maxω = Space_Reader(LC, RSS, 1)[1,1]
    n_m = size(Space_Reader(LC, RSS, 1),2)-1
    save("Capgraph_difference_[$n_m, $maxω].jld", "Diff", ω_difference)
    return ω_difference
end

CD = ω_difference(LC, RSS, GE)

function ω_states(LC::Matrix{Int}, RSS::Matrix{Int})

    ω = []
    for i in axes(RSS, 2)
        ω= vcat(ω, Space_Reader(LC, RSS, i)[1,1])
    end
    save("capacity_states_[3, 2].jld", "capacity", ω)
    return ω
end

CS = ω_states(LC, RSS)

function II_matrix(RSS::Matrix{Int}, MI::Vector, maxω)
    II_matrix = zeros(Float64, size(RSS,2), 1+maxω)
    @simd for i in 1:size(RSS, 2)
        II_matrix[i,1:(MI[i]+1)] = rand(Dirichlet((MI[i]+1),5),1)' 
    end
    save("II_matrix_[$maxω]_last.jld","Initial", II_matrix)
    return II_matrix
end

II = II_matrix(RSS, MI, maxω)

#ENTRY/EXIT STAGE

function IL_matrix(RSS::Matrix{Int}, n_m)
    IL_matrix = zeros(Float64, size(RSS, 2), 2^n_m)
    for i in axes(RSS, 2)
        IL_matrix[i, :] = rand(Dirichlet((2^n_m), 5), 1)
    end
    save("Link_strategy_[$n_m]_last.jld", "pr", IL_matrix)
    return IL_matrix
end

IL = IL_matrix(RSS, n_m)

function merger(LC, RSS, maxω)
    merger_states = zeros(Int, size(RSS, 2))
    b = findall(x->x==size(LC,1), RSS[2,:])
    
    for i in axes(RSS, 2)
        
        new_state = LC[RSS[1, i], :]+LC[RSS[2,i], :]
        if new_state[1] > maxω
            
            new_state[1] = maxω
        else
            nothing
        end

        for a in 2:size(new_state,1)
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
    save("Merger_states_[3, 2].jld", "states", merger_states)
    return merger_states

end
MS = merger(LC, n_m, maxω)

function RSS_graph(RSS::Matrix{Int}, n_m, maxω, LC::Matrix{Int})
    graphs = [0 0]
    for i in 1:2^n_m::Int, j in 1::Int:2^n_m::Int
        graphs = vcat(graphs, [i j])
    end
    graphs = graphs[1:end.!=1, :]
    #Equivalent link vectors
    equivalent_links = zeros(Int, 2^n_m, 1+maxω)
    graph_label = zeros(Int, 2^n_m, n_m)
    for i in 1:2^n_m
        equivalent_links[i, 1:end]=collect(i:2^n_m:((1+maxω)*2^n_m))'
        graph_label[i,:] = LC[i, 2:end]
    end
    #Define the Reduced state space in terms of link vectors
    for i in axes(RSS, 2)
        RSS[:, i] = [Tuple(findfirst(x->x==RSS[1,i], equivalent_links))[1], Tuple(findfirst(x->x==RSS[2,i], equivalent_links))[1]]
    end



    save("Graph_state_[$n_m, $maxω].jld", "graph", RSS)
    save("I_graph_state_[$n_m, $maxω].jld", "graph", RSS[1,:])
    save("O_graph_state_[$n_m, $maxω].jld", "graph", RSS[2,:])
    save("Link_label_[$n_m]_last.jld", "link", graph_label)
    return graph_label
end
Graph = RSS_graph(RSS, n_m, maxω, LC)
function EE_matrix(n_m)
    #consider all possible links
    e = [1, 0]::Vector{Int}

    for _=1:(n_m-1)::Int
        e = vcat(sparse(hcat(ones(Int, size(e,1)), e)), sparse(hcat(zeros(Int, size(e, 1)), e)))
        # new matrix for 'i+1' markets
    end
    entry_difference = zeros(Int, size(e,1), size(e,1))
    exit_difference = zeros(Int, size(e,1), size(e,1))
    for i in axes(e,1), j in axes(e, 1)
        entry_difference[i,j] = size(findall(x->x<0, e[i,:]-e[j,:]),1)
        exit_difference[i,j] = size(findall(x->x>0, e[i,:]-e[j,:]),1)
    end
    save("Entry_Diff_[$n_m]_last.jld", "diff", entry_difference)
    save("Exit_Diff_[$n_m]_last.jld", "diff", exit_difference)
    return entry_difference
end
EE = EE_matrix(n_m)
function equiv_check(CE, LC, RSS, i)
    b = zeros(size(CE[i,:], 1))
    for a in axes(CE[i, :], 1)
        b[a] = sum(Space_Reader(LC, RSS, i)[:, 1]-Space_Reader(LC, RSS, CE[i, a])[:,1])
    end
    return b
end








        
