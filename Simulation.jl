using Random, StatsBase, Distributions, ProgressMeter
Random.seed!(123)

n_m, maxω, iterations = 3,2, 1000

RSS = load("RSS_matrix_[$n_m, $maxω].jld", "States")
LC = load("LC_matrix_[$n_m, $maxω].jld", "Individual")
IP = load("Investment_policy_[$n_m, $maxω,$iterations]_3_2.jld", "Policies")
EP = load("Entry_Policy_[$n_m, $maxω,$iterations]_3_2.jld", "Policies")
MP = load("Merger_Policy_[$n_m, $maxω, $iterations]_3_2.jld", "Policies")
IGS = load("I_Graph_state_[$n_m, $maxω].jld", "graph")
OGS = load("O_Graph_state_[$n_m, $maxω].jld", "graph")


function states(LC)
    #Define the unique states for both firms
    uRSS_matrix = [0,0]
    for i in 1:size(LC, 1)::Int
        state = [i,i]
        for j in i+1:size(LC,1)::Int
            state = hcat(state, [i,j])
        end
        uRSS_matrix = hcat(uRSS_matrix, state)
    end
    uRSS_matrix = uRSS_matrix[:, 2:end]
    return uRSS_matrix
end

uRSS = states(LC)

function policy_norm(policy1)
    policy1 = policy1./sum(policy1)
    return policy1
end

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
    return merger_states

end

function state_identifier(cap_1::Int, cap_2::Int, link_1, link_2, LC, RSS)
    begin_1 = findfirst(x->x==cap_1, LC[:, 1])
    begin_2 = findfirst(x->x==cap_2, LC[:, 1])
    choice_1 = begin_1+link_1-1
    choice_2 = begin_2+link_2-1
    
    newstate = 0
    for i in axes(RSS,2)
        if RSS[1,i] == choice_1 && RSS[2,i] == choice_2
            newstate = i
        else 
            nothing
        end
    end
    return newstate
end

MS = merger(LC, RSS, maxω)

function simulation(time, uRSS, RSS, LC, IP, EP, MP, maxω, MS, state)
    merger = []
    investment = zeros(Int, 2)
    link_choice = zeros(Int, 2)
    cap_real = zeros(Int, 2)
    #state = rand(collect(1:size(uRSS,2)))
    states_visited = [state]
    invest = collect(0:1:size(IP,2)-1)
    for i in 1:time
        #MERGER STAGE
        premstate = states_visited[1]
        if RSS[1, premstate] == RSS[2, premstate]
            m_1 = policy_norm(MP[premstate])
            m_2 = policy_norm(MP[premstate])
        elseif RSS[1, premstate] == RSS[2, premstate+1]
            m_1 = policy_norm(MP[premstate])
            m_2 = policy_norm(MP[premstate+1])
        elseif RSS[1, premstate] == RSS[2, premstate-1]
            m_1 = policy_norm(MP[premstate])
            m_2 = policy_norm(MP[premstate-1])
        end
        prop = rand([1,2])
        if prop == 1
            mstate = sample([premstate, MS[premstate]],weights([1-m_1, m_1]))
        elseif prop == 2
            mstate = sample([premstate, MS[premstate]],weights([1-m_2, m_2]))
        end

        if mstate == premstate
            merger = vcat(0, merger)
        else
            merger = vcat(mstate, merger)
        end
        states_visited = vcat(mstate, states_visited)
        preeestate = states_visited[1]
        
        #ENTRY/EXIT STAGE
        if RSS[1, preeestate] == RSS[2, preeestate]
            E_1 = policy_norm(EP[preeestate, :])
            E_2 = policy_norm(EP[preeestate, :])
        elseif RSS[1, preeestate] == RSS[2, preeestate+1] 
            E_1 = policy_norm(EP[preeestate, :])
            E_2 = policy_norm(EP[preeestate+1, :])
        elseif RSS[1, preeestate] == RSS[2, preeestate-1]
            E_1 = policy_norm(EP[preeestate, :])
            E_2 = policy_norm(EP[preeestate-1, :])
        end
        
        link_1 = sample(collect(1:1:size(EP,2)), weights(E_1))
        link_2 = sample(collect(1:1:size(EP,2)), weights(E_2))
        newgraph = vcat(link_1, link_2)
        link_choice = hcat(newgraph, link_choice)
        cap_1, cap_2 = LC[RSS[1, preeestate], 1],LC[RSS[2, preeestate], 1]
        
        eestate = state_identifier(cap_1, cap_2, link_1, link_2, LC, RSS)
        states_visited= vcat(eestate, states_visited)
        
        #INVESTMENT STAGE
        preistate = states_visited[1]
        if RSS[1, preistate] == RSS[2, preistate]
            f1_policy = policy_norm(IP[preistate,:])
            cap_1 = LC[RSS[1,preistate], 1]
            f2_policy = policy_norm(IP[preistate,:])
            cap_2 = LC[RSS[1,preistate], 1]
        elseif RSS[1, preistate] == RSS[2, preistate+1]
            f1_policy = policy_norm(IP[preistate,:])
            cap_1 = LC[RSS[1,preistate], 1]
            f2_policy = policy_norm(IP[preistate+1,:])
            cap_2 = LC[RSS[1,preistate+1], 1]
        elseif RSS[1, preistate] == RSS[2, preistate-1]
            f1_policy = policy_norm(IP[preistate,:])
            cap_2 = LC[RSS[1,preistate], 1]
            f2_policy = policy_norm(IP[preistate-1,:])
            cap_2 = LC[RSS[1,preistate-1], 1]
        end
             

        x1 = sample(invest, weights(f1_policy))
        x2 = sample(invest, weights(f2_policy))

        nc_1 = rand(Binomial(x1, 0.45),1)[1]
        nc_2 = rand(Binomial(x2, 0.45),1)[1]

        state_invest = vcat(x1, x2)
        investment = hcat(state_invest, investment)

        d1 = rand(Binomial(cap_1, 0.3), 1)[1]
        d2 = rand(Binomial(cap_2, 0.3), 1)[1]
        
        cap_1 = cap_1-d1+nc_1
        cap_2 = cap_2-d2+nc_2
        
        if cap_1 > maxω
            cap_1 = maxω
        elseif cap_1 < 0
            cap_1 = 0
        end

        if cap_2 > maxω
            cap_2 = maxω
        elseif cap_2 < 0
            cap_2 = 0
        end
        cap_real = hcat(vcat(cap_1, cap_2), cap_real)
        istate = state_identifier(cap_1, cap_2, link_1, link_2, LC, RSS)
        states_visited= vcat(istate, states_visited)
    end
    
    return merger[1:end-1], investment[:, 1:end-1], link_choice[:, 1:end-1], cap_real[:, 1:end-1], states_visited[1:end]
end

function multiple_initial(iter, uRSS, RSS, LC, IP, EP, MP, maxω, MS, size, i_states)
    
    states_visited = zeros(Int, size, 3*iter+1)
    investment = zeros(Int, 2*size, iter)
    link_choice = zeros(Int, 2*size, iter)
    merger = zeros(Int, size, iter-1)
    cap_real = zeros(Int, 2*size, iter)

    @showprogress for i in 1:size
        A, B, C, D, E = simulation(iter, uRSS, RSS, LC, IP, EP, MP, maxω, MS, i_states[i])[[1,2,3,4,5]]
        states_visited[i,:] = E
        investment[(i-1)*2+1:i*2, :] = B
        link_choice[(i-1)*2+1:i*2, :] = C
        merger[i, :] = A 
        cap_real[(i-1)*2+1:i*2, :] = D
    end
    save("states_visited_[3,2]_le_lc.jld", "-", states_visited)
    save("investment_[3,2]_le_lc.jld", "-", investment)
    save("link_choice_[3,2]_le_lc.jld", "-", link_choice)
    save("merger_[3,2]_le_lc.jld", "-", merger)
    save("cap_real_[3,2]_le_lc.jld", "-", cap_real)

    return link_choice
end
i_states = sample(1:size(uRSS, 2), 20, replace = false)
multiple_initial(5000, uRSS, RSS, LC, IP, EP, MP, maxω, MS, 20, i_states)






