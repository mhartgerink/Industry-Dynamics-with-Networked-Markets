using JLD, Distributions, StatsPlots, CSV, TypedTables, PlotlyJS
RSS = load("RSS_matrix_[3, 2].jld", "States")

n_m, maxω, iterations = 3,2, 1000

LC = load("LC_matrix_[$n_m, $maxω].jld", "Individual")
GE = load("Graph_Equiv_[$n_m, $maxω].jld", "Equiv")
CE = load("Capacity_Equiv_[$n_m, $maxω].jld", "Equiv")
IP = load("Investment_policy_[$n_m, $maxω,$iterations]_3_2.jld", "Policies")
#EP = load("Entry_Policy_[$n_m, $maxω,$iterations]_he.jld", "Policies")

#SIMULATION PLOTS 
SimMerg = load("merger_[$n_m]_le_lc.jld", "-")
SimIn = load("investment_[$n_m,$maxω]_le_lc.jld", "-")
SimLc = load("link_choice_[$n_m,$maxω]_le_lc.jld", "-")
SimCr = load("cap_real_[$n_m,$maxω]_le_lc.jld", "-")
SimSt = load("states_visited_[$n_m]_le_lc.jld", "-")

#function link_choice(SimLc)
#    firm_1 = histogram(collect(Int, 1:4), SimLc[1, :], ylabel = "Frequency", xlabel = "Investment", xlims = (-0.2,2.2), ylims = (0,6500) , legend = false, title = "Firm 1")
#    firm_2 = histogram(collect(Int, 1:4), SimLc[2, :], ylabel = "Frequency", xlabel = "Investment", xlims = (-0.2,2.2), ylims = (0,6500), legend = false, title = "Firm 2")
#    #final = plot(firm_1, firm_2)
#    #savefig(final, "Investment_2.svg")
#    return final
#end
#
#
##link_choice(SimSt)
#
#function frequency_tables(SimMerg, SimIn, SimLc, SimCr, SimSt)
#    f1_av_invest = zeros(11)
#    f2_av_invest = zeros(11)
#    for i in 1:20
#        f1_freq_investment = zeros(Int, 11)
#        f2_freq_investment = zeros(Int, 11)
#        for j in 1:11
#            f1_freq_investment[j] = size(findall(x->x==j-1, SimIn[(i-1)*2+1, :]), 1)
#            f2_freq_investment[j] = size(findall(x->x==j-1, SimIn[i*2, :]), 1)
#        end
#        f1_av_invest = f1_av_invest + f1_freq_investment./(size(SimIn,2)*20)
#        f2_av_invest = f2_av_invest + f2_freq_investment./(size(SimIn,2)*20)
#    end
#
#    f1_av_cap = zeros(11)
#    f2_av_cap = zeros(11)
#    for i in 1:20
#        f1_freq_cap = zeros(Int, 11)
#        f2_freq_cap = zeros(Int, 11)
#        for j in 1:11
#            f1_freq_cap[j] = size(findall(x->x==j-1, SimCr[(i-1)*2+1, :]), 1)
#            f2_freq_cap[j] = size(findall(x->x==j-1, SimCr[i*2, :]), 1)
#        end
#        f1_av_cap = f1_av_cap +f1_freq_cap./(size(SimCr,2)*20)
#        
#        f2_av_cap = f2_av_cap + f2_freq_cap./(size(SimCr,2)*20)
#    end
#    merger_perc = 0
#    for i in 1:20
#        merger_perc = merger_perc + 1-(size(findall(x->x==0, SimMerg[i, :]), 1)/size(SimMerg,2))
#    end
#
#    f1_link = zeros(8)
#    f2_link = zeros(8)
#    for i in 1:20
#        f1_freq_link = zeros(Int, 8)
#        f2_freq_link = zeros(Int, 8)
#        for j in 1:8
#            f1_freq_link[j] = size(findall(x->x==j, SimLc[(i-1)*2+1, :]), 1)
#            f2_freq_link[j] = size(findall(x->x==j, SimLc[i*2, :]), 1)
#        end
#        f1_link = f1_link + f1_freq_link./(size(SimLc,2)*20)
#        f2_link = f2_link + f2_freq_link./(size(SimLc,2)*20)
#    end
#    
#    investment = Table( F1 = f1_av_invest',F2 = f2_av_invest')
#    
#    cap = Table(F1 = f1_av_cap', F2 = f2_av_cap') 
#    link = Table(F1 = f1_link, f2 = f2_link)
#    merger = merger_perc/20
#    #CSV.write("In_3_2_he_hc.csv", investment; delim = ",")
#    #CSV.write("Cap_3_2_he_hc.csv", cap; delim = ",")
#    #CSV.write("Link_3_2_he_hc.csv", link; delim = ",")
#    println(merger)
#    return merger
#end
#
##frequency_tables(SimMerg, SimIn, SimLc, SimCr, SimSt)
function simplotter(matrix)
    plot()
    for i in 1:2
        v_1 = convert(Vector{Float64}, reverse(matrix[(i-1)*2+1, :]))
        v_2 = convert(Vector{Float64}, reverse(matrix[i*2, :]))
        
        for j in axes(v_1, 1)
            if j < 15
                nothing
            else
                v_1[j] = sum(v_1[j-14:j])/15
                v_2[j] = sum(v_2[j-14:j])/15
            end
        end

        plot!(v_1; ylims = (0, 2), label = "$i:Firm 1")
        display(plot!(v_2; ylims = (0, 2), label = "$i:Firm 2", title = "15-period Average Capacity"))
    end
    return 
end
simplotter(SimCr)
##
#function visited_state(SimSt, RSS)
#    count = zeros(Int, size(RSS, 2))
#    for i in axes(RSS, 2)
#        count[i] = size(findall(x->x==i, SimSt), 1)
#    end
#    so = sort(count, rev=true)
#    m_visited = []
#    times_visited = so[1:10]
#    for i in 1:10
#        m_visited= vcat(m_visited, findall(x->x== so[i], count)[1])
#    end
#    m_visited = hcat(times_visited, m_visited)
#    println(m_visited)
#    return m_visited
#end
#
#visited_state(SimSt, RSS)
#
#function Space_Reader(LC, RSS::Matrix{Int}, index::Int)
#    state = RSS[:, index]::Vector{Int}
#    actual_state = Matrix(vcat(LC[state[1], :]', LC[state[2], :]'))
#    return actual_state
#end
##
#function investment(LC, RSS, GE, i, IP)
#    firm_1 = zeros(3,3)
#    firm_2 = zeros(3,3)
#    investment = collect(0:2)
#    g = GE[Tuple(findfirst(x->x==i, GE))[1], :]
#    
#    for i in axes(g,1)
#        cap_1 = Space_Reader(LC, RSS, g[i])[1,1]
#        cap_2 = Space_Reader(LC, RSS, g[i])[2,1]
#
#        policy = IP[g[i], :]
#        
#        if RSS[1, g[i]] == RSS[2, g[i]]
#            opp_pol = policy
#        elseif RSS[1, g[i]] == RSS[2, g[i]+1]
#            opp_pol = IP[g[i]+1, :]
#        elseif RSS[1, g[i]] == RSS[2, g[i]-1]
#            opp_pol = IP[g[i]-1, :]
#        end
#        average = sum(policy.*investment)
#        o_av = sum(opp_pol.*investment)
#        firm_1[cap_1+1, cap_2+1] = average
#        firm_2[cap_1+1, cap_2+1] = o_av
#    end
#
#
#    f1 = PlotlyJS.surface(z = firm_1, opacity = 1, contours_z=attr(
#        show=true,
#        usecolormap=true,
#        
#        project_z=true))
#    layout = Layout(title="Investment Policy",
#        
#        scene=attr(
#        xaxis_title="Firm 1 Capacity",
#        yaxis_title="Firm 2 Capacity",
#        zaxis_title="Investment"
#        )
#    )
#    b = PlotlyJS.contour(z=firm_1', 
#    contours=attr(showlabels=true, xaxis_title="Firm 1 Capacity",
#    yaxis_title="Firm 2 Capacity",
#    zaxis_title="Investment"))
#    
#    x = PlotlyJS.plot(f1, layout)
#    
#    PlotlyJS.savefig(x, "511.svg")
#    open("511.html", "w") do io
#        PlotlyBase.to_html(io, x.plot)
#    end
#    return x
#end
#investment(LC, RSS, GE,511, IP)

#
#function link_cap(GE, index, EP, cap_2)
#    firm_1 = zeros(3, 4)
#    g = GE[Tuple(findfirst(x->x==index, GE))[1], :]
#    states= []
#    for i in g
#        if Space_Reader(LC, RSS, i)[2,1] == cap_2
#            states = vcat(states, i)
#        else
#            nothing
#        end
#    end
#    println(states)
#    for i in axes(states,1)
#        firm_1[i, :] = EP[states[i], :]
#    end
#    f1 = PlotlyJS.heatmap(z = firm_1, y=["0","1","2"],x=["1", "2", "3", "4"],  opacity = 1, contours_z=attr(
#        show=true,
#        usecolormap=true,
#        
#        project_z=true))
#    layout = Layout(title="Entry-, Exit Policy",
#        
#        
#        yaxis_title="Firm 1 Capacity",
#        xaxis_title="Link Choice",
#        zaxis_title="Probability"
#    
#    )
#    
#    
#    x = PlotlyJS.plot(f1, layout)
#    PlotlyJS.savefig(x, "s12_le_lc_2.svg")
#    open("s12_le_lc_2.html", "w") do io
#        PlotlyBase.to_html(io, x.plot)
#    end
#
#
#    return x
#end
#link_cap(GE, 1, EP, 2)


#Space_Reader(LC, RSS, 1)
#findall(x->x==0, SimMerg[5,:])









