using Plots, JLD, Distributions, StatsPlots
#Select the MPE
n_m, maxω, iterations = 3,2, 1000

#Import the matrices, computed in MPE_calculator for the correct name.
diff = load("Progress_[$n_m, $maxω,$iterations]_3_2.jld", "Values")
EP_conv = load("Entry_Policy_conv_av_[$n_m, $maxω]_3_2.jld", "change")
IP_conv = load("Investment_Policy_conv_av_[$n_m, $maxω]_3_2.jld", "change")
MP_conv = load("Merger_Policy_conv_av_[$n_m, $maxω]_3_2.jld", "change")

#The function below plots the convergence of the value functions
function convergence_plotter(diff, iterations)
    minus = diff[1:end-1, :]
    difference = diff[2:end, :].-minus
    av_diff = zeros(Float64, size(difference, 1))
    std_diff = zeros(Float64, size(difference, 1))
    for i in axes(av_diff, 1)
        av_diff[i] = sum(difference[i,:])/size(diff, 1)
        std_diff[i] = std(difference[i,:])
    end
    #difference convergence for full ranges
    a = plot(av_diff;
    #title = "Convergence between 1:$iterations",
    xlabel = "Iteration", 
    label = "Averaged Difference"
    )

    b = plot(std_diff;
    #title = "Standard Deviation between 1:$iterations",
    xlabel = "Iteration", 
    label = "St. Dev. of Difference"
    )
    full = plot(a, b, layout = (2,1), label = ["Averaged Difference" "St. Dev. of Difference"], plot_title = "Convergence between 1-$iterations")
    StatsPlots.savefig(full, "3_2_fullConvergence_he_hc.svg")
    
    #difference convergence for 2:10
    c = plot(collect(2:10),av_diff[2:10];
    #title = "Convergence between 2-10",
    xlabel = "Iteration",
    label = "Averaged Difference"
    )

    d = plot(collect(2:10),std_diff[2:10];
    #title = "Convergence between 2-10",
    xlabel = "Iteration",
    label = "St. Dev. of Difference"
    )
    beg = plot(c, d, layout = (2,1), label = ["Averaged Difference" "St. Dev. of Difference"], plot_title = "Convergence between 2-10")
    StatsPlots.savefig(beg, "3_2_beginConvergence_he_hc.svg")
    

    #difference for end ranges
    e = plot(collect(size(av_diff,1)-50:size(av_diff,1)), av_diff[end-50:end]; 
    #title = "Convergence between 150-200", xlabel = "Iteration", label = "Averaged Difference"
    )
    start = iterations - 50
    f = plot(collect(size(std_diff,1)-50:size(std_diff,1)), std_diff[end-50:end]; 
    #title = "Convergence between $start-$iterations", xlabel = "Iteration", label = "St. Dev. of Difference"
    )
    en = plot(e, f, layout = (2,1), label = ["Averaged Difference" "St. Dev. of Difference"], plot_title = "Convergence between $start-$iterations")
    StatsPlots.savefig(en, "3_2_endConvergence_he_hc.svg")
    

    return nothing
end
convergence_plotter(diff, iterations)

#This function plots the convergence of the policy functions.
function policy_convergence(vector, string, vector2, string2, vector3, string3)
    a = plot(vector;
    xlabel = "Iterations",
    label = "Averaged Difference",
    ylims = (0,0.1),
    title = "$string"
    )
    #savefig(a, "le_he_$string.svg")

    b = plot(vector2;
    xlabel = "Iterations",
    label = "Averaged Difference",
    ylims = (0,0.1),
    title = "$string2"
    )
    #savefig(a, "le_he_$string2.svg")

    c = plot(vector3;
    xlabel = "Iterations",
    label = "Averaged Difference",
    ylims = (0,0.1),
    title = "$string3"
    )
    #savefig(a, "le_le_$string3.svg")
    final = plot(a, b, c, layout = (3,1))
    StatsPlots.savefig(final, "3_2_Policies_he_he_conv.svg")

    return nothing
end

policy_convergence(IP_conv, "Investment Policy Convergence", EP_conv, "Entry/Exit Policy Convergence", MP_conv, "Merger Policy Convergence")


