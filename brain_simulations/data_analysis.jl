
using DelimitedFiles
using Plots

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/brain_simulations/"
end

savepath = "/Users/mikkelme/Documents/Github/Simula_SummerProject/figures/"




function plot_solution_convergence(filname)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    m, n = size(data_cells)

    lc = data_cells[2:m, 1]
    Δus = data_cells[2:m, 2]
    Δps = data_cells[2:m, 3]
    Δpd = data_cells[2:m, 4]

    fig = plot(lc, Δus, xaxis=:log,  yaxis=:log, marker=:o, label="us")
    fig = plot!(lc, Δps, xaxis=:log,  yaxis=:log, marker=:o, label="ps")
    fig = plot!(lc, Δpd, xaxis=:log,  yaxis=:log, marker=:o, label="pd")
    xlabel!(fig, "Resolution (lc)")
    ylabel!(fig, "l²-norm with respect to finest resolution")
    savefig(fig, savepath*"solution_convergence.pdf")
    display(fig)
    


    
    
    
end







filename = path * "txt_files/sol_conv.txt"
plot_solution_convergence(filename)