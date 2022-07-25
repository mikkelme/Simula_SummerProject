
using DelimitedFiles
using Plots

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
readpath = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/txt_files/"
savepath = "/Users/mikkelme/Documents/Github/Simula_SummerProject/figures/"
if !ispath(path)
    path = "/home/mikkelme/Simula_SummerProject/brain_simulations/"
    readpath = "/home/mikkelme/Simula_SummerProject/brain_simulations/txt_files/"
    savepath = "/home/mikkelme/Simula_SummerProject/figures/"
end





function plot_solution_convergence(filename)
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
    savefig(fig, savepath*"solution_convergence.png")
    display(fig)
      
    
end





function plot_pressure_variance(filename)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    # Get info
    infile = open(filename, "r")
    lines = readlines(infile)
    
    
    # number of samples
    words =  split(lines[23], '=')
    @assert words[1] == "# num_samples " "Reading wrong line for num_samples"
    num_samples = parse(Int, words[2][2:length(words[2])])

    # number of radial line
    words =  split(lines[24], '=')
    @assert words[1] == "# num_rad_lines " "Reading wrong line for num_rad_lines"
    num_rad_lines = parse(Int, words[2][2:length(words[2])])
 

    @show num_samples
    @show num_rad_lines
    width = []
    max_var = []
    mean_var = []

    for i in 1:num_samples
        start_index = (i-1)*num_rad_lines + 1
        end_index = start_index + num_rad_lines - 1

        max, idx = findmax(data_cells[start_index:end_index, 4])
        append!(width, data_cells[start_index + idx, 3])
        append!(max_var, max)
        append!(mean_var, sum(data_cells[start_index:end_index, 4])/num_rad_lines)

    end
    @show width
    @show max_var
    @show mean_var
    @show filename 


    # fig = plot(width)

    maxfig = plot(width, max_var, marker = :o, label="max variance")
    meanfig = plot(width, mean_var,  marker=:x, label="mean variance") # Get different axis for this one 

    # xlabel!(fig, "Resolution (lc)")
    # ylabel!(fig, "l²-norm with respect to finest resolution")
    # savefig(fig, savepath*"solution_convergence.png")
    display(maxfig)
    display(meanfig)

      
    
end



# plot_solution_convergence(readpath*"sol_conv.txt")


plot_pressure_variance(readpath * "25_7_15_31_ps_radial_var.txt")
