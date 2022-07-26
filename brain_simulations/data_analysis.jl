
using DelimitedFiles
using Plots

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
savepath = "/Users/mikkelme/Documents/Github/Simula_SummerProject/figures/"
if !ispath(path)
    path = "/home/mikkelme/Simula_SummerProject/brain_simulations/"
    savepath = "/home/mikkelme/Simula_SummerProject/figures/"
end



function solution_convergence(filename; save = false)
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
    ylabel!(fig, "l²-norm w.r.t. reference solution")
    save && savefig(fig, savepath*"solution_convergence.png")
    display(fig)
      
    
end



function pressure_variance_vs_width(filename; save = false)
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
        @show i, start_index, end_index
        max, idx = findmax(data_cells[start_index:end_index, 4])
        append!(width, data_cells[start_index + idx - 1, 3])
        append!(max_var, max)
        append!(mean_var, sum(data_cells[start_index:end_index, 4])/num_rad_lines)

    end
    
    fig = plot(width*1e3, max_var, marker = :o, label = "Max variance", ylabel = "Max variance",color = :red, 
    legend = :topleft, grid = :off,left_margin = 5Plots.mm, right_margin = 18Plots.mm)
    fig = plot!(twinx(), width*1e3, mean_var, marker = :x, label = "Mean variance", legend = :topright, 
    box = :on, grid = :off, ylabel = "Mean variance",left_margin = 5Plots.mm, right_margin = 18Plots.mm)
    plot!([1.5], seriestype = :vline, label = "Typical brain width (1.5 mm)", color = :black, linestyle = :dash)
    xlabel!("Width [mm]")
    save && savefig(savepath*"ps_var_width.png")
    display(fig)
      
    
end



function pressure_variance_vs_angle(filename; save = false)
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


   
    fig = plot()
    for i in 1:num_samples
        start_index = (i-1)*num_rad_lines + 1
        end_index = start_index + num_rad_lines - 1

    
        # Get angle and corresponding variance
        x_cord = data_cells[start_index:end_index, 1]
        y_cord = data_cells[start_index:end_index, 2]
        angle = atan.(y_cord ./x_cord)
        angle = ifelse.(angle .> 0, angle, angle .+ pi) # fix domain of tan output
        width =  sum(data_cells[start_index:end_index, 3])/num_rad_lines
        plot!(angle, data_cells[start_index:end_index, 4], label = @sprintf("width = %.2e", width))


        
    end
    

    
    xlabel!("Angle [rad]")
    ylabel!("Variance [Pa²]")
    save && savefig(savepath*"ps_var_angle.png")
    display(fig)
   


end


function nflow_interface(filename; save = false)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    m, n = size(data_cells)


    width = data_cells[1:m, 1]
    nflow = data_cells[1:m, 2]
    
    fig = plot(width*1e3, nflow*1e3, marker = :o)
    xlabel!("Width [mm]")
    ylabel!("u × n̂ [mm/s]")
    save && savefig(savepath*"nflow.png")
    display(fig)




end

data_folder =  "data_26_7_12_52"
# data_folder =  "data_ssh1e-4"

readpath = path * data_folder * "/txt_files/"


# solution_convergence(path * "/txt_files/"*"solution_converge.txt"; save = true)
# pressure_variance_vs_width(readpath * "ps_radial_var.txt"; save = false)
# pressure_variance_vs_angle(readpath * "ps_radial_var.txt"; save = false)
# nflow_interface(readpath * "us_nflow.txt")
