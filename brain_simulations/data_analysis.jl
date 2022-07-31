
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
savepath = "/Users/mikkelme/Documents/Github/Simula_SummerProject/figures/"
if !ispath(path)
    path = "/home/mikkelme/Simula_SummerProject/brain_simulations/"
    savepath = "/home/mikkelme/Simula_SummerProject/figures/"
end



function plot_solution_convergence(filename; save = false)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    m, n = size(data_cells)

    lc = data_cells[2:m, 1]
    Δus = data_cells[2:m, 2]
    Δps = data_cells[2:m, 3]
    Δpd = data_cells[2:m, 4]

    plot(lc, Δus, xaxis=:log,  yaxis=:log, marker=:o, label="us")
    plot!(lc, Δps, xaxis=:log,  yaxis=:log, marker=:o, label="ps")
    plot!(lc, Δpd, xaxis=:log,  yaxis=:log, marker=:o, label="pd")
    xlabel!(fig, "Resolution (lc)")
    ylabel!(fig, "l²-norm w.r.t. reference solution")
    save && savefig(save)
    
end



function pressure_variance_vs_width(filename; label = "", var = "max", makefig = false, save = false)
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
 
    width = []
    max_var = []
    mean_var = []

    for i in 1:num_samples
        start_index = (i-1)*num_rad_lines + 1
        end_index = start_index + num_rad_lines - 1
        max, idx = findmax(data_cells[start_index:end_index, 4])
        append!(width, data_cells[start_index + idx - 1, 3])
        append!(max_var, max)
        append!(mean_var, sum(data_cells[start_index:end_index, 4])/num_rad_lines)
    end

    if var == "max"
        y_val = max_var
        ylabel = "max variance [Pa²]"
    elseif var == "mean"
        y_val = mean_var 
        ylabel = "mean variance [Pa²]"
    elseif var == "multi"
        ylabel = nothing
    else 
        @printf("keyword var = %s not understood", var)
        return
    end
    
    makefig && plot()
    plot!(width*1e3, y_val, marker = :o, legend = :topleft, label = label)

    # lm = 5Plots.mm
    # rm = 25Plots.mm
    # plot!(fig, width*1e3, max_var, marker = :o, label = "Max variance", ylabel = "Max variance", 
    # legend = :topleft, grid = :off,left_margin = lm, right_margin = rm, palette = :Dark2_5)
    # plot!(tw, width*1e3, mean_var, marker = :x, label = "Mean variance", legend = :topright, 
    # box = :on, grid = :off, ylabel = "Mean variance",left_margin = lm, right_margin = rm, palette = :Dark2_5)
    
    if save!=false
        plot!([1.5], seriestype = :vline, label = "Typical width (1.5 mm)", color = :black, linestyle = :dash)
        xlabel!("Width [mm]")
        !isnothing(ylabel) && ylabel!(ylabel)
        savefig(save)
    end
          
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
        plot!(angle, data_cells[start_index:end_index, 4], label = @sprintf("width = %.2f mm", width*1e3))
 
    end 
    xlabel!("Angle [rad]")
    ylabel!("Variance [Pa²]")
    save!=false && savefig(save)
end


function nflow_interface(filename; label = "", makefig = false, save = false)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    m, n = size(data_cells)


    width = data_cells[1:m, 1]
    nflow = data_cells[1:m, 2]
    nflow_sqr = data_cells[1:m, 3]
    nflow_abs = sqrt.(nflow_sqr)

    makefig && plot()
    plot!(width*1e3, nflow_abs*1e3, marker = :o, legend = :topleft, label = label)

    if save!=false
        plot!([1.5], seriestype = :vline, label = "Typical width (1.5 mm)", color = :black, linestyle = :dash)
        xlabel!("Width [mm]")
        ylabel!("Normal flow: " * L"\sqrt{\int (u_S \cdot \hat{n})^2 d\Gamma}" * " [mm/s]")
        savefig(save)
    end

end



# data_folder =  "data_flat_curve/"
# # data_folder =  "data_ssh1e-4"

data_folder = ["data_flat", "data_single_sine"]
# readpath = path * data_folder[2] * "/txt_files/"



# palette(:tab10)

pressure_variance_vs_width(path * data_folder[1] * "/txt_files/" * "ps_radial_var.txt"; label = "flat", var = "mean", makefig = true, save = false)
pressure_variance_vs_width(path * data_folder[2] * "/txt_files/" * "ps_radial_var.txt"; label = "sine", var = "mean", save = savepath * "test.png")

nflow_interface(path * data_folder[1] * "/txt_files/" * "us_nflow.txt"; label = "flat", makefig = true, save = false)
nflow_interface(path * data_folder[2] * "/txt_files/" * "us_nflow.txt"; label = "sine", save = savepath * "test2.png")


# pressure_variance_vs_angle(readpath * "ps_radial_var.txt"; save = false)
