
using DelimitedFiles
using Plots
using Printf
using LaTeXStrings

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mikkelme/Simula_SummerProject/brain_simulations/"
end



function plot_solution_convergence(filename; save = false)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    m, n = size(data_cells)

    lc = data_cells[2:m, 1]
    Δus = data_cells[2:m, 2]
    Δps = data_cells[2:m, 3]
    Δpd = data_cells[2:m, 4]

    plot(lc*1e3, Δus, xaxis=:log, legend = :topleft, yaxis=:log, marker=:o, label="us")
    plot!(lc*1e3, Δps, xaxis=:log, legend = :topleft, yaxis=:log, marker=:o, label="ps")
    plot!(lc*1e3, Δpd, xaxis=:log, legend = :topleft, yaxis=:log, marker=:o, label="pd")
    xlabel!("Resolution (lc [mm])")
    ylabel!("l²-norm w.r.t. reference solution")
    save!=false && savefig(save)
    
end



function plot_ps_var_vs_width(filename, fig; label = "", var = "max", save = false)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    
    # Get info
    infile = open(filename, "r")
    lines = readlines(infile)
    
    # number of samples
    words =  split(lines[23], '=')
    @assert words[1] == "# num_samples " "Reading wrong line for num_samples"
    num_samples = parse(Int, strip(words[2], ' '))

    # number of radial line
    words =  split(lines[24], '=')
    @assert words[1] == "# num_rad_lines " "Reading wrong line for num_rad_lines"
    num_rad_lines = parse(Int, strip(words[2], ' '))

    # number of radial line
    words =  split(lines[25], '=')
    @assert words[1] == "# width " "Reading wrong line for width"
    words[2] = strip(words[2], [' ','[',']'])
    width = [parse(Float64, strip(word, ' ')) for word in split(words[2], ',')]
    
    # width = []
    max_var = []
    mean_var = []

    for i in 1:num_samples
        start_index = (i-1)*num_rad_lines + 1
        end_index = start_index + num_rad_lines - 1
        max, idx = findmax(data_cells[start_index:end_index, 4])
        # append!(width, data_cells[start_index + idx - 1, 3])
        append!(max_var, max)
        append!(mean_var, sum(data_cells[start_index:end_index, 4])/num_rad_lines)
    end

    if var == "max"
        y_val = max_var
        ylabel = "max variance [Pa²]"
    elseif var == "mean"
        y_val = mean_var 
        ylabel = "mean variance [Pa²]"
    else 
        @printf("keyword var = \"%s\" not understood", var)
        return
    end
    
    plot!(fig, width*1e3, y_val, marker = :o, legend = :topleft, label = label)
    
    if save!=false
        plot!(fig, [1.5], seriestype = :vline, label = "Typical width (1.5 mm)", color = :black, linestyle = :dash)
        xlabel!(fig, "Width [mm]")
        ylabel!(fig, ylabel)
        savefig(fig, save)
    end
          
end



function plot_ps_var_vs_angle(filename; save = false)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    
    # Get info
    infile = open(filename, "r")
    lines = readlines(infile)
    
    # number of samples
    words =  split(lines[23], '=')
    @assert words[1] == "# num_samples " "Reading wrong line for num_samples"
    num_samples = parse(Int, strip(words[2], ' '))

    # number of radial line
    words =  split(lines[24], '=')
    @assert words[1] == "# num_rad_lines " "Reading wrong line for num_rad_lines"
    num_rad_lines = parse(Int, strip(words[2], ' '))

   
    plot()
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


function plot_nflow_interface(filename, fig; label = "", save = false)
    data_cells, header_cells =  readdlm(filename, ',', Float64, '\n'; header = true, comments = true)
    m, n = size(data_cells)


    width = data_cells[1:m, 1]
    nflow = data_cells[1:m, 2]
    nflow_sqr = data_cells[1:m, 3]
    nflow_abs = sqrt.(nflow_sqr)

    plot!(fig, width*1e3, nflow_abs*1e3, marker = :o, legend = :topleft, label = label)

    if save!=false
        plot!(fig, [1.5], seriestype = :vline, label = "Typical width (1.5 mm)", color = :black, linestyle = :dash)
        xlabel!(fig, "Width [mm]")
        ylabel!(fig, "Normal flow: " * L"\sqrt{\int (u_S \cdot \hat{n})^2 d\Gamma}" * " [mm/s]")
        savefig(fig, save)
    end

end

function standard_analyse(readpath, folder_name)
    plot_ps_var_vs_width(readpath *  "ps_radial_var.txt", plot(); var = "max", save = path * "data_" * folder_name * "/png_files/"  * "ps_maxvar_width.png")
    plot_ps_var_vs_width(readpath *  "ps_radial_var.txt", plot(); var = "mean", save = path * "data_" * folder_name * "/png_files/"  * "ps_meanvar_width.png")
    plot_ps_var_vs_angle(readpath *  "ps_radial_var.txt"; save = path * "data_" * folder_name * "/png_files/"  * "ps_var_angle.png")
    plot_nflow_interface(readpath * "us_nflow.txt", plot(), save = path * "data_" * folder_name * "/png_files/" * "us_nflow_abs.png")
end

function combinned_analyse(savename, folder_names, labels)
     fig1 = plot()  
     fig2 = plot()
     fig3 = plot()
    for (i, folder_name) in enumerate(folder_names)
        readpath = path * "data_" * folder_name * "/txt_files/"
        savepath = path *  "/png_files/" * savename
    
        save1 = i == length(folder_names) ? savepath *  "_ps_maxvar_width.png" : false
        save2 = i == length(folder_names) ? savepath * "_ps_meanvar_width.png" : false
        save3 = i == length(folder_names) ? savepath * "_us_nflow_abs.png" : false
    
        plot_ps_var_vs_width(readpath * "ps_radial_var.txt", fig1; label = labels[i], var = "max", save = save1)
        plot_ps_var_vs_width(readpath * "ps_radial_var.txt", fig2; label = labels[i], var = "mean", save = save2)
        plot_nflow_interface(readpath * "us_nflow.txt", fig3, label = labels[i], save = save3)

    end
 
 



end



# data_folder =  "data_flat_curve/"
# # data_folder =  "data_ssh1e-4"

# data_folder = ["data_flat", "data_single_sine"]
# readpath = path * data_folder[2] * "/txt_files/"



# palette(:tab10)



# pressure_variance_vs_angle(readpath * "ps_radial_var.txt"; save = false)
