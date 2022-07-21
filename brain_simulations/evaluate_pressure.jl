

include("../mesh_generators/create_brain.jl")
include("./brain_sim_2D.jl")


path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/brain_simulations/"
end


function open_file(filename, brain_param, brain_PDE)
    # outfile = open(filename, "w")
    write(outfile, "# 2D brain simulation | pressure variance in radial direction under decrease of fluid path width\n")
    write(outfile, "# --- Brain parameters --- # ")
    @printf(outfile,"# lc = %f\n", brain_param.lc)
    @printf(outfile,"# arcLen = %f\n", brain_param.arcLen[1]) 
    @printf(outfile,"# r_brain = %f\n", brain_param.r_brain)
    write(outfile, "# d_ratio: varying parameter")
    @printf(outfile,"# r_curv = %f\n", brain_param.r_curv)
    @printf(outfile,"# inner_perturb = %f\n", brain_param.inner_perturb)
    @printf(outfile,"# outer_perturb = %f\n", brain_param.outer_perturb)
    @printf(outfile,"# BS_p = %f\n", brain_param.BS_points[1])
    @printf(outfile,"# field_Lc_lim = (%f, %f) \n", brain_param.field_Lc_lim[1], brain_param.field_Lc_lim[1])
    @printf(outfile,"# field_Dist_lim = (%f, %f) \n", brain_param.field_Dist_lim[1], brain_param.field_Dist_lim[1])




 



    outfile = open("pressure_reduction.txt", "w")




    write(outfile, "# info\n")
end

function add_sample_to_file()
    outfile = open("pressure_reduction.txt", "w")
    write(outfile, "# info\n")


    write(outfile, "mean_pos_x, mean_pos_y, rad_len, variance\n")
    for i in 1:num_rad_lines
        @printf(outfile,"%0.8f, %0.8f, %0.8f, %0.16f\n", mean_pos[i][1], mean_pos[i][2], rad_len[i], var[i]) 
       
    end

    
    close(outfile)

end

# per sample

# sample 1 
# mean_pos_x, mean_pos_y, rad_len,  

function evaluate_pressure(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; filename = "pressure_reduction.txt")
    @assert start_width <= brain_param.r_brain ["start width must be less than radial brain size"]
    @assert num_samples > 1  ["must have at least 2 sample points"]

    # outfile = open(filename, "w")
  


    width_diff = end_width - start_width 
    for i in LinRange(0, 1, num_samples)
        width = start_width + width_diff * i
        brain_param.d_ratio = width/brain_param.r_brain  
        # model, pgs_dict = create_brain(brain_param)
        # ush, psh, pdh = brain_PDE(model, pgs_dict, PDE_param)
        # mean_pos, rad_len, var = evaluate_radial_var(brain_param, psh, num_rad_lines)
        # @show mean_pos
        # @show rad_len
        # @show var

        return
        # return mean_pos, rad_len, var

    

    end






end












# --- Brain Model --- # 
lc = 0.1
arcLen = (5, 0)
r_brain = 2
d_ratio = 0.25
r_curv = 15 # cm
# inner_perturb(x, y) = 0.2 * cos(pi * abs(x) / 0.5) 
# outer_perturb(x, y) = 0.2 * cos(pi * abs(x) / 2)  
inner_perturb(x, y) = 0.0
outer_perturb(x, y) = 0.0
BS_points = (arcLen[1]*20, arcLen[2]*10)
field_Lc_lim = [1 / 2, 1]
field_Dist_lim = [0.1, 0.5]
brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)

# --- PDE parameters --- #
μ = 0.8e-3  # Fluid viscosity 
Κ = 1e-16   # Permeability in porous brain
α(x) = 1*μ/sqrt(Κ)

ps0(x) = x[1] < 0 ? 10 : 0 # go by g amplitude
fs0(x) = VectorValue(0.0, 0.0)
fd0(x) = 0.0 
∇pd0(x) = VectorValue(0.0, 0.0) # Zero flux
PDE_param = Dict(:μ => μ, :Κ => Κ, :α => α, :fs0 => fs0, :fd0 => fd0, :ps0 => ps0, :∇pd0 => ∇pd0) 


# --- Evaluation parameter --- #
start_width = 0.5
end_width = 0.2
num_samples = 3
num_rad_lines = 5

mean_pos, rad_len, var = evaluate_pressure(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines)

