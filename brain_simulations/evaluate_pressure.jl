

include("../mesh_generators/create_brain.jl")
include("./brain_sim_2D.jl")


path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/brain_simulations/"
end


function open_file(filename, brain_param, PDE_param, width, num_samples, num_rad_lines)
    # outfile = open(filename, "w")

    outfile = open(path*"pressure_reduction.txt", "w")
    write(outfile, "# 2D brain simulation: Pressure variance in radial direction under decrease of fluid path width\n")
    write(outfile, "#\n# --- Brain parameters --- #\n")
    @printf(outfile,"# lc = %e\n", brain_param.lc)
    @printf(outfile,"# arcLen = %f\n", brain_param.arcLen[1]) 
    @printf(outfile,"# r_brain = %f\n", brain_param.r_brain)
    @printf(outfile, "# d_ratio = [")
    [@printf(outfile,"%f, ", width[i]/brain_param.r_curv) for i in 1:length(width)-1]
    @printf(outfile, "%f]\n", last(width)/brain_param.r_curv)

    @printf(outfile,"# r_curv = %f\n", brain_param.r_curv)
    # @printf(outfile,"# inner_perturb = %f\n", brain_param.inner_perturb)
    # @printf(outfile,"# outer_perturb = %f\n", brain_param.outer_perturb)
    @printf(outfile,"# BS_points = %d\n", brain_param.BS_points[1])
    @printf(outfile,"# field_Lc_lim = (%f, %f) \n", brain_param.field_Lc_lim[1], brain_param.field_Lc_lim[1])
    @printf(outfile,"# field_Dist_lim = (%f, %f) \n", brain_param.field_Dist_lim[1], brain_param.field_Dist_lim[1])

    write(outfile, "#\n# --- PDE parameters --- #\n")
    @printf(outfile,"# μ = %e\n", PDE_param[:μ])
    @printf(outfile,"# Κ = %e\n", PDE_param[:Κ])
    # @printf(outfile,"# α(x) = %f\n", PDE_param[:α])
    # @printf(outfile,"# ps0(X) = %f\n", PDE_param[:ps0])
    # @printf(outfile,"# ∇pd0(x) = %f\n", PDE_param[:∇pd0])
    
    # Print body of function
 
    write(outfile, "#\n# --- Sampling --- #\n")
    @printf(outfile,"# num_samples = %d\n", num_samples)
    @printf(outfile,"# num_rad_lines = %d\n", num_rad_lines)
    @printf(outfile, "# width = [")
    [@printf(outfile,"%.3e, ", width[i]) for i in 1:length(width)-1]
    @printf(outfile, "%.3e]\n", last(width))


    write(outfile, "#\n# --- Data --- #\n")
    write(outfile, "# mean_pos_x, mean_pos_y, rad_len, variance\n")


    return outfile
   
end

function add_sample_to_file(outfile, data, num_rad_lines)
    mean_pos, rad_len, var = data # unpack
    for i in 1:num_rad_lines
        @printf(outfile,"%e, %e, %e, %e\n", mean_pos[i][1], mean_pos[i][2], rad_len[i], var[i]) 
    end
end

# per sample

# sample 1 
# mean_pos_x, mean_pos_y, rad_len,  


function evaluate_radial_var(brain_param, psh, num_lines; degree = 2, view = false)
    rad_model, _ =  create_radial_lines(brain_param, num_lines; view = view)
    ip = Interpolable(psh)
  
    mean_pos = []
    rad_len = []
    var =[]
    for i in 1:num_lines
      # --- Get line triangulation and measure --- #
      L = Triangulation(rad_model, tags=[i])
      dL = Measure(L, 2)
  
      # --- Calculate metrics --- #
      len = sum(∫(1)*dL)                              # Line length
      mean_p = sum(∫(ip)*dL)/len                      # Mean pressure
      append!(mean_pos, [sum(∫( identity )*dL)/len])  # Mean position
      rel_diff = (psh - mean_p)/mean_p                # Relative difference to mean
      sq_diff = Interpolable((rel_diff)*(rel_diff))   # squared relative difference
      append!(var, sum(∫(sq_diff)*dL) / len)          # Variance
      append!(rad_len, len)
    end
    
    return mean_pos, rad_len, var
    
  end
  

function eval_pressure_var(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; filename = "pressure_reduction.txt")
    @assert start_width <= brain_param.r_brain ["start width must be less than radial brain size"]
    @assert num_samples > 1  ["must have at least 2 sample points"]

    
    width_diff = end_width - start_width 
    width = [start_width + width_diff * i for i in LinRange(0, 1, num_samples)]
    outfile = open_file(filename, brain_param, PDE_param, width, num_samples, num_rad_lines)
    for i in 1:num_samples
        @printf("i = %d/%d | width = %e\n", i, num_samples, width[i])
        brain_param.d_ratio = width[i]/brain_param.r_brain  
        model, pgs_dict = create_brain(brain_param)
        ush, psh, pdh = brain_PDE(model, pgs_dict, PDE_param)
        mean_pos, rad_len, var = evaluate_radial_var(brain_param, psh, num_rad_lines)
        add_sample_to_file(outfile, [mean_pos, rad_len, var], num_rad_lines)
   
        # # @show mean_pos
        # # @show rad_len
        # # @show var
        # return mean_pos, rad_len, var
        
        
        
    end
    close(outfile)


end





# --- Brain Model (length unit in meter) --- # 
lc = 1e-3 
arcLen = (100e-3, 0)
r_brain = 10e-3  
d_ratio = 1.5e-3/r_brain
r_curv = 150e-3 
inner_perturb(x, y) = 0.0 
outer_perturb(x, y) = 0.0 
BS_points = (20, 0) 
field_Lc_lim = [1 / 2, 1]
field_Dist_lim = [1e-3, 5e-3] 
brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)

# --- PDE parameters --- #
μ = 0.8e-3  # Cerobrospinal fluid viscosity [Pa * s]
Κ = 1e-16   # Permeability in brain parenchyma [m^2] 
α(x) = 1*μ/sqrt(Κ) # Slip factor on Γ [Pa * s / m]
ps0(x) = x[1] < 0 ? 10 : 0 # go by g amplitude
∇pd0(x) = VectorValue(0.0, 0.0) # Zero flux
PDE_param = Dict(:μ => μ, :Κ => Κ, :α => α, :ps0 => ps0, :∇pd0 => ∇pd0) 


# --- Evaluation parameter --- #
start_width = 5e-3
end_width = 0.5e-3
num_samples = 3
num_rad_lines = 10


# model, pgs_dict = create_brain(brain_param; view=true, write=false)
# ush, psh, pdh = brain_PDE(model, pgs_dict, PDE_param; write = true)

eval_pressure_var(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines)