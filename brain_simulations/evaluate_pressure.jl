

include("../mesh_generators/create_brain.jl")
include("./brain_sim_2D.jl")


function write_to_file()


end

function evaluate_pressure(brain_param, PDE_param, start_width, end_width, num_samples; filename = "pressure_reduction.txt")
    @assert start_width <= brain_param.r_brain ["start width must be less than radial brain size"]
    @assert num_samples > 1  ["must have at least 2 sample points"]

    open(filename, "w")

    width_diff = end_width - start_width 
    for i in LinRange(0, 1, num_samples)
        width = start_width + width_diff * i
        brain_param.d_ratio = width/brain_param.r_brain  
        model, pgs_dict = create_brain(brain_param)
        ush, psh, pdh = brain_PDE(model, pgs_dict, PDE_param)
        mean_pos, rad_len, var = evaluate_radial_var(num_lines)
        @show mean_pos
        @show rad_len
        @show var


        return mean_pos, rad_len, var

    

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


mean_pos, rad_len, var = evaluate_pressure(brain_param, PDE_param, 0.5, 0.2, 3)

# # --- Run simulation --- #
# model, pgs_dict = create_brain(brain_params; view=true, write=false)
# ush, psh, pdh = brain_PDE(model, pgs_dict, PDE_params; write = true)


# num_lines = 5
# evaluate_radial_var(num_lines)