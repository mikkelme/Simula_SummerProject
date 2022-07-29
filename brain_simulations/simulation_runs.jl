# Script for the various simulations I make
include("./evaluate.jl")
include("./data_analysis.jl")

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
savepath = "/Users/mikkelme/Documents/Github/Simula_SummerProject/figures/"
if !ispath(path)
    path = "/home/mikkelme/Simula_SummerProject/brain_simulations/"
    savepath = "/home/mikkelme/Simula_SummerProject/figures/"
end



# --- Default values (globally defined) --- #

# Brain Model [length unit: meter] 
lc = 1e-3 
arcLen = (100e-3, 0)
r_brain = 10e-3  
d_ratio = 1.5e-3/r_brain
r_curv = 50e-3 
inner_perturb = "(x,z) -> 0.0"  
outer_perturb = "(x,z) -> 0.0"  
BS_points = (1000, 0) 
field_Lc_lim = [1 / 2, 1]
field_Dist_lim = [1e-3, 5e-3] 

# PDE parameters 
μ = 0.8e-3  # Cerobrospinal fluid viscosity [Pa * s]
Κ = 1e-16   # Permeability in brain parenchyma [m^2] 
α = "(x) -> 1*μ/sqrt(Κ)" # Slip factor on Γ [Pa * s / m]
ps0 = "(x) -> x[1] < 0 ? 1*133.3224 : 0." # 1*mmHg [Pa]
∇pd0 = "(x) -> VectorValue(0.0, 0.0)" # Zero flux


# function run_solution_convergence(brain_param, PDE_param)
#     start_lc = 1e-2
#     end_lc = 5e-5
#     num_samples = 10
#     solution_convergence_vs_lc(brain_param, start_lc, end_lc, num_samples; filename = "solv_conv.txt")
# end


function run_flat()
    # Settings
    folder_name = "flat"
    start_width = 5e-3
    end_width = 0.5e-3
    num_samples = 10
    num_rad_lines = 100

    # Run
    brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
    eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; folder_name = folder_name)
    
    # Analyse     
    readpath = path * "data_" * folder_name * "/txt_files/"
    pressure_variance_vs_width(readpath *  "ps_radial_var.txt"; save = path * "data_" * folder_name * "/png_files/"  * "ps_var_width.png")
    pressure_variance_vs_angle(readpath * "ps_radial_var.txt"; save = path * "data_" * folder_name * "/png_files/" * "ps_var_angle.png")
    nflow_interface(readpath * "us_nflow.txt", save = path * "data_" * folder_name * "/png_files/" * "us_nflow_abs.png")
end

function run_single_sine()
    # Settings
    folder_name = "single_sine"
    start_width = 5e-3
    end_width = 0.5e-3
    num_samples = 10
    num_rad_lines = 100


    ω = 15.0/(arcLen[1] * (1 - d_ratio * r_brain / r_curv)) # 15 periods in total
    inner_perturb = @sprintf("(x,z) -> 0.3e-3 * cos(2*pi * abs(x) * %f)", ω)
    @show inner_perturb
    
    # Run
    brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
    eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; folder_name = folder_name)
    
    # Analyse     
    readpath = path * "data_" * folder_name * "/txt_files/"
    pressure_variance_vs_width(readpath *  "ps_radial_var.txt"; save = path * "data_" * folder_name * "/png_files/"  * "ps_var_width.png")
    pressure_variance_vs_angle(readpath * "ps_radial_var.txt"; save = path * "data_" * folder_name * "/png_files/"  * "ps_var_angle.png")
    nflow_interface(readpath * "us_nflow.txt", save = path * "data_" * folder_name * "/png_files/"  * "us_nflow_abs.png")
end




# function run_single_sine()
#     ...
# end

function run_inner_sines()
    # Settings
    start_width = 5e-3
    end_width = 0.5e-3
    num_samples = 10
    num_rad_lines = 100

    folder_name = "inner_sine"
    sines = ["(x,z) -> 0.3e-3 * cos(pi * abs(x) / 2e-3)", "(x,z) -> 0.3e-3 * cos(pi * abs(x) / 2e-3)"]

    # Run
    for (i, inner_perturb) in enumerate(sines)
        @show inner_perturb

        # brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
        # PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
        # eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; folder_name = folder_name)
    end
    # # Analyse     
    # readpath = path * "data_" * folder_name * "/txt_files/"
    # pressure_variance_vs_width(readpath *  "ps_radial_var.txt"; save = savepath * folder_name * "_" * "ps_var_width.png")
    # pressure_variance_vs_angle(readpath * "ps_radial_var.txt"; save = savepath * folder_name * "_" * "ps_var_angle.png")
    # nflow_interface(readpath * "us_nflow.txt", save = savepath * folder_name * "_" * "us_nflow_abs.png")
end



# run_flat()
run_single_sine()
# run_innter_sine()


# View model
# brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
# PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
# model, pgs_dict = create_brain(brain_param; view=true, write=false) 
