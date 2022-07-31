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




function run_solution_convergence()
    # Settings 
    start_lc = 1e-2
    end_lc = 5e-5
    num_samples = 10
    filename = "test_conv"

    # Run
    brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
    solution_convergence_vs_lc(brain_param, PDE_param, start_lc, end_lc, num_samples; filename = filename * ".txt")
    plot_solution_convergence(path * "/txt_files/" * filename * ".txt"; save = path * "/png_files/" * filename * ".png")
end


function run_flat()
    # Settings
    folder_name = "flat"
    start_width = 5e-3
    end_width = 0.5e-3
    num_samples = 4
    num_rad_lines = 100

    # Run
    brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
    eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; folder_name = folder_name)
    
    # Standard analyse     
    readpath = path * "data_" * folder_name * "/txt_files/"
    standard_analyse(readpath, folder_name)

end

function run_single_inner_negsine()
    # Settings
    folder_name = "single_inner_negsine"
    start_width = 5e-3
    end_width = 0.5e-3
    num_samples = 5
    num_rad_lines = 100

    A = 0.3e-3
    f = 15.0 # 15 periods in total
    ω = 2*pi * f/(arcLen[1] * (1 - d_ratio * r_brain / r_curv)) # 15 periods in total
    inner_perturb = @sprintf("(x,z) -> %f * sin(abs(x) * %f) * fld(mod2pi(abs(x) * %f),pi) ", A , ω, ω)
    
    # Run
    brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
    eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; folder_name = folder_name)
    
    # Analyse     
    readpath = path * "data_" * folder_name * "/txt_files/"
    standard_analyse(readpath, folder_name)
end


function run_single_inner_sine()
    # Settings
    folder_name = "single_inner_sine"
    start_width = 5e-3
    end_width = 0.5e-3
    num_samples = 5
    num_rad_lines = 100

    A = 0.3e-3
    f = 15.0 # 15 periods in total
    ω = 2*pi * f/(arcLen[1] * (1 - d_ratio * r_brain / r_curv)) 
    inner_perturb = @sprintf("(x,z) -> %f * sin(abs(x) * %f)", A , ω)

    # Run
    brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
    eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; folder_name = folder_name)
    
    # Analyse     
    readpath = path * "data_" * folder_name * "/txt_files/"
    standard_analyse(readpath, folder_name)
end


# Complete run_inner_sines for running multiple sine function as inner_perturb
# Then make plotter function for plotting multiple lc_decrease_simulation in one plot. 
# Then you can make combinned plots for different frequency different amplitude, or different wave form (thus 3 plots)

function run_inner_negsines()
    # Settings
    start_width = 5e-3
    end_width = 0.5e-3
    num_samples = 3
    num_rad_lines = 10

    folder_name = "inner_negsines"
    A = 0.3e-3
    freq = [5.0, 15.0, 25.0]
    ω(f) = 2*pi * f/(arcLen[1] * (1 - d_ratio * r_brain / r_curv)) 
    # inner_perturb = @sprintf("(x,z) -> %f * sin(abs(x) * %f)", A , ω)

    negsines = [@sprintf("(x,z) -> %f * sin(abs(x) * %f) * fld(mod2pi(abs(x) * %f),pi) ", A , ω(f), ω(f)) for f in freq]
    
    # # Run
    # for (i, inner_perturb) in enumerate(negsines)
    #     i_folder_name = folder_name * @sprintf("_f%s", freq[i])

    #     # Run
    #     brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    #     PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
    #     eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; folder_name = i_folder_name)
        
    #     # Analyse     
    #     readpath = path * "data_" * i_folder_name * "/txt_files/"
    #     standard_analyse(readpath, i_folder_name)
               
    # end

    folder_names = [folder_name * @sprintf("_f%s", f) for f in freq]
    labels = [@sprintf("f = %s", f) for f in freq]
    combinned_analyse(folder_names, labels)
    # # Combinned analyse   
    # fig1 = plot()  
    # fig2 = plot()
    # fig3 = plot()


    # for (i, inner_perturb) in enumerate(negsines)
    #     i_folder_name = folder_name * @sprintf("_f%s", freq[i])
    #     readpath = path * "data_" * i_folder_name * "/txt_files/"

       
    #     save1 = i == length(negsines) ? savepath * "test1.png" : false
    #     save2 = i == length(negsines) ? savepath * "test2.png" : false
    #     save3 = i == length(negsines) ? savepath * "test3.png" : false

  
    #     plot_ps_var_vs_width(readpath * "ps_radial_var.txt", fig1; label = @sprintf("f = %.2f", freq[i]), var = "mean", save = save1)
    #     plot_ps_var_vs_width(readpath * "ps_radial_var.txt", fig2; label = @sprintf("f = %.2f", freq[i]), var = "max", save = save2)
    #     plot_nflow_interface(readpath * "us_nflow.txt", fig3, label = @sprintf("f = %.2f", freq[i]), save = save3)

    # end



  
end


# --- View model --- #
# A = 0.3e-3
# f = 15.0
# ω = 2*pi * f/(arcLen[1] * (1 - d_ratio * r_brain / r_curv)) # 15 periods in total
# d_ratio = 0.5e-3 / r_brain 
# inner_perturb = @sprintf("(x,z) -> %f * sin(abs(x) * %f) * fld(mod2pi(abs(x) * %f),pi) ", A , ω, ω)
# brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
# model, pgs_dict = create_brain(brain_param; view=true, write=false) 


# --- Simulaiton runs --- #
# run_solution_convergence()
# run_flat()
# run_single_inner_sine()
# run_single_inner_negsine()
run_inner_negsines()


