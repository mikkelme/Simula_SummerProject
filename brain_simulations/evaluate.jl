using Dates
using Printf


include("../mesh_generators/create_brain.jl")
include("./brain_sim_2D.jl")
include("../mesh_generators/generate_radial_lines.jl")
include("../mesh_generators/generate_center_line.jl")
include("./interface_evaluation.jl")


path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mikkelme/Simula_SummerProject/brain_simulations/" # SSH
end



function create_file(path, filename, title, brain_param, PDE_param)    
    # Create file
    outfile = open(path*filename, "w") 
    
    # Write general info of system
    @printf(outfile,"# %s\n", title)
    
    # Geometry parameters
    write(outfile, "#\n# --- Brain parameters (fixed) --- #\n")
    @printf(outfile,"# lc = %e \n", brain_param.lc)
    @printf(outfile,"# arcLen = %f m\n", brain_param.arcLen[1]) 
    @printf(outfile,"# r_brain = %f m\n", brain_param.r_brain)
    @printf(outfile,"# d_ratio = %f\n", brain_param.d_ratio)
    @printf(outfile,"# r_curv = %f m\n", brain_param.r_curv)
    @printf(outfile,"# inner_perturb = %s\n", brain_param.inner_perturb_body)
    @printf(outfile,"# outer_perturb = %s\n", brain_param.outer_perturb_body)
    @printf(outfile,"# BS_points = %d\n", brain_param.BS_points[1])
    @printf(outfile,"# field_Lc_lim = (%f, %f) \n", brain_param.field_Lc_lim[1], brain_param.field_Lc_lim[2])
    @printf(outfile,"# field_Dist_lim = (%f, %f) m \n", brain_param.field_Dist_lim[1], brain_param.field_Dist_lim[2])
    
    # PDE parameters
    write(outfile, "#\n# --- PDE parameters (fixed=) --- #\n")
    @printf(outfile,"# μ = %e Pa*s\n", PDE_param.μ)
    @printf(outfile,"# Κ = %e m^2\n", PDE_param.Κ)
    @printf(outfile,"# α = %s\n", PDE_param.α_body)
    @printf(outfile,"# ps0 = %s\n", PDE_param.ps0_body)
    @printf(outfile,"# ∇pd0 = %s\n", PDE_param.∇pd0_body)
    
    return outfile
    
end


function add_sample_to_file(path, filename, data, vector_length = 0)
    outfile = open(path * filename, "a") # open in append mode
    data_len = length(data)
    if vector_length == 0 
        for i in 1:data_len-1
            @printf(outfile, "%e, ", data[i]) 
        end
        @printf(outfile, "%e\n", data[data_len]) 
    else 
        for j in 1:vector_length
            for i in 1:data_len-1
                @printf(outfile, "%e, ", data[i][j]) 
            end
            @printf(outfile, "%e\n", data[data_len][j] ) 
        end
    end
    close(outfile) 
end


function calc_solution_diff(brain_param, coarse_res, fine_res; degree = 2, num_nearest_vertices = 4)
    # Unpack inputs
    coarse_us, coarse_ps, coarse_pd = coarse_res
    fine_us, fine_ps, fine_pd = fine_res

    # Set searchmethod
    sm = searchmethod=KDTreeSearch(num_nearest_vertices=num_nearest_vertices)

    #--- Stokes domain ---#
    Stokes_CL, _ = create_centerline(brain_param, "S"; view = false)
    SL = Triangulation(Stokes_CL)
    dSL = Measure(SL, degree) 
    Slen = sum(∫(1)*dSL)   # Length of centerline
    # writevtk(SL, path*"SL")

    # us
    us_diff =  Interpolable(coarse_us, searchmethod = sm) - fine_us
    Δus_l2norm = sqrt(sum(∫(Interpolable(us_diff ⋅ us_diff, searchmethod = sm)) * dSL))/Slen 

    # #ps
    ps_diff =  Interpolable(coarse_ps, searchmethod = sm) - fine_ps
    Δps_l2norm = sqrt(sum(∫(Interpolable(ps_diff * ps_diff, searchmethod = sm)) * dSL))/Slen 

    #--- Darcy domain ---#
    Darcy_CL, _ = create_centerline(brain_param, "D"; view = false)
    DL = Triangulation(Darcy_CL)
    dDL = Measure(DL, degree) 
    Dlen = sum(∫(1)*dDL)   # Length of centerline
    # writevtk(DL, path*"DL")

    # pd
    pd_diff =  Interpolable(coarse_pd, searchmethod = sm) - fine_pd
    Δpd_l2norm = sqrt(sum(∫(Interpolable(pd_diff * pd_diff, searchmethod = sm)) * dDL))/Dlen 

    return [Δus_l2norm, Δps_l2norm, Δpd_l2norm]


    
end

function evaluate_radial_var(brain_param, ps, num_lines; degree = 2, view = false, num_nearest_vertices = 4)
    rad_model, _ =  create_radial_lines(brain_param, num_lines; view = view)
    
    # Set searchmethod
    sm = searchmethod=KDTreeSearch(num_nearest_vertices=num_nearest_vertices)
    
    ip = Interpolable(ps, searchmethod = sm)
    
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
        # rel_diff = (ps - mean_p)/mean_p                # Relative difference to mean
        diff = ps - mean_p               # difference to mean
        sq_diff = Interpolable((diff)*(diff), searchmethod = sm)   # squared relative difference
        append!(var, sum(∫(sq_diff)*dL) / len)          # Variance
        append!(rad_len, len)
    end
    return mean_pos, rad_len, var
end

function compute_nflow(brain_param, us, Γ; degree = 2)
    dΓ = Measure(Γ, degree)
    n̂Γ = get_normal_vector(Γ) 
    nflow = sum(∫(us.⁺ ⋅ n̂Γ.⁺)dΓ) / sum(∫(1)dΓ)
    nflow_sqr = sum(∫((us.⁺ ⋅ n̂Γ.⁺)*(us.⁺ ⋅ n̂Γ.⁺))dΓ) / sum(∫(1)dΓ)
    return nflow, nflow_sqr
end
  

function solution_convergence_vs_lc(brain_param, PDE_param, start_lc, end_lc, num_samples; logrange = true, filename = "sol_conv.txt")
  
    title = "2D brain simulation: l²-norm difference between final lc and previous lc respectively"
    brain_param.lc = NaN
    outfile = create_file(path * "txt_files/", filename, title, brain_param, PDE_param)

    lc = logrange ?  10 .^(range(log10(start_lc),stop=log10(end_lc),length=num_samples)) : LinRange(start_lc, end_lc, num_samples)

    write(outfile, "#\n# --- Sampling --- #\n")
    @printf(outfile,"# num_samples = %d\n", num_samples)
    @printf(outfile, "# lc = [")
    [@printf(outfile,"%.3e, ", lc[1 + length(lc) - i]) for i in 1:length(lc)-1]
    @printf(outfile, "%.3e]\n", lc[1])

    write(outfile, "#\n# --- Data --- #\n")
    write(outfile, "lc, l²_Δus, l²_Δps, l²_Δpd\n")
    close(outfile)

    brain_param.lc = last(lc)
    model, pgs_dict = create_brain(brain_param)
    ref_us, ref_ps, ref_pd, ref_ΩS, ref_ΩD, ref_Γ  = brain_PDE(model, pgs_dict, PDE_param)
    add_sample_to_file(path * "txt_files/", filename, [last(lc), 0.0, 0.0, 0.0], 0)
  
    # Convergence
    for i in 2:num_samples
        idx = 1 + num_samples - i
        @printf("i = %d/%d | lc = %e\n", i-1, num_samples-1, lc[i])
        brain_param.lc = lc[idx]
        model, pgs_dict = create_brain(brain_param)
        us, ps, pd, ΩS, ΩD, Γ = brain_PDE(model, pgs_dict, PDE_param)
        l2_diff = calc_solution_diff(brain_param, [us, ps, pd], [ref_us, ref_ps, ref_pd])
        add_sample_to_file(path * "txt_files/", filename, [lc[idx], l2_diff...], 0)
    end
    close(outfile)
end


function eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines; folder_name = nothing, p_filename = "ps_radial_var.txt", nflow_filename = "us_nflow.txt")
    @assert start_width <= brain_param.r_brain ["start width must be less than radial brain size"]
    @assert num_samples > 1  ["must have at least 2 sample points"]
    
    # Create folder for data
    folder_name = isnothing(folder_name) ?  "data_" * Dates.format(now(), "d_m_HH_MM") * "/" : "data_" * folder_name * "/"
    mkdir(path * folder_name)
    mkdir(path *folder_name * "txt_files")
    mkdir(path *folder_name * "vtu_files")
    mkdir(path *folder_name * "png_files")

    
    brain_param.d_ratio = NaN
    width = LinRange(start_width, end_width, num_samples)
    
    
    # --- Pressure variance --- #
    title = "2D brain simulation: Stokes pressure variance in radial direction as a function of fluid path width"
    outfile = create_file(path * folder_name * "txt_files/", p_filename, title, brain_param, PDE_param)

    write(outfile, "#\n# --- Sampling --- #\n")
    @printf(outfile,"# num_samples = %d\n", num_samples)
    @printf(outfile,"# num_rad_lines = %d\n", num_rad_lines)
    @printf(outfile, "# width = [")
    [@printf(outfile,"%.3e, ", width[i]) for i in 1:length(width)-1]
    @printf(outfile, "%.3e]\n", last(width))
    write(outfile, "#\n# --- Data --- #\n")
    write(outfile, "mean_pos_x, mean_pos_y, rad_len, variance\n")
    close(outfile)

    # --- Normal flow on interface --- #
    title = "2D brain simulation: Stokes normal velocity on interface as a function of fluid path width"
    outfile = create_file(path * folder_name * "txt_files/", nflow_filename, title, brain_param, PDE_param)
    write(outfile, "#\n# --- Sampling --- #\n")
    @printf(outfile,"# num_samples = %d\n", num_samples)
    @printf(outfile, "# width = [")
    [@printf(outfile,"%.3e, ", width[i]) for i in 1:length(width)-1]
    @printf(outfile, "%.3e]\n", last(width))
    write(outfile, "#\n# --- Data --- #\n")
    write(outfile, "width, u×n̂, (u×n̂)²\n")
    close(outfile)
    

    for i in 1:num_samples
        @printf("i = %d/%d | width = %e\n", i, num_samples, width[i])

        # Create geometry
        brain_param.d_ratio = width[i]/brain_param.r_brain  
        model, pgs_dict = create_brain(brain_param)

        # Solve PDE
        us, ps, pd, ΩS, ΩD, Γ = brain_PDE(model, pgs_dict, PDE_param; write = (path * folder_name * "vtu_files/", @sprintf("%.2e", width[i]))) # or write = false
        
        # Perform evaluations 
        mean_pos, rad_len, var = evaluate_radial_var(brain_param, ps, num_rad_lines)
        nflow, nflow_sqr = compute_nflow(brain_param, us, Γ)
        plot_nflow_profile(brain_param, us, Γ, path * folder_name * "png_files/",  @sprintf("%.2e", width[i])) 
        
        # Write to file
        p_data = [getindex.(mean_pos, 1), getindex.(mean_pos, 2), rad_len, var]
        nflow_data = [width[i], nflow, nflow_sqr]
        add_sample_to_file(path * folder_name * "txt_files/", p_filename, p_data, num_rad_lines)
        add_sample_to_file(path * folder_name * "txt_files/", nflow_filename, nflow_data, 0)

    end
    
   

end







# # --- Brain Model [length unit: meter] --- # 
# lc = 1e-3 
# arcLen = (100e-3, 0)
# r_brain = 10e-3  
# d_ratio = 1.5e-3/r_brain
# r_curv = 50e-3 
# inner_perturb = "(x,z) -> 0.3e-3 * cos(pi * abs(x) / 2e-3) "  # 10 - 20 waves pr. cm is a good number
# outer_perturb = "(x,z) -> 0.0 "  # 10 - 20 waves pr. cm is a good number
# BS_points = (1000, 0) 
# field_Lc_lim = [1 / 2, 1]
# field_Dist_lim = [1e-3, 5e-3] 
# brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)

# # --- PDE parameters --- #
# μ = 0.8e-3  # Cerobrospinal fluid viscosity [Pa * s]
# Κ = 1e-16   # Permeability in brain parenchyma [m^2] 
# α = "(x) -> 1*μ/sqrt(Κ)" # Slip factor on Γ [Pa * s / m]
# ps0 = "(x) -> x[1] < 0 ? 1*133.3224 : 0." # 1*mmHg [Pa]
# ∇pd0 = "(x) -> VectorValue(0.0, 0.0)" # Zero flux
# PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)


# model, pgs_dict = create_brain(brain_param; view=false, write=false)
# ush, psh, pdh, ΩS, ΩD, Γ = brain_PDE(model, pgs_dict, PDE_param; write = (path * "data_testing/", @sprintf("%.2e", lc)))
# ush, psh, pdh, ΩS, ΩD, Γ = brain_PDE(model, pgs_dict, PDE_param; write = false)



# # --- Evaluations --- #



# start_width = 5e-3
# end_width = 0.5e-3
# num_samples = 5
# num_rad_lines = 100
# eval_decreasing_lc(brain_param, PDE_param, start_width, end_width, num_samples, num_rad_lines)



# start_lc = 1e-2
# end_lc = 1e-3
# num_samples = 3
# solution_convergence_vs_lc(brain_param, start_lc, end_lc, num_samples)



# mean_pressure_vs_lc(brain_param, start_lc, end_lc, num_samples, logrange = true)


