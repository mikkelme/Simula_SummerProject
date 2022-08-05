using Dates
using Printf


include("../mesh_generators/create_brain.jl")
include("./brain_sim_3D.jl")
include("./evaluate.jl")



path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mikkelme/Simula_SummerProject/brain_simulations/" # SSH
end




function eval_decreasing_lc_3D(brain_param, PDE_param, start_width, end_width, num_samples; folder_name = nothing, nflow_filename = "us_nflow.txt")
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
    

    # --- Normal flow on interface --- #
    title = "3D brain simulation: Stokes normal velocity on interface as a function of fluid path width"
    outfile = create_file(path * folder_name * "txt_files/", nflow_filename, title, brain_param, PDE_param)
    write(outfile, "#\n# --- Sampling --- #\n")
    @printf(outfile,"# num_samples = %d\n", num_samples)
    @printf(outfile, "# width = [")
    [@printf(outfile,"%.3e, ", width[i]) for i in 1:length(width)-1]
    @printf(outfile, "%.3e]\n", last(width))
    write(outfile, "#\n# --- Data --- #\n")
    write(outfile, "width, u⋅n̂, (u⋅n̂)², u⋅u \n")
    close(outfile)
    

    for i in 1:num_samples
        @printf("i = %d/%d | width = %e\n", i, num_samples, width[i])

        # Create geometry
        brain_param.d_ratio = width[i]/brain_param.r_brain  
        model, pgs_dict = create_brain(brain_param)

        # Solve PDE
        us, ps, pd, ΩS, ΩD, Γ = brain_PDE_3D(model, pgs_dict, PDE_param; write = (path * folder_name * "vtu_files/", @sprintf("%.2e", width[i]))) # or write = false
        
        # Perform evaluations 
        nflow, nflow_sqr, flow_sqr = compute_nflow(brain_param, us, Γ)
        
        # Write to file
        nflow_data = [width[i], nflow, nflow_sqr, flow_sqr]
        add_sample_to_file(path * folder_name * "txt_files/", nflow_filename, nflow_data, 0)

    end
    
   

end





