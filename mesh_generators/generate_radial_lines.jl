using Gridap
using GridapGmsh
using GridapGmsh: gmsh
include("./brain_mesh_utils.jl")


function create_perturbed_arc(angle, r, perturbation_func, BS_points)
    pointTags = []

    for i in LinRange(0, 1, BS_points)
        phi = - angle / 2 + angle * i
        x_arcLen = r * phi
        r_pert = r + perturbation_func(x_arcLen, 0)
        append!(pointTags, gmsh.model.occ.addPoint(spherical_to_cartesian(r_pert, 0, phi)...))
    end

    vertex = [pointTags[1], last(pointTags)]
    BSpline = gmsh.model.occ.addBSpline(pointTags)
    gmsh.model.occ.synchronize()
    return vertex, BSpline
end

function get_perturbed_point(r, angle, perturbation_func)      
    x_arcLen = r * angle # x-arclen for circle-line coordinate 
    r_pert = r + perturbation_func(x_arcLen, 0) # perturbed radius
    gmsh.model.occ.addPoint(spherical_to_cartesian(r_pert, 0, angle)...)
end


function perturbed_radial_line(param::model_params, angle)
    # Dividing line radius
    rD = param.r_curv - param.d_ratio * param.r_brain 

    # Inner (P1) and outer (P2) point
    P1 = get_perturbed_point(rD, angle, param.inner_perturb) 
    P2 = get_perturbed_point(param.r_curv, angle, param.outer_perturb) 
   
    # Connect with line
    line = gmsh.model.occ.addLine(P1, P2)
    return line

end




function create_radial_lines(param::model_params, num_lines; view = false)
    gmsh.initialize(["", "-clmax", string(param.lc)])
    gmsh.option.setNumber("Mesh.SaveAll", 1)  # For direct wiring

    # Get info
    arcLen = param.arcLen[1]  
    θ = arcLen / param.r_curv   
    Δθ =  θ/(num_lines+1) 

    # Create lines 
    lines = [perturbed_radial_line(param, Δθ*i - θ/2) for i in 1:num_lines] # Make choice
    gmsh.model.occ.synchronize()

    # Add physical groups: left| 1 → num_points |right
    [gmsh.model.addPhysicalGroup(1, [lines[i]]) for i in 1:num_lines]
    gmsh.model.mesh.generate(1)


    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end


    model, pgs_dict = direct_wiring(gmsh)
    gmsh.finalize()

    return model, pgs_dict



end



# # --- Brain Model --- # 
# lc = 0.1
# arcLen = (5, 0)
# r_brain = 2
# d_ratio = 0.5
# # r_curv = 10
# inner_perturb(x, y) = 0.2 * cos(pi * abs(x) / 0.5) 
# outer_perturb(x, y) = 0.2 * cos(pi * abs(x) / 2)  

# # inner_perturb(x, y) = 0
# # outer_perturb(x, y) = 0
# BS_points = (arcLen[1]*20, arcLen[2]*10)
# field_Lc_lim = [1 / 2, 1]
# field_Dist_lim = [0.1, 0.5]
# brain_params = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)


# num_lines = 5
# create_radial_lines(brain_params, num_lines, view = true)