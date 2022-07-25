using Gridap
using GridapGmsh
using GridapGmsh: gmsh
include("./brain_mesh_utils.jl")


# function cl_perturbed_arc(param::model_params, domain; safety = 0)
#     pointTags = []
#     rI = param.r_curv - param.r_brain                   # Inner radius  
#     rD = param.r_curv - param.d_ratio * param.r_brain   # Radius for dividing line
#     arcLen = param.arcLen[1]                            # Length of outer arc
#     angle = arcLen / param.r_curv                       # Angle span y -> x -axis


#     # if domain == "S"
#     #     println("#---Stokes---#")
#     #     # r_cl = param.r_curv - (param.r_curv - rD)/2
#     #     r_cl = param.r_curv - (param.d_ratio * param.r_brain)/2 
        
#     # elseif domain == "D"
#     #     println("#---Darcy---#")
#     #     r_cl = rI + (rD-rI)/2
        
#     # else 
#     #     @printf("Domain: %s, is not defined in create_centerline", domain)
#     #     return
#     # end

#     r_SCL = param.r_curv - (param.d_ratio * param.r_brain)/2 
#     r_DCL = rI + (rD-rI)/2
    
#     Stokes_perturbation_func(x,z) = (param.inner_perturb(x,z) + param.outer_perturb(x,z))/2
#     Darcy_perturbation_func(x,z) = param.inner_perturb(x,z)/2
    
#     for i in LinRange(0 + safety, 1 - safety, param.BS_points[1])
#         phi = - angle / 2 + angle * i
#         x_arcLen = r_cl * phi
#         r_pert = domain == "S" ?  r_SCL + Stokes_perturbation_func(x_arcLen, 0) : r_DCL + Darcy_perturbation_func(x_arcLen, 0)
#         append!(pointTags, gmsh.model.occ.addPoint(spherical_to_cartesian(r_pert, 0, phi)...))
#     end
  
#      vertex = [pointTags[1], last(pointTags)]
#      BSpline = gmsh.model.occ.addBSpline(pointTags)
#      gmsh.model.occ.synchronize()
  
#      return vertex, BSpline
    



# end


function cl_perturbed_arc(param::model_params, domain)
    pointTags = []
    rI = param.r_curv - param.r_brain                   # Inner radius  
    rD = param.r_curv - param.d_ratio * param.r_brain   # Radius for dividing line
    arcLen = param.arcLen[1]                            # Length of outer arc
    angle = arcLen / param.r_curv                       # Angle span y -> x -axis

    if domain == "S"
        r_cl = param.r_curv - (param.r_curv - rD)/2
    elseif domain == "D"
        r_cl = rI + (rD-rI)/2 
    else 
        @printf("Domain: %s, is not defined in create_centerline", domain)
        return
    end

        perturbation_func(x,z) = domain == "S" ? (param.inner_perturb(x,z) + param.outer_perturb(x,z))/2 : (param.inner_perturb(x,z))/2
    
    for i in LinRange(0, 1, param.BS_points[1])
        phi = - angle / 2 + angle * i
        x_arcLen = r_cl * phi
        r_pert = r_cl + perturbation_func(x_arcLen, 0)
        append!(pointTags, gmsh.model.occ.addPoint(spherical_to_cartesian(r_pert, 0, phi)...))
    end

    vertex = [pointTags[1], last(pointTags)]
    BSpline = gmsh.model.occ.addBSpline(pointTags)
    gmsh.model.occ.synchronize()

    return vertex, BSpline
end


function create_centerline(param::model_params, domain = "S"; view = false)
    gmsh.initialize(["", "-clmax", string(param.lc)])
    gmsh.option.setNumber("Mesh.SaveAll", 1)  # For direct wiring

    vertex, BSpline = cl_perturbed_arc(param, domain)
    gmsh.model.addPhysicalGroup(0, [vertex[1]]) # center line (left point)
    gmsh.model.addPhysicalGroup(1, [BSpline])   # center line (arc)
    gmsh.model.addPhysicalGroup(0, [vertex[2]]) # center line (right point)
    gmsh.model.mesh.generate(1)

    # @show gmsh.model.mesh.getNodes()
    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

   
    model, pgs_dict = direct_wiring(gmsh)
    gmsh.finalize()

    return model, pgs_dict




end



# pert_func(x,z) = 0.0
# lc = 0.1
# arcLen = 5
# r_brain = 2
# d_ratio = 0.5
# r_curv = 50
# BS_points = 10

