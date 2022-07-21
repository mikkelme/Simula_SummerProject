using Gridap
using GridapGmsh
using GridapGmsh: gmsh
include("./brain_mesh_utils.jl")


function cl_perturbed_arc(param::model_params, distort = 0)
    pointTags = []
   
    arcLen = param.arcLen[1]  
    r_cl = param.r_curv - (param.d_ratio * param.r_brain)/2 + distort
    angle = arcLen / param.r_curv   
    @show angle*r_cl
    perturbation_func(x,z) = (param.inner_perturb(x,z) + param.outer_perturb(x,z))/2
  
    angle = arcLen / param.r_curv  
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

function create_centerline(param::model_params; view = false)
    gmsh.initialize(["", "-clmax", string(param.lc)])
    gmsh.option.setNumber("Mesh.SaveAll", 1)  # For direct wiring

    vertex, BSpline = cl_perturbed_arc(param)
    gmsh.model.addPhysicalGroup(0, [vertex[1]]) # center line (left point)
    gmsh.model.addPhysicalGroup(1, [BSpline])   # center line (arc)
    gmsh.model.addPhysicalGroup(0, [vertex[2]]) # center line (right point)


    vertex, BSpline = cl_perturbed_arc(param, 1)
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

