using Gridap
using GridapGmsh
using GridapGmsh: gmsh
using Printf
include("./create_brain_3D.jl")
include("./create_brain_2D.jl")
include("./DiscreteModel_utils.jl")

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/mesh_generators/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/mesh_generators/"
end


function create_brain(param::model_params; view=true, write=false)
    gmsh.initialize(["", "-clmax", string(param.lc)])

    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.option.setNumber("Mesh.MedImportGroupsOfNodes", 1)

    gmsh.model.add("brain")
    param.arcLen[2] == 0 ? create_brain_2D(param) : create_brain_3D(param)

    # View and finalize
    if view
        # Shift rotation center 
        gmsh.option.setNumber("General.RotationCenterGravity", 0)
        gmsh.option.setNumber("General.RotationCenterZ", param.r_curv-param.r_brain/2)
        
        # Rotate view
        gmsh.option.setNumber("General.Trackball", 0)
        gmsh.option.setNumber("General.RotationX", -90) # Used if trackball = 0

        gmsh.fltk.initialize()
        gmsh.option.setNumber("General.Trackball", 1) # Reenable for better rotation control
        gmsh.fltk.run()
    end

    write && gmsh.write(path * "brain.msh")
    
    model, pgs_dict = direct_wiring(gmsh)
    gmsh.finalize()

    return model, pgs_dict


end


if abspath(PROGRAM_FILE) == @__FILE__

    lc = 1
    arcLen = (5, 0)
    r_brain = 2
    d_ratio = 0.5
    r_curv = 50
    inner_perturb(x, y) = 0.2 * cos(pi * abs(x) / 0.5) + 0.2 * cos(pi * abs(y) / 0.5)
    outer_perturb(x, y) = 0.2 * cos(pi * abs(x) / 2)  + 0.2 * cos(pi * abs(y) / 1)
    BS_points = (arcLen[1]*2, arcLen[2]*10)
    field_Lc_lim = [1 / 2, 1]
    field_Dist_lim = [0.1, 0.5]

    param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)


    create_brain(param; view=false, write=false)
end
