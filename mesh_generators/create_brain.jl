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


function pgs_tags(pgs_dict, tags)
    return [pgs_dict[tag] for tag in tags]
end


function create_brain(param::model_params; view=true, write=false)
    gmsh.initialize(["", "-clmax", string(param.lc)])
    gmsh.option.setNumber("Mesh.SaveAll", 1)  # For direct wiring
    gmsh.model.add("brain")
    param.arcLen[2] == 0 ? create_brain_2D_xy(param) : create_brain_3D_xy(param)

    # View and finalize
    if view
        # Change rotation center 
        gmsh.option.setNumber("General.RotationCenterGravity", 0)
        gmsh.option.setNumber("General.RotationCenterY", param.r_curv-param.r_brain/2)
        
        # Shift view
        gmsh.option.setNumber("General.TranslationY", -1/2*(param.r_curv-param.r_brain/2)) 

        # View
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    write && gmsh.write(path * "brain.msh")    
    model, pgs_dict = direct_wiring(gmsh)
    gmsh.finalize()

    return model, pgs_dict


end


if abspath(PROGRAM_FILE) == @__FILE__
    lc = 0.5
    arcLen = (5, 0)
    r_brain = 2
    d_ratio = 0.5
    r_curv = 50
    inner_perturb(x, z) = 0.2 * cos(pi * abs(x) / 0.5) + 0.2 * cos(pi * abs(z) / 0.5)
    outer_perturb(x, z) = 0.2 * cos(pi * abs(x) / 2)  + 0.2 * cos(pi * abs(z) / 1)
    # inner_perturb(x, y) = 0
    # outer_perturb(x, y) = 0
    BS_points = (arcLen[1]*20, arcLen[2]*10)
    field_Lc_lim = [1 / 2, 1]
    field_Dist_lim = [0.1, 0.5]

    param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    create_brain(param; view=true, write=false)
end



