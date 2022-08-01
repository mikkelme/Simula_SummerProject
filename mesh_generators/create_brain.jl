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


function create_brain(param::model_params; view=false, write=false)
    gmsh.initialize(["", "-clmax", string(param.lc)])
    # gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 1)  # For direct wiring
    # gmsh.option.setNumber("Mesh.MedImportGroupsOfNodes", 1)
    gmsh.model.add("brain")
    param.arcLen[2] == 0 ? create_brain_2D(param) : create_brain_3D(param)

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

    write!=false && gmsh.write(write)    

    model, pgs_dict = direct_wiring(gmsh)
    gmsh.finalize()
    # model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/mesh_generators/brain.msh") # Test the mesh
    return model, pgs_dict
    

end


if abspath(PROGRAM_FILE) == @__FILE__
    lc = 2e-4
    arcLen = (100e-3, 0)
    r_brain = 10e-3  
    d_ratio = 1.5e-3/r_brain
    r_curv = 50e-3 
    A = 1e-3
    λ = 10*1e-3
    ω(λ) = 2*pi/λ      

    inner_perturb = @sprintf("(x,z) -> %f * sin(abs(x + z) * %f - pi/2) * fld(mod2pi(abs(x + z) * %f - pi/2),pi) ", A , ω(λ), ω(λ))
    outer_perturb = "(x,z) -> 0.0"  

    # outer_perturb = @sprintf("(x,z) -> %f * cos(abs(x) * %f)", 0.5e-3 , ω(3.0))

    BS_points = (1000, 1000) 
    field_Lc_lim = [1 / 2, 1]
    field_Dist_lim = [1e-3, 5e-3] 


    # 2D brain example
    param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    create_brain(param; view=true, write = false)


    # 3D brain example (Work in progress)
    # arcLen = (100e-3, 20e-3)
    # param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
    # create_brain(param; view=true, write=false)
end





