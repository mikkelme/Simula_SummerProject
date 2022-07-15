include("./create_brain_3D.jl")
include("./create_brain_2D.jl")



function create_brain(param::model_params, view=true, write=false)
    gmsh.initialize(["", "-clmax", string(param.lc)])
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

    if write
        gmsh.write("brain.msh")
    end 

    gmsh.finalize()



end



if abspath(PROGRAM_FILE) == @__FILE__

    lc = 0.5
    arcLen = (5, 0)
    r_brain = 2
    d_ratio = 0.5
    r_curv = 50
    inner_perturb(x, y) = 0.2 * cos(pi * abs(x) / 0.5) + 0.2 * cos(pi * abs(y) / 0.5)
    outer_perturb(x, y) = 0.2 * cos(pi * abs(x) / 2)  + 0.2 * cos(pi * abs(y) / 1)
    BS_points = (arcLen[1]*20, arcLen[2]*10)
    field_Lc_lim = [1 / 2, 1]
    field_Dist_lim = [0.1, 0.5]

    param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)


    create_brain(param, true, true)
end