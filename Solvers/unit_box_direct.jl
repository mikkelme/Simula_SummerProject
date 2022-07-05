using Gridap
using GridapGmsh
using GridapGmsh: gmsh


include("./DiscreteModel_utils.jl")


function create_unit_box(lc, view=false)
    # Physical tags
    # (8) 3 (7)
    #  4     2
    # (5) 1 (6)

    gmsh.initialize(["", "-clscale", string(lc)])

    A = gmsh.model.occ.addPoint(0, 0, 0)
    B = gmsh.model.occ.addPoint(1, 0, 0)
    C = gmsh.model.occ.addPoint(1, 1, 0)
    D = gmsh.model.occ.addPoint(0, 1, 0)

    P = [A, B, C, D]
    line = []

    for i in 1:4
        from, to = P[1+(i-1)%4], P[1+i%4]
        append!(line, gmsh.model.occ.addLine(from, to))

        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(1, [line[i]], i)
    end


    # Mark vertices too
    gmsh.model.addPhysicalGroup(0, [A], 5)
    gmsh.model.addPhysicalGroup(0, [B], 6)
    gmsh.model.addPhysicalGroup(0, [C], 7)
    gmsh.model.addPhysicalGroup(0, [D], 8)

    loop = gmsh.model.occ.addCurveLoop([line[1], line[2], line[3], line[4]])
    surf = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setTransfiniteSurface(surf)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    # surf_group = gmsh.model.addPhysicalGroup(2, [surf]) # This must be here for some reason when writing



    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    
    
    model, pgs_dict = direct_wiring(gmsh)
    gmsh.finalize()

    return model, pgs_dict
end



if abspath(PROGRAM_FILE) == @__FILE__
    lc = 0.1
    lc = 2
    create_unit_box(lc, true)
end



