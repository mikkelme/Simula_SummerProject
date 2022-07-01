using GridapGmsh: gmsh


function create_unit_box(lc, view=false)
    path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/"

    gmsh.initialize(["", "-clscale", string(lc)])


    A = gmsh.model.occ.addPoint(0, 0, 0)
    B = gmsh.model.occ.addPoint(1, 0, 0)
    C = gmsh.model.occ.addPoint(1, 1, 0)
    D = gmsh.model.occ.addPoint(0, 1, 0)

    P = [A, B, C, D]
    line = []
    # sides = []
    for i in 0:3
        # println(1 + i % 4 , " ", 1 + (i + 1) % 4)
        append!(line, gmsh.model.occ.addLine(1 + i % 4, 1 + (i + 1) % 4))
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(1, [line[i+1]])
    end


    loop = gmsh.model.occ.addCurveLoop([line[1], line[2], line[3], line[4]])
    surf = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setTransfiniteSurface(surf)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    surf_group = gmsh.model.addPhysicalGroup(2, [surf])


    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end
    gmsh.write(path * "unit_box.msh")
    gmsh.finalize()

end


if abspath(PROGRAM_FILE) == @__FILE__

    lc = 0.1
    lc = 2
    create_unit_box(lc, true)
end


