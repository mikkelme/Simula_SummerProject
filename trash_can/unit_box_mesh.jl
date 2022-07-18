using GridapGmsh: gmsh


function create_unit_box(lc, view=false)
    path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/"
    if !ispath(path)
        path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/PDE_stuff/"
    end
    @show path

    gmsh.initialize(["", "-clscale", string(lc)])


    A = gmsh.model.occ.addPoint(0, 0, 0, 0., 20)
    B = gmsh.model.occ.addPoint(1, 0, 0)
    C = gmsh.model.occ.addPoint(1, 1, 0)
    D = gmsh.model.occ.addPoint(0, 1, 0)

    P = [A, B, C, D]
    line = []
    # sides = []
    for i in 0:3
        # println(1 + i % 4 , " ", 1 + (i + 1) % 4)
        l = gmsh.model.occ.addLine(P[1 + i % 4], P[1 + (i + 1) % 4])
        append!(line, l)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(1, [l], 10 + i+1)
    end
    gmsh.model.setPhysicalName(1, 11, "line")

    
    gmsh.model.occ.synchronize()

    # Mark vertices too
    gmsh.model.addPhysicalGroup(0, [A], 50)
    gmsh.model.addPhysicalGroup(0, [B], 60)
    gmsh.model.addPhysicalGroup(0, [C], 70)
    gmsh.model.addPhysicalGroup(0, [D], 80)        
    
    # gmsh.model.occ.synchronize()
    gmsh.model.setPhysicalName(0, 80, "UR")
    # println("this code is up to date (10.19)")


    loop = gmsh.model.occ.addCurveLoop([line[1], line[2], line[3], line[4]])
    surf = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setTransfiniteSurface(surf)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    surf_group = gmsh.model.addPhysicalGroup(2, [surf], 101)


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


