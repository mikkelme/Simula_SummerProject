using Gridap
using GridapGmsh
using GridapGmsh: gmsh
using Printf

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/brain_simulations/"
end

function simple_model()

    view = false
    write = true

    y_shift = 2.0
    gmsh.initialize(["", "-clmax", "0.1"])
    gmsh.model.add("brain")


    # Points
    Origo = gmsh.model.occ.addPoint(0.0, 0.0, 0.0)
    A = gmsh.model.occ.addPoint(-1.0, y_shift + 0.0, 0.0)
    B = gmsh.model.occ.addPoint(1.0, y_shift + 0.0, 0.0)
    C = gmsh.model.occ.addPoint(1.0, y_shift + 0.5, 0.0)
    D = gmsh.model.occ.addPoint(1.0, y_shift + 1.0, 0.0)
    E = gmsh.model.occ.addPoint(-1.0, y_shift + 1.0, 0.0)
    F = gmsh.model.occ.addPoint(-1.0, y_shift + 0.5, 0.0)


    # Lines
    AB = gmsh.model.occ.addLine(A, B)
    BC = gmsh.model.occ.addLine(B, C)
    CD = gmsh.model.occ.addLine(C, D)
    # CF = gmsh.model.occ.addLine(C, F)
    # DE = gmsh.model.occ.addLine(D, E)
    EF = gmsh.model.occ.addLine(E, F)
    FA = gmsh.model.occ.addLine(F, A)

    # Arcs
    CF = gmsh.model.occ.addCircleArc(C, Origo, F)
    DE = gmsh.model.occ.addCircleArc(D, Origo, E)


    Loop1 = gmsh.model.occ.addCurveLoop([AB, BC, CF, FA])
    Loop2 = gmsh.model.occ.addCurveLoop([-CF, CD, DE, EF])
    Surf1 = gmsh.model.occ.addPlaneSurface([Loop1])
    Surf2 = gmsh.model.occ.addPlaneSurface([Loop2])
    gmsh.model.occ.synchronize()


    # Physical groups (similar to unit_box)
    gmsh.model.addPhysicalGroup(1, [AB], 1)
    gmsh.model.addPhysicalGroup(1, [FA], 2)
    gmsh.model.addPhysicalGroup(1, [BC], 3)
    gmsh.model.addPhysicalGroup(1, [CF], 4)
    gmsh.model.addPhysicalGroup(1, [EF], 5)
    gmsh.model.addPhysicalGroup(1, [CD], 6)
    gmsh.model.addPhysicalGroup(1, [DE], 7)

    gmsh.model.addPhysicalGroup(0, [A], 8)
    gmsh.model.addPhysicalGroup(0, [B], 9)
    gmsh.model.addPhysicalGroup(0, [F],10)
    gmsh.model.addPhysicalGroup(0, [C],11)
    gmsh.model.addPhysicalGroup(0, [E], 12)
    gmsh.model.addPhysicalGroup(0, [D], 13)

    gmsh.model.addPhysicalGroup(2, [Surf1], 100)
    gmsh.model.addPhysicalGroup(2, [Surf2], 200)


    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    write && gmsh.write(path * "brain.msh")

    gmsh.finalize()



end



simple_model()