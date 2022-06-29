using GridapGmsh: gmsh


gmsh.initialize(["", "-clmax", "1"])

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
    # append!(sides, line[i+1])
end
# boundary = gmsh.model.addPhysicalGroup(1, sides)
# gmsh.model.setPhysicalName(1, boundary, "sides")


loop = gmsh.model.occ.addCurveLoop([line[1], line[2], line[3], line[4]])
surf = gmsh.model.occ.addPlaneSurface([loop])
gmsh.model.occ.synchronize()

surf_group = gmsh.model.addPhysicalGroup(2, [surf])
# gmsh.model.setPhysicalName(2, surf_group, "surface")


gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)

gmsh.fltk.initialize()
gmsh.fltk.run()

gmsh.write("dummy.msh")
gmsh.finalize()