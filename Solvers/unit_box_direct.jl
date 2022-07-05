using Gridap
using GridapGmsh
using GridapGmsh: gmsh


include("./GmshDiscreteModels.jl")
path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/Solvers/"

function create_unit_box(lc, view=false)
    gmsh.initialize(["", "-clscale", string(lc)])


    A = gmsh.model.occ.addPoint(0, 0, 0)
    B = gmsh.model.occ.addPoint(1, 0, 0)
    C = gmsh.model.occ.addPoint(1, 1, 0)
    D = gmsh.model.occ.addPoint(0, 1, 0)

    P = [A, B, C, D]
    line = []
    # sides = []
    for i in 1:4
        from, to = P[1+(i-1)%4], P[1+i%4]
        append!(line, gmsh.model.occ.addLine(from, to))

        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(1, [line[i]], i)
    end
    # gmsh.model.setPhysicalName(1, 11, "line")


    # gmsh.model.occ.synchronize()

    # Mark vertices too
    gmsh.model.addPhysicalGroup(0, [A], 5)
    gmsh.model.addPhysicalGroup(0, [B], 6)
    gmsh.model.addPhysicalGroup(0, [C], 7)
    gmsh.model.addPhysicalGroup(0, [D], 8)

    # gmsh.model.occ.synchronize()
    # gmsh.model.setPhysicalName(0, 80, "UR")
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

    direct_wiring(gmsh)

    # gmsh.finalize()
end



function direct_wiring(gmsh; renumber=true)
    renumber && gmsh.model.mesh.renumberNodes()
    renumber && gmsh.model.mesh.renumberElements()

    Dc = GridapGmsh._setup_cell_dim(gmsh)
    Dp = GridapGmsh._setup_point_dim(gmsh, Dc)
    node_to_coords = GridapGmsh._setup_node_coords(gmsh, Dp)
    nnodes = length(node_to_coords)
    vertex_to_node, node_to_vertex = _setup_nodes_and_vertices(gmsh, node_to_coords)
    grid, cell_to_entity = _setup_grid(gmsh, Dc, Dp, node_to_coords, node_to_vertex)
    cell_to_vertices = _setup_cell_to_vertices(Gridap.Geometry.get_cell_node_ids(grid), node_to_vertex, nnodes)
    grid_topology = Gridap.Geometry.UnstructuredGridTopology(grid, cell_to_vertices, vertex_to_node)
    labeling = GridapGmsh._setup_labeling(gmsh, grid, grid_topology, cell_to_entity, vertex_to_node, node_to_vertex)
    gmsh.finalize()

    Gridap.Geometry.UnstructuredDiscreteModel(grid, grid_topology, labeling)
end

# Gridapp.Geometryget_cell_node_ids


# if abspath(PROGRAM_FILE) == @__FILE__
#     lc = 0.1
#     lc = 2
#     create_unit_box(lc, true)
# end



lc = 0.1
lc = 2
create_unit_box(lc, false)