using Gridap
using GridapGmsh
using Printf
using Plots




path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/Solvers/"

include("./unit_box_direct.jl")



function poisson_solver(model, f0, g0, h0)
    # (8) 3 (7)
    #  4     2
    # (5) 1 (6)

    labels = get_face_labeling(model)
    plus_tags = get_tags_from_group(path * "unit_box.msh", [1, 2, 3, 4, 5, 6, 8])
    minus_tags = get_tags_from_group(path * "unit_box.msh", [8])

    add_tag_from_tags!(labels, "plus", plus_tags)
    add_tag_from_tags!(labels, "minus", minus_tags)



    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)

    V = TestFESpace(model, reffe, labels=labels, conformity=:H1, dirichlet_tags=["plus", "minus"])
    U = TrialFESpace(V, [5, -5])

    # Integration 
    degree = 2
    Ω = Triangulation(model) # Integration mesh of the domain Omega
    dΩ = Measure(Ω, degree) # Gauss-like quadrature (how to connect point in numerical integration)



    # --- Weak form --- #
    a(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ
    b(v) = ∫(v * f0) * dΩ

    # --- Solve --- #
    op = AffineFEOperator(a, b, U, V)
    ls = LUSolver()
    solver = LinearFESolver(ls)
    uh = solve(solver, op)


    writevtk(Ω, path * "poisson_solver_results", cellfields=["uh" => uh])

end


# function my_GmshDiscreteModel(mshfile; renumber=true)
#     # @check_if_loaded
#     if !isfile(mshfile)
#         error("Msh file not found: $mshfile")
#     end

#     gmsh.initialize()
#     gmsh.option.setNumber("General.Terminal", 1)
#     gmsh.option.setNumber("Mesh.SaveAll", 1)
#     gmsh.option.setNumber("Mesh.MedImportGroupsOfNodes", 1)
#     gmsh.open(mshfile)

#     renumber && gmsh.model.mesh.renumberNodes()
#     renumber && gmsh.model.mesh.renumberElements()

#     Dc = 2
#     Dp = 2

#     node_to_coords = _setup_node_coords(gmsh, Dp)
#     # nnodes = length(node_to_coords)
#     # vertex_to_node, node_to_vertex = _setup_nodes_and_vertices(gmsh, node_to_coords)
#     # grid, cell_to_entity = _setup_grid(gmsh, Dc, Dp, node_to_coords, node_to_vertex)
#     # cell_to_vertices = _setup_cell_to_vertices(get_cell_node_ids(grid), node_to_vertex, nnodes)
#     # grid_topology = UnstructuredGridTopology(grid, cell_to_vertices, vertex_to_node)
#     # labeling = _setup_labeling(gmsh, grid, grid_topology, cell_to_entity, vertex_to_node, node_to_vertex)
#     gmsh.finalize()

#     # UnstructuredDiscreteModel(grid, grid_topology, labeling)
# end





function get_tags_from_group(mshfile, pgs_tags)
    gmsh.initialize()
    gmsh.open(mshfile)
    pgs = gmsh.model.getPhysicalGroups()
    tags = Vector{Int64}()

    for i in 1:length(pgs_tags)
        pgs_tag = pgs_tags[i]
        for j in 1:length(pgs)

            if pgs[j][2] == pgs_tag
                append!(tags, j)
            end
        end
        length(tags) < i && @printf("Physical group tag: %i, not found\n", pgs_tag)
        length(tags) > i && @printf("Multiple candidates found for physical group tag: %i\n", pgs_tag)
        @assert(length(tags) == i)

    end

    gmsh.finalize()


    return tags

end


###############
# MS
# u0(x) = cos(π * x[1] * x[2])
# f0(x) = π^2 * (x[1]^2 + x[2]^2) * cos(π * x[1] * x[2])
h0(x) = VectorValue(-π * x[2] * sin(π * x[1] * x[2]), -π * x[1] * sin(π * x[1] * x[2]))

f0(x) = 0
g0(x) = 5


# create_unit_box(2, false)
model = GmshDiscreteModel(path * "unit_box.msh")

const D3 = 3
const POINT = 15
const UNSET = 0


function wiring(mshfile; renumber=true)
    @check_if_loaded
    if !isfile(mshfile)
        error("Msh file not found: $mshfile")
    end

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.option.setNumber("Mesh.MedImportGroupsOfNodes", 1)
    gmsh.open(mshfile)

    renumber && gmsh.model.mesh.renumberNodes()
    renumber && gmsh.model.mesh.renumberElements()

    Dc = _setup_cell_dim(gmsh)
    Dp = _setup_point_dim(gmsh, Dc)
    node_to_coords = _setup_node_coords(gmsh, Dp)
    nnodes = length(node_to_coords)
    vertex_to_node, node_to_vertex = _setup_nodes_and_vertices(gmsh, node_to_coords)
    grid, cell_to_entity = _setup_grid(gmsh, Dc, Dp, node_to_coords, node_to_vertex)
    cell_to_vertices = _setup_cell_to_vertices(get_cell_node_ids(grid), node_to_vertex, nnodes)
    grid_topology = UnstructuredGridTopology(grid, cell_to_vertices, vertex_to_node)
    labeling = _setup_labeling(gmsh, grid, grid_topology, cell_to_entity, vertex_to_node, node_to_vertex)
    gmsh.finalize()

    UnstructuredDiscreteModel(grid, grid_topology, labeling)
end


wiring(path * "unit_box.msh")

# tags = get_tags_from_group(path * "unit_box.msh", [5, 6, 7, 8])
# println(tags)

# gmsh.initialize()
# gmsh.open(path * "unit_box.msh")
# dim_to_group_to_name = _setup_dim_to_group_to_name(gmsh)
# pgs = gmsh.model.getPhysicalGroups()
# println(pgs)
# # gmsh.fltk.initialize()
# # gmsh.fltk.run()
# gmsh.finalize()
# my_GmshDiscreteModel(path * "unit_box.msh")





# poisson_solver(model, f0, g0, h0)