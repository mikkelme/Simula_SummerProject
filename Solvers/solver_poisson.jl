using Gridap
using GridapGmsh
using Printf
using Plots




path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/Solvers/"

include("./unit_box_direct.jl")



function poisson_solver(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags)

    labels = get_face_labeling(model)
    neumann_conditions = !isempty(neumann_tags)

    # all_tags = collect(keys(pgs_dict)) # Only works with no irrelevant tags
    # dirichlet_tags = [pgs_dict[tag] for tag in dirichlet_tags] # is going to be dircihlet 
    # neumann_tags = filter(x -> x ∉ dirichlet_tags, collect(keys(pgs_dict))) # is going to be neumann

    dirichlet_tags = [pgs_dict[tag] for tag in dirichlet_tags]
    neumann_tags = [pgs_dict[tag] for tag in neumann_tags]


    add_tag_from_tags!(labels, "diri", dirichlet_tags)




    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)

    V = TestFESpace(model, reffe, labels=labels, conformity=:H1, dirichlet_tags=["diri"])
    U = TrialFESpace(V, g0)

    # Integration 
    degree = 2
    Ω = Triangulation(model) # Integration mesh of the domain Omega
    dΩ = Measure(Ω, degree) # Gauss-like quadrature (how to connect point in numerical integration)

    if neumann_conditions
        Γ = [BoundaryTriangulation(model, tags=tag) for tag in neumann_tags]
        ν = [get_normal_vector(Γ[i]) for i in 1:length(neumann_tags)]
        dΓ = [Measure(Γ[i], degree) for i in 1:length(neumann_tags)]
        # h = [neumann[tag] for tag in neumann_tags]

   

    end


    # --- Weak form --- #
    a(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ
    b(v) = neumann_conditions ? ∫(v * f0) * dΩ + sum([∫(v * (h0 ⋅ ν[i])) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v * f0) * dΩ

    # b(v) = ∫(v * f0) * dΩ

    # --- Solve --- #
    op = AffineFEOperator(a, b, U, V)
    ls = LUSolver()
    solver = LinearFESolver(ls)
    uh = solve(solver, op)


    writevtk(Ω, path * "poisson_solver_results", cellfields=["uh" => uh])

end





###############
# MS
u0(x) = cos(π * x[1] * x[2])
f0(x) = π^2 * (x[1]^2 + x[2]^2) * cos(π * x[1] * x[2])
h0(x) = VectorValue(-π * x[2] * sin(π * x[1] * x[2]), -π * x[1] * sin(π * x[1] * x[2]))


# f0(x) = 0
# g0(x) = 5


model, pgs_dict = create_unit_box(0.01, false)
# create_unit_box(2, false)
# model = GmshDiscreteModel(path * "unit_box_new.msh")

all_btags = [1, 2, 3, 4, 5, 6, 7, 8]
dirichlet_tags = [1, 2, 4, 5, 6, 7, 8]
neumann_tags = filter(x -> x ∉ dirichlet_tags, collect(keys(pgs_dict))) 




#-------> Do convergence test tomorrow and check that everything is good <----------#
poisson_solver(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags)





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







# function get_tags_from_group(mshfile, pgs_tags)
#     gmsh.initialize()
#     gmsh.open(mshfile)
#     pgs = gmsh.model.getPhysicalGroups()
#     tags = Vector{Int64}()

#     for i in 1:length(pgs_tags)
#         pgs_tag = pgs_tags[i]
#         for j in 1:length(pgs)

#             if pgs[j][2] == pgs_tag
#                 append!(tags, j)
#             end
#         end
#         length(tags) < i && @printf("Physical group tag: %i, not found\n", pgs_tag)
#         length(tags) > i && @printf("Multiple candidates found for physical group tag: %i\n", pgs_tag)
#         @assert(length(tags) == i)

#     end

#     gmsh.finalize()


#     return tags

# end