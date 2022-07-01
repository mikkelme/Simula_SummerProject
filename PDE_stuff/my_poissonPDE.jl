using Gridap
using GridapGmsh
using Printf


# Equation to solve
# Find scalar field u such that 
# -Δu = -∇²u = f, f ∈ Ω, 
# u = g, on boundary
# ∇u⋅n = h, 
# u g 

# model = GmshDiscreteModel("dummy.msh")
# writevtk(model, "model")


function poisson(model, f, dirichlet, neumann, MMS=nothing)
    path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/"

    dirichlet_conditions = !isempty(dirichlet)
    neumann_conditions = !isempty(neumann)


    order = 1
    reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, order) # type of FE interpolation

    # Test space and trial space
    if dirichlet_conditions

        labels = get_face_labeling(model)

        dirichlet_tags = collect(keys(dirichlet))
        dirichlet_names = Vector{String}()
        for tags in dirichlet_tags
            name = "$tags"
            add_tag_from_tags!(labels, name, tags)
            push!(dirichlet_names, name)
        end



        V = TestFESpace(model, reffe, labels=labels, conformity=:H1, dirichlet_tags=dirichlet_names)
        U = TrialFESpace(V, [dirichlet[tag] for tag in dirichlet_tags])
    else
        V = TestFESpace(model, reffe; conformity=:H1, constraint=:zeromean)
        U = TrialFESpace(V)
    end



    # Integration 
    degree = 2
    Ω = Triangulation(model) # Integration mesh of the domain Omega
    dΩ = Measure(Ω, degree) # Gauss-like quadrature (how to connect point in numerical integration)


    # neumann_tags = neumann_conditions ? collect(keys(neumann)) : Vector{Int}()
    if neumann_conditions
        neumann_tags = collect(keys(neumann))


        Γ = [BoundaryTriangulation(model, tags=tag) for tag in neumann_tags]
        dΓ = [Measure(Γ[i], degree) for i in 1:length(neumann_tags)]

        h = [neumann[tag] for tag in neumann_tags]



    end



    # --- Weak form --- #
    a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dΩ
    b(v) = neumann_conditions ? ∫(v ⋅ f) * dΩ + sum([∫(v ⋅ h[i]) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v ⋅ f) * dΩ
    # b(v) = neumann_conditions ? ∫(v * f) * dΩ + ∫(v * h) * dΓ : ∫(v * f) * dΩ

    # b(v) = ∫(v * f) * dΩ + ∫(v * h) * dΓ




    # --- Solve --- #
    op = AffineFEOperator(a, b, U, V)
    ls = LUSolver()
    solver = LinearFESolver(ls)
    uh = solve(solver, op)


    # --- Compare to manufactured solution --- #
    if MMS == nothing
        println("No manufactured solution (MMS)")
        writevtk(Ω, path * "poisson_results", cellfields=["uh" => uh])
    else
        println("Using (MMS)")

        error = MMS - uh
        l2norm = sqrt(sum(∫(error ⋅ error) * dΩ))
        @printf("l2 norm = %e \n", l2norm)

        writevtk(Ω, path * "poisson_results", cellfields=["uh" => uh, "error" => error])

    end


end




model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/dummy.msh")

# u(x) = [sin(pi * x[1]), sin(pi * x[2])]
# f(x) = sin(pi * x[1]) * pi^2 + sin(pi * x[2]) * pi^2


u(x) = VectorValue(x[1]^2, x[2])
f(x) = VectorValue(-2, 0)

h1(x) = VectorValue(x[1] + x[2], 0)
h2(x) = 100 * x[2]

# dirichlet = Dict(1 => u, 2 => u, 3 => u, 4 => u)
# dirichlet = Dict(1 => -1, 2 => 1, 3 => -1, 4 => 1)


dirichlet = Dict([1, 2, 3, 4] => u)


# neumann = Dict([1,2] => h1, 3 =>  h2)
neumann = Dict([4] => h1)



# dirichlet = Dict()
# neumann = Dict()
poisson(model, f, dirichlet, neumann, u)
