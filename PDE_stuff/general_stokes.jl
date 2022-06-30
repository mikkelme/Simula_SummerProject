using Gridap
using GridapGmsh



function stokes(model, f, dirichlet, MMS=nothing)
    path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/"
    dirichlet_conditions = !isempty(dirichlet)
    neumann_conditions = !isempty(neumann)


    # Define reference FE (Q2/P1(disc) pair)
    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order - 1; space=:P)


    if dirichlet_conditions

        labels = get_face_labeling(model)

        dirichlet_tags = collect(keys(dirichlet))
        dirichlet_names = Vector{String}()
        for tags in dirichlet_tags
            name = "$tags"
            add_tag_from_tags!(labels, name, tags)
            push!(dirichlet_names, name)
        end


        # Define test FESpaces
        V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=dirichlet_names)
        Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
        Y = MultiFieldFESpace([V, Q])

        # Define trial FESpaces from Dirichlet values

        # println([dirichlet[tag] for tag in dirichlet_tags])
        # U = TrialFESpace(V, [dirichlet[tag] for tag in dirichlet_tags])
        u0 = VectorValue(0, 0)
        u1 = VectorValue(1, 0)
        U = TrialFESpace(V, [u0, u1])

        P = TrialFESpace(Q)
        X = MultiFieldFESpace([U, P])

    else
        println("Handle lack of dirichlet conditions")
        exit()
    end


    # Define triangulation and integration measure
    degree = order
    Ωₕ = Triangulation(model)
    dΩ = Measure(Ωₕ, degree)

    μ = 1

    # New equations
    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))
    a((u, p), (v, q)) = ∫(2 * μ * ϵ(u) ⊙ ϵ(v) - p * (∇ ⋅ v) - q * (∇ ⋅ u)) * dΩ
    b((v, q)) = ∫(v ⋅ f) * dΩ #+   ∫(v ⋅ h) * dΓ



    # # Old equations
    # a((u, p), (v, q)) = ∫(∇(v) ⊙ ∇(u) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ
    # b((v, q)) = ∫(v ⋅ f)dΩ

    # Build affine FE operator
    op = AffineFEOperator(a, b, X, Y)

    # Solve
    uh, ph = solve(op)


    if MMS == nothing
        println("No manufactured solution (MMS)")
        writevtk(Ωₕ, path * "general_results", order=2, cellfields=["uh" => uh, "ph" => ph])
    else
        println("Using (MMS)")

        error = MMS - uh
        println()
        writevtk(Ωₕ, path * "general_results", order=2, cellfields=["uh" => uh, "ph" => ph, "error" => error])

    end


    # # Export results to vtk
    # writevtk(Ωₕ, path * "general_results", order=2, cellfields=["uh" => uh, "ph" => ph])
end



model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/dummy.msh")


u(x) = VectorValue(x[1]^2, -2 * x[1] * x[2])
f = VectorValue(2.0, 0.0)

dirichlet = Dict([1,2,3,4] => u)
# dirichlet = Dict([1,2,3,4] => VectorValue(0.0, 0.0))


# dirichlet = Dict([1] => u)


stokes(model, f, dirichlet, u)



# Problems to investigate friday
#---> Only seems to accept to groups of tags
#---> Check manufactured solution
#---> Include possibility for neumann conditions
