using Gridap
using GridapGmsh
using Printf
# include("./unit_box_mesh.jl")


function stokes(model, f, dirichlet, neumann, MMS=nothing)
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
        U = TrialFESpace(V, [dirichlet[tag] for tag in dirichlet_tags])
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


    if neumann_conditions
        neumann_tags = collect(keys(neumann))
        Γ = [BoundaryTriangulation(model, tags=tag) for tag in neumann_tags]
        dΓ = [Measure(Γ[i], degree) for i in 1:length(neumann_tags)]
        h = [neumann[tag] for tag in neumann_tags]
        
    end
    
    
    
    # New equations
    μ = 1
    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))
    a((u, p), (v, q)) = ∫(2 * μ * ϵ(u) ⊙ ϵ(v) - p * (∇ ⋅ v) - q * (∇ ⋅ u)) * dΩ
    b((v, q)) = neumann_conditions ? ∫(v ⋅ f) * dΩ + sum([∫(v ⋅ h[i]) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v ⋅ f) * dΩ
    # b((v, q)) = ∫(v ⋅ f) * dΩ #+   ∫(v ⋅ h) * dΓ
    
    
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
        l2norm = sqrt(sum(∫(error ⋅ error) * dΩ))
        @printf("l2 norm = %e \n", l2norm)

        writevtk(Ωₕ, path * "general_results", order=2, cellfields=["uh" => uh, "ph" => ph, "error" => error])

    end


end


# function error_conv()
#     lc = 0.1
#     num_points = 3
#     for i in 1:num_points
#         create_unit_box(lc)
#         model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/unit_box.msh")
#         stokes(model, f, dirichlet, neumann, MMS=nothing)


# end



model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/unit_box.msh")


# u(x) = VectorValue(x[1]^2, -2 * x[1] * x[2])
# f(x) = VectorValue(2.0, 0.0)

# h1(x) = VectorValue(x[1] + x[2], 0)

# u(x) = VectorValue(cos(pi*x[1]), pi*sin(pi*x[1])*x[2])
# f(x) = VectorValue(-pi^2*cos(pi*x[1]), -pi^3*sin(pi*x[1]*x[2]))




u0(x) = VectorValue(sin(π * x[2]), cos(π * x[1]))
p0(x) = sin(π * (x[1] + x[2]))
# μ = 1
f(x) = VectorValue(-2 * π * cos(π * (x[1] + x[2])), π^2 * (cos(π * x[1]) + sin(π * x[2])))


σ(x) = [sin(π * (x[1] + x[2])) π*(cos(π * x[1])+sin(π * x[2])); π*(cos(π * x[1])+sin(π * x[2])) sin(π * (x[1] + x[2]))]


h1 = σ([0.0, -1.0]) * [0.0, -1.0]
h2 = σ([1.0, 0.0]) * [1.0, 0.0]
h3 = σ([0.0, 1.0]) * [0.0, 1.0]
h4 = σ([-1.0, 0.0]) * [-1.0, 0.0]


dirichlet = Dict([1, 2, 3, 4] => u0)
# neumann = Dict(1 => h1, 2 => h2, 3 => h3, 4 => h4)
neumann = Dict(1 => h1)

# neumann = Dict()
stokes(model, f, dirichlet, neumann, u0)



