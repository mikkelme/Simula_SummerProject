using Gridap
using GridapGmsh
using Printf
using Plots

include("./unit_box_mesh.jl")


function stokes(model, f, dirichlet, neumann, MMS=nothing, write=false)
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
        if write
            writevtk(Ωₕ, path * "general_results", order=2, cellfields=["uh" => uh, "ph" => ph])
        end
    else
        println("Using (MMS)")
        u_error = MMS[1] - uh
        p_error = MMS[2] - ph
        u_l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
        p_l2norm = sqrt(sum(∫(p_error ⋅ p_error) * dΩ))

        @printf("u: l2 norm = %e \n", u_l2norm)
        @printf("p: l2 norm = %e \n", p_l2norm)


        if write
            writevtk(Ωₕ, path * "general_results", order=2, cellfields=["uh" => uh, "ph" => ph, "u_error" => u_error, "p_error" => u_error])
        end
        return u_l2norm, p_l2norm
    end


end


function error_conv(f, dirichlet, neumann, MMS)


    lc_start = 2
    num_points = 4

    norm = zeros(Float64, num_points, 2)
    lc = zeros(Float64, num_points)
    for p in 1:num_points
        lc[p] = lc_start * (1 / 2)^(p - 1)
        println(p, " ", lc[p])
        create_unit_box(lc[p])
        model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/unit_box.msh")
        if p < num_points
            norm[p, :] .= stokes(model, f, dirichlet, neumann, MMS)
        else
            norm[p, :] .= stokes(model, f, dirichlet, neumann, MMS, true)
        end
    end


    X = ones(num_points, 2)
    X[:, 2] = log.(lc)
    y = log.(norm[:,1])
    p = plot(lc, norm[:, 1], xaxis=:log, yaxis=:log)
    println("---------")
    println("u l2 norm: ", norm[:, 1])
    println("p l2 norm: ", norm[:, 2])
    
    b = inv(transpose(X) * X) * transpose(X)*y
    println("Linear fit: ", b)
    display(p)
    return



end





# u(x) = VectorValue(x[1]^2, -2 * x[1] * x[2])
# f(x) = VectorValue(2.0, 0.0)

# h1(x) = VectorValue(x[1] + x[2], 0)

# u(x) = VectorValue(cos(pi*x[1]), pi*sin(pi*x[1])*x[2])
# f(x) = VectorValue(-pi^2*cos(pi*x[1]), -pi^3*sin(pi*x[1]*x[2]))



# model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/unit_box.msh")

u0(x) = VectorValue(sin(π * x[2]), cos(π * x[1]))
p0(x) = sin(π * (x[1] + x[2]))

# f(x) = VectorValue(-2 * π * cos(π * (x[1] + x[2])), π^2 * (cos(π * x[1]) + sin(π * x[2])))
f(x) = VectorValue(π^2*sin(π*x[2]) - π*cos(π*(x[1]+x[2])), π^2*cos(π*x[1]) - π*cos(π*(x[1]+x[2])))



σ(x) = [sin(π * (x[1] + x[2])) π*(cos(π * x[1])+sin(π * x[2])); π*(cos(π * x[1])+sin(π * x[2])) sin(π * (x[1] + x[2]))]


h1 = VectorValue(σ([0.0, -1.0]) * [0.0, -1.0])
h2 = VectorValue(σ([1.0, 0.0]) * [1.0, 0.0])
h3 = VectorValue(σ([0.0, 1.0]) * [0.0, 1.0])
h4 = VectorValue(σ([-1.0, 0.0]) * [-1.0, 0.0])


dirichlet = Dict([1, 2, 3, 4] => u0)
# neumann = Dict(1 => h1, 2 => h2, 3 => h3, 4 => h4)
neumann = Dict()

# neumann = Dict()
# stokes(model, f, dirichlet, neumann, (u0, p0))
error_conv(f, dirichlet, neumann, (u0, p0))




