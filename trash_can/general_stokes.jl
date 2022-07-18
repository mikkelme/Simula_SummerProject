using Gridap
using GridapGmsh
using Printf
using Plots

include("./unit_box_mesh.jl")


function stokes(model, f, dirichlet, neumann, MS=nothing, write=false)
    path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/"
    if !ispath(path)
        path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/PDE_stuff/"
    end
    @show path

    dirichlet_conditions = !isempty(dirichlet)
    neumann_conditions = !isempty(neumann)


    # Define reference FE (Q2/P1(disc) pair)
    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order - 1, space=:P)

    labels = get_face_labeling(model)

    dirichlet_tags = collect(keys(dirichlet))
    dirichlet_names = Vector{String}()
    for tags in dirichlet_tags
        name = "$tags"
        add_tag_from_tags!(labels, name, tags)
        push!(dirichlet_names, name)
    end
    @show dirichlet_names

    if !neumann_conditions

        # Define test FESpaces
        V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=dirichlet_names)
        Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
        Y = MultiFieldFESpace([V, Q])

        # Define trial FESpaces from Dirichlet values
        U = TrialFESpace(V, [dirichlet[tag] for tag in dirichlet_tags])
        P = TrialFESpace(Q)
        X = MultiFieldFESpace([U, P])

    else
        # Define test FESpaces
        V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=dirichlet_names)
        Q = TestFESpace(model, reffeₚ, conformity=:L2)

        Y = MultiFieldFESpace([V, Q])

        # Define trial FESpaces from Dirichlet values
        U = TrialFESpace(V, [dirichlet[tag] for tag in dirichlet_tags])
        P = TrialFESpace(Q)
        X = MultiFieldFESpace([U, P])

    end

    # Define triangulation and integration measure
    degree = order
    Ωₕ = Triangulation(model)
    dΩ = Measure(Ωₕ, degree)


    if neumann_conditions
        neumann_tags = collect(keys(neumann))
        @show neumann_tags
        Γ = [BoundaryTriangulation(model, tags=tag) for tag in neumann_tags]
        dΓ = [Measure(Γ[i], degree) for i in 1:length(neumann_tags)]
        h = [neumann[tag] for tag in neumann_tags]
    end

    # New equations
    μ = 1.0
    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))
    a((u, p), (v, q)) = ∫(2 * μ * ϵ(u) ⊙ ϵ(v) - p * (∇ ⋅ v) - q * (∇ ⋅ u)) * dΩ
    b((v, q)) = neumann_conditions ? ∫(v ⋅ f) * dΩ + sum([∫(v ⋅ h[i]) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v ⋅ f) * dΩ


    # println(h[1]([0.0, 0.0]))
    # return
    # # For babuska
    # a((u, λ), (v, μ)) = ∫( ∇u⋅∇v + λ*v + μ*u) * dΩ
    # b((v, q)) = ∫(v ⋅ f) * dΩ - ∫(g * λ) * dΓ 


    # # Old equations
    # a((u, p), (v, q)) = ∫(∇(v) ⊙ ∇(u) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ
    # b((v, q)) = ∫(v ⋅ f)dΩ

    # Build affine FE operator
    op = AffineFEOperator(a, b, X, Y)

    # Solve
    uh, ph = solve(op)

    @show norm(collect(op.op.vector), 2)

    dΩ = Measure(Ωₕ, degree + 2)
    if MS == nothing
        println("No manufactured solution (MS)")
        if write
            writevtk(Ωₕ, path * "general_results", order=2, cellfields=["uh" => uh, "ph" => ph])
        end
    else
        println("Using (MS)")
        u_error = MS[1] - uh
        p_error = MS[2] - ph

        u_l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
        p_l2norm = sqrt(sum(∫(p_error ⋅ p_error) * dΩ))

        @printf("u: l2 norm = %e \n", u_l2norm)
        @printf("p: l2 norm = %e \n", p_l2norm)


        if write
            writevtk(Ωₕ, path * "general_results", order=2, cellfields=["uh" => uh, "ph" => ph, "u_error" => u_error, "p_error" => p_error])
        end
        return u_l2norm, p_l2norm
    end


end


function error_conv(f, dirichlet, neumann, MS)


    lc_start = 2
    num_points = 5

    norm = zeros(Float64, num_points, 2)
    lc = zeros(Float64, num_points)
    for p in 1:num_points
        lc[p] = lc_start * (1 / 2)^(p - 1)
        println(p, " ", lc[p])
        create_unit_box(lc[p])

        path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/unit_box.msh"
        if !ispath(path)
            path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/PDE_stuff/unit_box.msh"
        end
        @show path

        model = GmshDiscreteModel(path)
        if p < num_points
            norm[p, :] .= stokes(model, f, dirichlet, neumann, MS)
        else
            norm[p, :] .= stokes(model, f, dirichlet, neumann, MS, true)
        end

        # Put the actual mesh size
        Ω = Triangulation(model)
        lc[p] = minimum(sqrt.(collect(get_cell_measure(Ω))))
    end


    X = ones(num_points, 2)
    X[:, 2] = log.(lc)
    p = plot(lc, norm[:, 1], xaxis=:log, yaxis=:log, label="Velocity", marker=:x)
    p = plot!(lc, norm[:, 2], xaxis=:log, yaxis=:log, label="Pressure", marker=:o)

    println("---------")
    println("mesh size", lc)
    println("u l2 norm: ", norm[:, 1])
    println("p l2 norm: ", norm[:, 2])

    for i in 1:2
        y = log.(norm[:, i])
        b = inv(transpose(X) * X) * transpose(X) * y
        slope = b[2]
        println("Convergence rate: ", slope)
    end
    display(p)
    return


end






# model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/unit_box.msh")


# MS 1
u0(x) = VectorValue(sin(π * x[2]), cos(π * x[1]))
p0(x) = sin(π * (x[1] + x[2]))

f(x) = VectorValue(π^2 * sin(π * x[2]) + π * cos(π * (x[1] + x[2])), π^2 * cos(π * x[1]) + π * cos(π * (x[1] + x[2])))
σ(x) = [-p0(x) π*(cos(π * x[2])-sin(π * x[1])); π*(cos(π * x[2])-sin(π * x[1])) -p0(x)]



# MS 2
# u0(x) = VectorValue(x[2]^2, x[1]^2)
# p0(x) = (x[1] - 1 / 2)^3 + (x[2] - 1 / 2)^3

# σ(x) = [-p0(x) 2*(x[1]+x[2]); 2*(x[1]+x[2]) -p0(x)]
# f(x) = VectorValue(3 * (x[1] - 1 / 2)^2 - 2, 3 * (x[2] - 1 / 2)^2 - 2)




# (8)  3  (7)
#  4       2
# (5)  1  (6)

# NOTE: 5:8 are marked corner vertices; with only 1:4 the corners would be
# missing
# NOTE: isempty(dirichlet) is not expected to work; some dirichlet boundary
# should always be preset. 

h1(x) = VectorValue(σ(x) * [0.0, -1.0])
h2(x) = VectorValue(σ(x) * [1.0, 0.0])
h3(x) = VectorValue(σ(x) * [0.0, 1.0])
h4(x) = VectorValue(σ(x) * [-1.0, 0.0])

# This one works
# dirichlet = Dict([1, 2, 3, 4, 5, 6, 7, 8] => u0)
# neumann = Dict()



# # FIXME
# h1(x) = VectorValue([1000.0, 1000.0])
dirichlet = Dict([2, 3, 4, 5, 6, 7, 8] => u0)
neumann = Dict(1 => h1)



# path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/unit_box.msh"
# lc = 1
# create_unit_box(lc)
# model = GmshDiscreteModel(path)

# Ω = Triangulation(model)
# lc = minimum(sqrt.(collect(get_cell_measure(Ω))))

# stokes(model, f, dirichlet, neumann, (u0, p0), true)
error_conv(f, dirichlet, neumann, (u0, p0))

# x = [0.0, 0.0]
# println(h1(x))
# println([-p0(x), 2(x[1] + x[2])])



