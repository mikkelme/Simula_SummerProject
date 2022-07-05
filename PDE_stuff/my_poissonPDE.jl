using Gridap
using GridapGmsh
using Printf
using Plots

include("./unit_box_mesh.jl")


# Equation to solve
# Find scalar field u such that 
# -Δu = -∇²u = f, f ∈ Ω, 
# u = g, on boundary
# ∇u⋅n = h, 
# u g 

# model = GmshDiscreteModel("dummy.msh")
# writevtk(model, "model")
function create_square(n)
    # (3) 6 (4)
    #  7     8
    # (1) 5 (2)

    domain = (0, 1, 0, 1)
    partition = (n, n)
    model = CartesianDiscreteModel(domain, partition)
    return model
end


function poisson(model, f, dirichlet, neumann, MMS=nothing, write=false)
    path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/"

    dirichlet_conditions = !isempty(dirichlet)
    neumann_conditions = !isempty(neumann)


    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order) 
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



    # Integration 
    degree = 2
    Ω = Triangulation(model) # Integration mesh of the domain Omega
    dΩ = Measure(Ω, degree) # Gauss-like quadrature (how to connect point in numerical integration)


    # neumann_tags = neumann_conditions ? collect(keys(neumann)) : Vector{Int}()
    if neumann_conditions
        neumann_tags = collect(keys(neumann))
        Γ = [BoundaryTriangulation(model, tags=tag) for tag in neumann_tags]
        ν = [get_normal_vector(Γ[i]) for i in 1:length(neumann_tags)]
        dΓ = [Measure(Γ[i], degree) for i in 1:length(neumann_tags)]
        h = [neumann[tag] for tag in neumann_tags]


    end



    # --- Weak form --- #
    a(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ
    b(v) = neumann_conditions ? ∫(v * f) * dΩ + sum([∫(v * (h[i] ⋅ ν[i])) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v * f) * dΩ
    
    # b(v) = neumann_conditions ? ∫(v * f) * dΩ + sum([∫(v * h[i]) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v * f) * dΩ
    # b(v) = ∫(v * f) * dΩ


    # tag = 3

    # Γ = BoundaryTriangulation(model, tags=[tag])
    # dΓ = Measure(Γ, degree)
    # # h = neumann[tag]
    # h(x) = VectorValue([1000.0, -1000.0])

    # b(v) = ∫(v ⋅ f) * dΩ + ∫(v ⋅ h) * dΓ
    # b(v) = ∫(v ⋅ f) * dΩ 


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
        if write
            writevtk(Ω, path * "poisson_results", cellfields=["uh" => uh])
        end
    else
        println("Using (MMS)")

        error = MMS - uh
        l2norm = sqrt(sum(∫(error ⋅ error) * dΩ))
        @printf("l2 norm = %e \n", l2norm)

        if write
            writevtk(Ω, path * "poisson_results", cellfields=["uh" => uh, "error" => error])
        end

        return l2norm
    end


end




function error_conv(solver, f0, dirichlet, neumann, u0)


    lc_start = 2
    num_points = 5

    norm = zeros(Float64, num_points)
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
        norm[p] = solver(model, f0, dirichlet, neumann, u0)


        # Put the actual mesh size
        Ω = Triangulation(model)
        lc[p] = minimum(sqrt.(collect(get_cell_measure(Ω))))
    end


    X = ones(num_points, 2)
    X[:, 2] = log.(lc)
    p = plot(lc, norm, xaxis=:log, yaxis=:log, label="Velocity", marker=:x)

    println("---------")
    println("mesh size", lc)
    println("u l2 norm: ", norm)

    y = log.(norm)
    b = inv(transpose(X) * X) * transpose(X) * y
    slope = b[2]
    println("Convergence rate: ", slope)

    # display(p)
    return


end




# function error_conv(solver, f0, u0, dirichlet, neumann)
#     path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/"


#     n_start = 4
#     num_points = 5

#     norm = zeros(Float64, num_points)
#     n = zeros(Float64, num_points)
#     for p in 1:num_points

#         n[p] = n_start^p
#         model = create_square(n[p])
#         norm[p] = poisson(model, f0, dirichlet, neumann, u0)

#     end


#     X = ones(num_points, 2)
#     X[:, 2] = log.(n)
#     p = plot(n, norm, xaxis=:log, yaxis=:log, label="Velocity", marker=:x)

#     println("---------")
#     println("u l2 norm: ", norm)

#     y = log.(norm)
#     b = inv(transpose(X) * X) * transpose(X) * y
#     slope = b[2]
#     println("Convergence rate: ", slope)

#     # display(p)
#     return


# end


# model = GmshDiscreteModel("/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_stuff/unit_box.msh")


# u0(x) = VectorValue(sin(π * x[1]), sin(π * x[2]))
# f(x) = VectorValue(π^2 * sin(π * x[1]), π^2 * sin(π * x[2]))

# grad_u(x) = [π*cos(π * x[1]) 0.0; 0.0 π*cos(π * x[2])]

# h1(x) = VectorValue(grad_u(x) * [0.0, -1.0])
# h2(x) = VectorValue(grad_u(x) * [1.0, 0.0])
# h3(x) = VectorValue(grad_u(x) * [0.0, 1.0])
# h4(x) = VectorValue(grad_u(x) * [-1.0, 0.0])



u0(x) = cos(π * x[1] * x[2])
f0(x) = π^2 * (x[1]^2 + x[2]^2) * cos(π * x[1] * x[2])
h0(x) = VectorValue(-π * x[2] * sin(π * x[1] * x[2]), -π * x[1] * sin(π * x[1] * x[2]))



# h1(x) = -π * x[2] * sin(π * x[1] * x[2])


# h1(x) = VectorValue( π*cos( π * x[1] ), 0)

# u(x) = VectorValue(x[1]^2, x[2])
# f(x) = VectorValue(-2, 0)

# h2(x) = 100 * x[2]

# h3(x) = VectorValue([1000.0, -1000.0])
# dirichlet = Dict([1, 2, 3, 4, 5, 6, 7, 8] => u0)

dirichlet = Dict([1, 2, 3, 4, 5, 6, 7] => u0)
neumann = Dict(8 => h0)
# neumann = Dict()



# dirichlet = Dict()
# poisson(model, f, dirichlet, neumann, u)
error_conv(poisson, f0, dirichlet, neumann, u0)




