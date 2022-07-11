using Gridap
using GridapGmsh
using Printf
using Plots

include("./unit_box_direct.jl")
path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_Solvers/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/PDE_Solvers/"
end
@show path


"""
- Δu + u = f on Ω
u = g on ΓD
∂uₙ = h on ΓN
u = g on Γ    
"""


function poisson_solver(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags; write=false)
    labels = get_face_labeling(model)
    neumann_conditions = !isempty(neumann_tags)

    dirichlet_tags = [pgs_dict[tag] for tag in dirichlet_tags]
    neumann_tags = [pgs_dict[tag] for tag in neumann_tags]
    add_tag_from_tags!(labels, "diri", dirichlet_tags) # is this used????


    Ω = Triangulation(model)
    ΓD = BoundaryTriangulation(model, tags=dirichlet_tags)
    ΓN = BoundaryTriangulation(model, tags=neumann_tags)

    order = 2
    Velm = ReferenceFE(lagrangian, Float64, order)
    Melm = ReferenceFE(lagrangian, Float64, order)

    # What is going on here with Ω instead of "model"
    V = TestFESpace(Ω, Velm, conformity=:H1)
    M = TestFESpace(ΓD, Velm, conformity=:H1)
    W = MultiFieldFESpace([V, M])

    degree = 3 # might introduce some more flexibility here
    dΩ = Measure(Ω, order)
    dΓD = Measure(ΓD, order)
    dΓN = Measure(ΓN, order)
    ν = get_normal_vector(ΓN)



    # Bilinear form
    # Scalar
    a((u, λ), (v, μ)) = ∫(∇(u) ⋅ ∇(v)) * dΩ + ∫(u * v) * dΩ - ∫(λ * v) * dΓD - ∫(u * μ) * dΓD
    b((v, μ)) = neumann_conditions ? ∫(v * f0) * dΩ - ∫(g0 * μ) * dΓD + ∫(v * (h0 ⋅ ν)) * dΓN : ∫(v * f0) * dΩ - ∫(g0 * μ) * dΓD

    # Vector 
    # a((u, λ), (v, μ)) = ∫(∇(u) ⊙ ∇(v)) * dΩ + ∫(u ⋅ v) * dΩ - ∫(λ ⋅ v) * dΓD - ∫(u ⋅ μ) * dΓD
    # b((v, μ)) = neumann_conditions ? ∫(v ⋅ f0) * dΩ - ∫(g0 ⋅ μ) * dΓD + ∫(v ⋅ (h0 ⋅ ν)) * dΓN : ∫(v ⋅ f0) * dΩ - ∫(g0 ⋅ μ) * dΓD





    # --- Solve --- #
    op = AffineFEOperator(a, b, W, W)
    uh, λh = solve(op)

    # Why do I need this <----------#
    Γ = get_triangulation(λh)
    dΓ = Measure(Γ, order)
    ν = get_normal_vector(Γ)

    # --- Check with manufactured solution --- #
    u_error = g0 - uh
    λ_error = λh - (h0 ⋅ ν)  # On surface we should be approximating ∇u⋅ν

    u_l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
    λ_l2norm = sqrt(sum(∫(λ_error ⋅ λ_error) * dΓ))

    @printf("u: l2 norm = %e \n", u_l2norm)
    @printf("λ: l2 norm = %e \n", λ_l2norm)

    write && writevtk(Ω, path * "poisson_babuska_results", cellfields=["uh" => uh, "u_error" => u_error])

    return [u_l2norm, λ_l2norm]
end


# function error_conv(solver, f0, g0, h0, dirichlet_tags, neumann_tags; lc_start=2, num_points=5, show_plot=false)
#     l2norm = zeros(Float64, num_points)
#     lc = zeros(Float64, num_points)

#     # --- Decrease mesh size and get l2norm --- #
#     write = false
#     for p in 1:num_points
#         lc[p] = lc_start * (1 / 2)^(p - 1)
#         p == num_points && (write = true)
#         model, pgs_dict = create_unit_box(lc[p])
#         l2norm[p] = solver(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags; write)

#         # Get the actual mesh size
#         Ω = Triangulation(model)
#         lc[p] = minimum(sqrt.(collect(get_cell_measure(Ω))))
#     end

#     println("Mesh size: ", lc)
#     println("l2 norm: ", l2norm)


#     # ---  Convergence rate --- #
#     X = ones(num_points, 2)
#     X[:, 2] = log.(lc)

#     y = log.(l2norm)
#     b = inv(transpose(X) * X) * transpose(X) * y
#     slope = b[2]
#     println("Convergence rate: ", slope)

#     if show_plot
#         fig = plot(lc, l2norm, xaxis=:log, yaxis=:log, label="Velocity", marker=:x)
#         display(fig)
#     end

# end



function error_conv(solver, f0, g0, h0, dirichlet_tags, neumann_tags; lc_start=2, num_points=5, show_plot=false)
    l2norm = zeros(Float64, num_points, 2)
    lc = zeros(Float64, num_points)
    slope = zeros(Float64, 2)

    # --- Decrease mesh size and get l2norm --- #
    write = false
    for p in 1:num_points
        lc[p] = lc_start * (1 / 2)^(p - 1)
        p == num_points && (write = true)
        model, pgs_dict = create_unit_box(lc[p])
        l2norm[p, :] .= solver(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags; write)

        # Get the actual mesh size
        Ω = Triangulation(model)
        lc[p] = minimum(sqrt.(collect(get_cell_measure(Ω))))
    end

    println("mesh size: ", lc)
    println("u: l2 norm = ", l2norm[:, 1])
    println("λ: l2 norm =  ", l2norm[:, 2])


    # ---  Convergence rate --- #
    X = ones(num_points, 2)
    X[:, 2] = log.(lc)
    for i in 1:2
        y = log.(l2norm[:, i])
        b = inv(transpose(X) * X) * transpose(X) * y
        slope[i] = b[2]

    end

    println("u: Convergence rate: ", slope[1])
    println("λ: Convergence rate: ", slope[2])

    if show_plot
        fig = plot(lc, l2norm[:, 1], xaxis=:log, yaxis=:log, label="Velocity", marker=:x)
        fig = plot!(lc, l2norm[:, 2], xaxis=:log, yaxis=:log, label="Multiplier", marker=:o)
        display(fig)
    end
end


# MS (scalar)
u0(x) = cos(π * x[1] * x[2])
f0(x) = π^2 * (x[1]^2 + x[2]^2 + 1 / π^2) * cos(π * x[1] * x[2])
h0(x) = VectorValue(-π * x[2] * sin(π * x[1] * x[2]), -π * x[1] * sin(π * x[1] * x[2])) # <---- Unsure of this



all_btags = [1, 2, 3, 4, 5, 6, 7, 8]
dirichlet_tags = [1, 4, 5, 6, 7, 8]
neumann_tags = filter(x -> x ∉ dirichlet_tags, all_btags)


# model, pgs_dict = create_unit_box(2, false)
# poisson_solver(model, pgs_dict, f0, u0, h0, dirichlet_tags, neumann_tags; write=true)
error_conv(poisson_solver, f0, u0, h0, dirichlet_tags, neumann_tags; show_plot=false)