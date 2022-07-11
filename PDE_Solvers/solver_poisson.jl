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


function poisson_solver(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags; write=false)
    neumann_conditions = !isempty(neumann_tags)
    
    labels = get_face_labeling(model)
    dirichlet_tags = [pgs_dict[tag] for tag in dirichlet_tags]
    neumann_tags = [pgs_dict[tag] for tag in neumann_tags]
    add_tag_from_tags!(labels, "diri", dirichlet_tags)


    order = 2
    reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)


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
    end


    # --- Weak form --- #
    # Scalar
    # a(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ
    # b(v) = neumann_conditions ? ∫(v * f0) * dΩ + sum([∫(v * (h0 ⋅ ν[i])) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v * f0) * dΩ

    # Vector
    a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dΩ
    b(v) = neumann_conditions ? ∫(v ⋅ f0) * dΩ + sum([∫(v ⋅ (h0 ⋅ ν[i])) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v ⋅ f0) * dΩ


    # --- Solve --- #
    op = AffineFEOperator(a, b, U, V)
    ls = LUSolver()
    solver = LinearFESolver(ls)
    uh = solve(solver, op)

    # --- Check with manufactured solution --- #
    error = g0 - uh
    l2norm = sqrt(sum(∫(error ⋅ error) * dΩ))
    @printf("l2 norm = %e \n", l2norm)

    if write
        writevtk(Ω, path * "poisson_results", cellfields=["uh" => uh, "error" => error])
    end

    return l2norm


end



function error_conv(solver, f0, g0, h0, dirichlet_tags, neumann_tags; lc_start=2, num_points=5, show_plot=false)
    l2norm = zeros(Float64, num_points)
    lc = zeros(Float64, num_points)

    # --- Decrease mesh size and get l2norm --- #
    write = false
    for p in 1:num_points
        lc[p] = lc_start * (1 / 2)^(p - 1)
        p == num_points && (write = true)
        model, pgs_dict = create_unit_box(lc[p])
        l2norm[p] = solver(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags; write)

        # Get the actual mesh size
        Ω = Triangulation(model)
        lc[p] = minimum(sqrt.(collect(get_cell_measure(Ω))))
    end

    println("Mesh size: ", lc)
    println("l2 norm: ", l2norm)


    # ---  Convergence rate --- #
    X = ones(num_points, 2)
    X[:, 2] = log.(lc)

    y = log.(l2norm)
    b = inv(transpose(X) * X) * transpose(X) * y
    slope = b[2]
    println("Convergence rate: ", slope)

    if show_plot
        fig = plot(lc, l2norm, xaxis=:log, yaxis=:log, label="Velocity", marker=:x)
        display(fig)
    end

end




# MS (scalar)
# u0(x) = cos(π * x[1] * x[2])
# f0(x) = π^2 * (x[1]^2 + x[2]^2) * cos(π * x[1] * x[2])
# h0(x) = VectorValue(-π * x[2] * sin(π * x[1] * x[2]), -π * x[1] * sin(π * x[1] * x[2]))

# MS (vector)
u0(x) = VectorValue(sin(π * x[1]), sin(π * x[2]))
f0(x) = VectorValue(π^2 * sin(π * x[1]), π^2 * sin(π * x[2]))
h0(x) = TensorValue(π * cos(π * x[1]), 0.0, 0.0, π * cos(π * x[2]))




# Boundary tags
all_btags = [1, 2, 3, 4, 5, 6, 7, 8]
dirichlet_tags = [1, 2, 3, 4, 5, 6, 7, 8]
neumann_tags = filter(x -> x ∉ dirichlet_tags, all_btags)


# model, pgs_dict = create_unit_box(2, false)
# poisson_solver(model, pgs_dict, f0, u0, h0, dirichlet_tags, neumann_tags; write = true)
error_conv(poisson_solver, f0, u0, h0, dirichlet_tags, neumann_tags; show_plot=false)



