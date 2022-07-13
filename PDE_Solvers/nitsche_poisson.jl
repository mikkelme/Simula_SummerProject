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
- Δu = f0 on Ω
u = g0 on ΓD
∂uₙ = h0 on ΓN    
"""

function nitsche_poisson(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags; write=false)

    dirichlet_tags = [pgs_dict[tag] for tag in dirichlet_tags]
    neumann_tags = [pgs_dict[tag] for tag in neumann_tags]

    Ω = Triangulation(model)
    ΓD = BoundaryTriangulation(model, tags=dirichlet_tags)
    ΓN = BoundaryTriangulation(model, tags=neumann_tags)

    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)


    δV = TestFESpace(Ω, reffeᵤ, conformity=:H1)
    # V = δV # Trial space same as test space
    V = TrialFESpace(δV) # Trial space same as test space

    degree = 2 * order
    dΩ = Measure(Ω, degree)
    dΓD = Measure(ΓD, degree)
    dΓN = Measure(ΓN, degree)

    # Nitsche 
    n̂D = get_normal_vector(ΓD)
    n̂N = get_normal_vector(ΓN)
    γ = 10^order # Nitsche penalty parameter

    # And finally we need the notian of cell diameter
    h_ΓD = get_array(∫(1) * dΓD)
    h = CellField(lazy_map(h -> h, h_ΓD), ΓD)


    # # Scalar 
    # aN(u, v) = ∫(γ / h * u * v) * dΓD - ∫((∇(u) ⋅ n̂D) * v) * dΓD - ∫((∇(v) ⋅ n̂D) * u) * dΓD
    # bN(v) = ∫(γ / h * g0 * v) * dΓD - ∫((∇(v) ⋅ n̂D) * g0) * dΓD
    # a(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ + aN(u, v)
    # b(v) = ∫(v * f0) * dΩ + bN(v) + ∫(v * (h0 ⋅ n̂N)) * dΓN

    # Vector 
    aN(u, v) = ∫(γ / h * u ⋅ v) * dΓD - ∫((∇(u) ⋅ n̂D) ⋅ v) * dΓD - ∫((∇(v) ⋅ n̂D) ⋅ u) * dΓD
    bN(v) = ∫(γ / h * g0 ⋅ v) * dΓD - ∫((∇(v) ⋅ n̂D) ⋅ g0) * dΓD

    a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dΩ + aN(u, v)
    b(v) = ∫(v ⋅ f0) * dΩ + bN(v) + ∫(v ⋅ (h0 ⋅ n̂N)) * dΓN

    op = AffineFEOperator(a, b, V, δV)
    uh = solve(op)

    # --- Check with manufactured solution --- #
    u_error = g0 - uh
    l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
    @printf("l2 norm = %e \n", l2norm)

    write && writevtk(Ω, path * "vtu_files/" * "nitsche_poisson_results", cellfields=["uh" => uh, "u_error" => u_error])

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



# # MS (scalar)
# u0(x) = cos(π * x[1] * x[2])
# f0(x) = π^2 * (x[1]^2 + x[2]^2) * cos(π * x[1] * x[2])
# h0(x) = VectorValue(-π * x[2] * sin(π * x[1] * x[2]), -π * x[1] * sin(π * x[1] * x[2]))

# MS (vector)
u0(x) = VectorValue(sin(π * x[1]), sin(π * x[2]))
f0(x) = VectorValue(π^2 * sin(π * x[1]), π^2 * sin(π * x[2]))
h0(x) = TensorValue(π * cos(π * x[1]), 0.0, 0.0, π * cos(π * x[2]))

# Boundary tags
all_btags = [1, 2, 3, 4, 5, 6, 7, 8]
dirichlet_tags = [1, 2, 3, 5, 6, 7, 8]
neumann_tags = filter(x -> x ∉ dirichlet_tags, all_btags)


# model, pgs_dict = create_unit_box(2, false)
# nitsche_poisson(model, pgs_dict, f0, u0, h0, dirichlet_tags, neumann_tags; write = true)
error_conv(nitsche_poisson, f0, u0, h0, dirichlet_tags, neumann_tags; show_plot=false)

