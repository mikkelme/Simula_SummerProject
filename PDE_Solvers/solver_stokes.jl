using Gridap
using GridapGmsh
using Printf
using Plots


include("../mesh_generators/unit_box_direct.jl")

path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_Solvers/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/PDE_Solvers/"
end
@show path




function stokes_solver(model, pgs_dict, f0, g0, h0, dirichlet_tags, neumann_tags; write=false)
    labels = get_face_labeling(model)
    neumann_conditions = !isempty(neumann_tags)

    dirichlet_tags = [pgs_dict[tag] for tag in dirichlet_tags]
    neumann_tags = [pgs_dict[tag] for tag in neumann_tags]


    add_tag_from_tags!(labels, "diri", dirichlet_tags)


    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order - 1, space=:P)

    V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["diri"])
    !neumann_conditions ? Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean) : Q = TestFESpace(model, reffeₚ, conformity=:L2)
    Y = MultiFieldFESpace([V, Q])

    U = TrialFESpace(V, g0[1])
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U, P])



    # Integration 
    degree = order
    Ωₕ = Triangulation(model)
    dΩ = Measure(Ωₕ, degree)


    if neumann_conditions
        Γ = [BoundaryTriangulation(model, tags=tag) for tag in neumann_tags]
        ν = [get_normal_vector(Γ[i]) for i in 1:length(neumann_tags)]
        dΓ = [Measure(Γ[i], degree) for i in 1:length(neumann_tags)]
    
    end



    μ = 1.0
    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))
    a((u, p), (v, q)) = ∫(2 * μ * ϵ(u) ⊙ ϵ(v) - p * (∇ ⋅ v) - q * (∇ ⋅ u)) * dΩ
    b((v, q)) = neumann_conditions ? ∫(v ⋅ f0) * dΩ + sum([∫(v ⋅ (h0 ⋅ ν[i])) * dΓ[i] for i in 1:length(neumann_tags)]) : ∫(v ⋅ f0) * dΩ



    # --- Solve --- #
    op = AffineFEOperator(a, b, X, Y)
    uh, ph = solve(op)
    # ls = LUSolver()
    # solver = LinearFESolver(ls)
    # uh = solve(solver, op)


    # --- Check with manufactured solution --- #
    u_error = g0[1] - uh
    p_error = g0[2] - ph

    u_l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
    p_l2norm = sqrt(sum(∫(p_error ⋅ p_error) * dΩ))

    @printf("u: l2 norm = %e \n", u_l2norm)
    @printf("p: l2 norm = %e \n", p_l2norm)

    write && writevtk(Ωₕ, path * "vtu_files/" * "stokes_results", order=2, cellfields=["uh" => uh, "ph" => ph, "u_error" => u_error, "p_error" => p_error])
  

    return u_l2norm, p_l2norm


end


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
    println("p: l2 norm =  ", l2norm[:, 2])


    # ---  Convergence rate --- #
    X = ones(num_points, 2)
    X[:, 2] = log.(lc)
    for i in 1:2
        y = log.(l2norm[:, i])
        b = inv(transpose(X) * X) * transpose(X) * y
        slope[i] = b[2]
      
    end

    println("u: Convergence rate: ", slope[1])
    println("p: Convergence rate: ", slope[2])

    if show_plot
        fig = plot(lc, l2norm[:, 1], xaxis=:log, yaxis=:log, label="Velocity", marker=:x)
        fig = plot!(lc, l2norm[:, 2], xaxis=:log, yaxis=:log, label="Pressure", marker=:o)
        display(fig)
    end

end






# MS (scalar)
# u0(x) = cos(π * x[1] * x[2])
# f0(x) = π^2 * (x[1]^2 + x[2]^2) * cos(π * x[1] * x[2])
# h0(x) = VectorValue(-π * x[2] * sin(π * x[1] * x[2]), -π * x[1] * sin(π * x[1] * x[2]))

# MS (vector)
u0(x) = VectorValue(sin(π * x[2]), cos(π * x[1]))
p0(x) = sin(π * (x[1] + x[2]))
f0(x) = VectorValue(π^2 * sin(π * x[2]) + π * cos(π * (x[1] + x[2])), π^2 * cos(π * x[1]) + π * cos(π * (x[1] + x[2])))
σ(x) = TensorValue(-p0(x), π * (cos(π * x[2]) - sin(π * x[1])), π * (cos(π * x[2]) - sin(π * x[1])), -p0(x))



# Boundary tags
all_btags = [1, 2, 3, 4, 5, 6, 7, 8]
dirichlet_tags = [3, 4, 5, 7, 8]
neumann_tags = filter(x -> x ∉ dirichlet_tags, all_btags)


# model, pgs_dict = create_unit_box(2, false)
# stokes_solver(model, pgs_dict, f0, (u0, p0), σ, dirichlet_tags, neumann_tags; write=true)
error_conv(stokes_solver, f0, (u0, p0), σ, dirichlet_tags, neumann_tags; show_plot = false)






