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



function nitsche_stokes_solver(model, pgs_dict, f0, g0; write=false)

    # --- Boundary tags --- #
    Γ1_tags = pgs_tags(pgs_dict, [1, 3, 5, 6, 7, 8])
    Γ2_tags = pgs_tags(pgs_dict, [2, 4])
    neumann_tags = Vector{Int}()
    u0 = g0[1]; p0 = g0[2]

    # --- Triangulation and spaces --- #
    Ω = Triangulation(model)
    ΓD = BoundaryTriangulation(model, tags=Γ2_tags)
    ΓN = BoundaryTriangulation(model, tags=neumann_tags)

    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order - 1)

    δV = TestFESpace(Ω, reffeᵤ, conformity=:H1, dirichlet_tags = Γ1_tags)
    δQ = TestFESpace(model, reffeₚ, conformity=:H1)

    V = TrialFESpace(δV, u0)
    Q = TrialFESpace(δQ)

    δW = MultiFieldFESpace([δV, δQ])
    W = MultiFieldFESpace([V, Q])

    degree = 2 * order
    dΩ = Measure(Ω, degree)
    dΓD = Measure(ΓD, degree)
    dΓN = Measure(ΓN, degree)

    # Nitsche 
    n̂D = get_normal_vector(ΓD)
    t̂D = TensorValue(0, -1, 1, 0) ⋅ n̂D # Tangent
    γ = 10^order # Nitsche penalty parameter

    h_ΓD = get_array(∫(1) * dΓD)
    h = CellField(lazy_map(h -> h, h_ΓD), ΓD)
    
    # --- Weak formulation --- #
    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))
    σ(u,p) = 2 * μ * ε(u) - p * TensorValue(1, 0, 0, 1)

    gΓ2 = u0 × n̂D #

    aN((u, p), (v, q)) = ∫(γ/h * (u⋅t̂D ) * (v ⋅ t̂D)) * dΓD  - ∫( ((n̂D ⋅ σ(u,p)) ⋅ t̂D) * (v ⋅t̂D) ) * dΓD - ∫( ((n̂D ⋅ σ(v,q)) ⋅ t̂D) * (u⋅t̂D ) ) * dΓD
    bN((v, q)) = ∫(p0 * (-v ⋅ n̂D)) * dΓD + ∫(gΓ2 * (γ/h * (v ⋅ t̂D) - ((n̂D ⋅ σ(v,q)) ⋅ t̂D))) * dΓD

    

    a((u, p), (v, q)) = ∫(2 * μ * ϵ(u) ⊙ ϵ(v) - p * (∇ ⋅ v) - q * (∇ ⋅ u)) * dΩ + aN((u, p), (v, q))
    b((v, q)) = ∫(v ⋅ f0) * dΩ + bN((v, q)) # + ∫(v ⋅ (h0 ⋅ n̂N)) * dΓN 
    
    
    # --- Solve --- #
    op = AffineFEOperator(a, b, W, δW)
    uh, ph = solve(op)

    # --- Check with manufactured solution --- #
    u_error = g0[1] - uh
    p_error = g0[2] - ph
    
    u_l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
    p_l2norm = sqrt(sum(∫(p_error ⋅ p_error) * dΩ))

    @printf("u: l2 norm = %e \n", u_l2norm)
    @printf("p: l2 norm = %e \n", p_l2norm)

    write && writevtk(Ω, path * "vtu_files/" * "nitsche_stokes_results", order=2, cellfields=["uh" => uh, "ph" => ph, "u_error" => u_error, "p_error" => p_error])
    
    return u_l2norm, p_l2norm

end



function error_conv(solver, f0, g0; lc_start=2, num_points=5, show_plot=false)
    l2norm = zeros(Float64, num_points, 2)
    lc = zeros(Float64, num_points)
    slope = zeros(Float64, 2)

    # --- Decrease mesh size and get l2norm --- #
    write = false
    for p in 1:num_points
        lc[p] = lc_start * (1 / 2)^(p - 1)
        p == num_points && (write = true)
        model, pgs_dict = create_unit_box(lc[p])
        l2norm[p, :] .= solver(model, pgs_dict, f0, g0; write)

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
        fig = plot(lc, l2norm[:, 2], xaxis=:log, yaxis=:log, label="Pressure", marker=:x)
        display(fig)
    end
end



μ = 1.0

# # MS 1 (u × n̂ = 0 on Γ2)
# u0(x) = VectorValue(sin(π * x[2]), sin(π * x[1]))
# p0(x) = cos(π * x[1] * x[2])
# f0(x) = VectorValue(μ * π^2 * sin(π * x[2]) - π * x[2] * sin(π*x[1]*x[2]), μ * π^2 * sin(π * x[1]) - π  * x[1] * sin(π*x[1]*x[2]))

# # Visualization of quadratic shape of velocity field (ignore l2 norms here)
# u0(x) = VectorValue(0.0, 0.0) # u = 0 on boundary
# pL = 10.0; pR = 1.0
# p0(x) = (1 - x[1]) * pL + x[1] * pR # Pressure drop from left to right: pL → pR 
# f0(x) = VectorValue(0.0, 0.0)


# MS 2 (u × n̂ ≂̸ 0 on Γ2)
u0(x) = VectorValue(sin(π * x[2]), cos(π * x[1]))
p0(x) = sin(π * (x[1] + x[2]))
f0(x) = VectorValue(μ*π^2 * sin(π * x[2]) + π * cos(π * (x[1] + x[2])), μ*π^2 * cos(π * x[1]) + π * cos(π * (x[1] + x[2])))

g0 = (u0, p0)
# model, pgs_dict = create_unit_box(2, false)
# nitsche_stokes_solver(model, pgs_dict, f0, g0; write=true)
error_conv(nitsche_stokes_solver, f0, g0)



