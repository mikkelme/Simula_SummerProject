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




function babuska_stokes_solver(model, pgs_dict, f0, g0; write=false)

    # --- Boundary conditions --- #
    Γ1_tags = pgs_tags(pgs_dict, [1, 3, 5, 6, 7, 8])
    Γ2_tags = pgs_tags(pgs_dict, [2, 4])
    corner_tags = pgs_tags(pgs_dict, [5, 6, 7, 8])

    # --- Triangulation and spaces --- #
    Ω = Triangulation(model)
    Γ2 = BoundaryTriangulation(model, tags=Γ2_tags)

    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order - 1)
    reffeₘᵤₗ  = ReferenceFE(lagrangian, Float64, order)


    δV = TestFESpace(Ω, reffeᵤ, conformity=:H1, dirichlet_tags=Γ1_tags) # Velocity
    δM = TestFESpace(Ω, reffeₚ, conformity=:H1) # Pressure  
    δΛ = TestFESpace(Γ2, reffeₘᵤₗ, conformity=:H1, dirichlet_tags=corner_tags) # Multiplier

    V = TrialFESpace(δV, g0[1])
    M = TrialFESpace(δM)
    Λ = TrialFESpace(δΛ, 0)

    δW = MultiFieldFESpace([δV, δM, δΛ]) # Test spaces
    W = MultiFieldFESpace([V, M, Λ]) # Trial spaces

    degree = 2 * order
    dΩ = Measure(Ω, degree)
    dΓ2 = Measure(Γ2, degree)
    nΓ2 = get_normal_vector(Γ2)

    # --- Weak formulation --- #
    # μ = 1.0
    # ρ = 1.0
    p0 = g0[2]


    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))
    t = TensorValue(0, -1, 1, 0) ⋅ nΓ2 # Tangent
    c(u, λ) = ∫((u ⋅ t) * λ) * dΓ2


    a((u, p, λ), (v, q, η)) = ∫(2 * μ * ϵ(u) ⊙ ϵ(v) - p * (∇ ⋅ v) - q * (∇ ⋅ u)) * dΩ - c(v, λ) - c(u, η)
    b((v, q, η)) = ∫(ρ * v ⋅ f0) * dΩ - ∫(p0 * nΓ2 ⋅ v) * dΓ2

    # --- Solve --- #
    op = AffineFEOperator(a, b, W, δW)
    uh, ph, λh = solve(op)


    # --- Check with manufactured solution --- #
    u_error = g0[1] - uh
    p_error = g0[2] - ph

    σ(u,p) = 2 * μ * ε(u) - p*TensorValue(1, 0, 0, 1)
    τ(u, p) = ((σ(u,p) ⋅ nΓ2) · t)
    λ_error = λh - τ(uh, ph)

    u_l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
    p_l2norm = sqrt(sum(∫(p_error ⋅ p_error) * dΩ))
    λ_l2norm = sqrt(sum(∫(λ_error ⋅ λ_error) * dΓ2))

    @printf("u: l2 norm = %e \n", u_l2norm)
    @printf("p: l2 norm = %e \n", p_l2norm)
    @printf("λ: l2 norm = %e \n", λ_l2norm)


    write && writevtk(Ω, path * "vtu_files/" * "babuska_stokes_results", order=2, cellfields=["uh" => uh, "ph" => ph, "u_error" => u_error, "p_error" => p_error])

    return u_l2norm, p_l2norm, λ_l2norm

end


function error_conv(solver, f0, g0; lc_start=2, num_points=5, show_plot=false)
    l2norm = zeros(Float64, num_points, 3)
    lc = zeros(Float64, num_points)
    slope = zeros(Float64, 3)

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
    println("λ: l2 norm =  ", l2norm[:, 3])


    # ---  Convergence rate --- #
    X = ones(num_points, 2)
    X[:, 2] = log.(lc)
    for i in 1:3
        y = log.(l2norm[:, i])
        b = inv(transpose(X) * X) * transpose(X) * y
        slope[i] = b[2]

    end

    println("u: Convergence rate: ", slope[1])
    println("p: Convergence rate: ", slope[2])
    println("λ: Convergence rate: ", slope[3])

    if show_plot
        fig = plot(lc, l2norm[:, 1], xaxis=:log, yaxis=:log, label="Velocity", marker=:x)
        fig = plot(lc, l2norm[:, 2], xaxis=:log, yaxis=:log, label="Pressure", marker=:x)
        fig = plot!(lc, l2norm[:, 3], xaxis=:log, yaxis=:log, label="Multiplier", marker=:o)
        display(fig)
    end
end






μ = 1.0
ρ = 1.0

# # MS 1 (Simple)
# u0(x) = VectorValue(sin(π * x[2]), 0)
# p0(x) = sin(π * (x[1] + x[2]))
# f0(x) = VectorValue(π^2 * μ / ρ * sin(π * x[2]) + π / ρ * cos(π * (x[1] + x[2])), π / ρ * cos(π * (x[1] + x[2])))

# MS 2 
u0(x) = VectorValue(sin(π * x[2]), sin(π * x[1]))
p0(x) = cos(π * x[1] * x[2])
f0(x) = VectorValue(π^2 * μ / ρ * sin(π * x[2]) - π / ρ * x[2] * sin(π*x[1]*x[2]), π^2 * μ / ρ * sin(π * x[1])  - π / ρ * x[1] * sin(π*x[1]*x[2]))




g0 = (u0, p0)
# model, pgs_dict = create_unit_box(0.5, false)
# babuska_stokes_solver(model, pgs_dict, f0, g0; write = true)
error_conv(babuska_stokes_solver, f0, g0)
