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


function pgs_tags(pgs_dict, tags)
    return [pgs_dict[tag] for tag in tags]
end




function babuska_stokes_solver(model, pgs_dict, f0, g0, write=false)

    # Boundary conditions
    λ_diri_tags = pgs_tags(pgs_dict, [3]) # goes into Multiplier 
    λ_corner_tags = pgs_tags(pgs_dict, [7, 8])
    strong_diri_tags = pgs_tags(pgs_dict, [1, 2, 4, 5, 6, 7, 8]) # goes intro trial space
    neumann_tags = Vector{Int}()


    # Triangulation and spaces
    Ω = Triangulation(model)
    ΓD = BoundaryTriangulation(model, tags=λ_diri_tags)
    ΓN = BoundaryTriangulation(model, tags=neumann_tags)

    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order - 1, space=:P)
    reffeₘ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)

    δV = TestFESpace(Ω, reffeᵤ, conformity=:H1, dirichlet_tags=strong_diri_tags) # Velocity
    δM = TestFESpace(Ω, reffeₚ, conformity=:L2, constraint=:zeromean) # Pressure
    δΛ = TestFESpace(ΓD, reffeₘ, conformity=:H1, dirichlet_tags=λ_corner_tags) # Multiplier

    V = TrialFESpace(δV, g0)
    M = TrialFESpace(δP)
    Λ = TrialFESpace(δM, VectorValue(0, 0))

    δW = MultiFieldFESpace([δV, δM, δΛ]) # Test spaces
    W = MultiFieldFESpace([V, M, Λ]) # Trial spaces


    degree = 2 * order #
    dΩ = Measure(Ω, degree)
    dΓD = Measure(ΓD, degree)
    dΓN = Measure(ΓN, degree)
    nD = get_normal_vector(ΓD)

    # Weak formulation 
    μ = 1.0
    ρ = 1.0 # Hardcoded for now
    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))

    t(n) = VectorValue(n[2], -n[1])
    c(u, λ) = ∫((u ⋅ t(nD)) * λ) * dΓD

    a((u, p, λ), (v, q, η)) = ∫(2 * μ * ϵ(u) ⊙ ϵ(v) - p * (∇ ⋅ v) ± q * (∇ ⋅ u) + c(v, λ)) * dΩ ± c(u, η)
    # a((u, p, λ), (v, q, η)) = ∫(2 * μ * ϵ(u) ⊙ ϵ(v) - p * (∇ ⋅ v) ± q * (∇ ⋅ u) + c(v, λ)) * dΩ ± c(u, η)
    b((v, q, η)) = ∫(ρ * v ⋅ f0) * dΩ - ∫(p0 * nD ⋅ v) * dΓD


    # --- Solve --- #
    op = AffineFEOperator(a, b, W, δW)
    uh, ph, λh = solve(op)


    # --- Check with manufactured solution --- #
    u_error = g0[1] - uh
    p_error = g0[2] - ph

    u_l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
    p_l2norm = sqrt(sum(∫(p_error ⋅ p_error) * dΩ))

    @printf("u: l2 norm = %e \n", u_l2norm)
    @printf("p: l2 norm = %e \n", p_l2norm)

    wrtie && writevtk(Ωₕ, path * "babuska_stokes_results", order=2, cellfields=["uh" => uh, "ph" => ph, "u_error" => u_error, "p_error" => p_error])

end


# MS (vector)
u0(x) = VectorValue(sin(π * x[2]), cos(π * x[1]))
p0(x) = sin(π * (x[1] + x[2]))
f0(x) = VectorValue(π^2 * sin(π * x[2]) + π * cos(π * (x[1] + x[2])), π^2 * cos(π * x[1]) + π * cos(π * (x[1] + x[2])))
σ(x) = TensorValue(-p0(x), π * (cos(π * x[2]) - sin(π * x[1])), π * (cos(π * x[2]) - sin(π * x[1])), -p0(x))
g0 = (u0, p0)


model, pgs_dict = create_unit_box(1, false)
babuska_stokes_solver(model, pgs_dict, f0, g0, h0)
