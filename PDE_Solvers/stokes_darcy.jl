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


"""
      ---ΓS---
ΓS → |   ΩS   | ← ΓS
      ----Γ---
ΓD → |   ΩD   | ← ΓD
      ---ΓD---

"""
function stokes_darcy_solver(model, pgs_dict, f0, g0, params; write = false)
    # Unpack input
    fs0, fd0 = f0
    us0, ps0, pd0, gΓ = g0
    μ, Κ, α = [params[key] for key in [:μ, :Κ, :α]]

    # --- Boundary tags --- #
    ΩS_tags =  pgs_tags(pgs_dict, [200])
    ΩD_tags =  pgs_tags(pgs_dict, [100])
    ΓS_tags = pgs_tags(pgs_dict, [5, 6, 7, 10, 11, 12, 13])
    ΓD_tags = pgs_tags(pgs_dict, [1, 2, 3, 8, 9])

    
    # --- Triangulation and spaces --- #
    Ω = Triangulation(model)
    ΩS = Triangulation(model, tags=ΩS_tags)
    ΩD = Triangulation(model, tags=ΩD_tags)
    Γ = InterfaceTriangulation(ΩS, ΩD)


    ΓS = BoundaryTriangulation(model, tags=ΓS_tags)
    ΓD = BoundaryTriangulation(model, tags=ΓD_tags)

    order = 2
    ref_us = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    ref_ps = ReferenceFE(lagrangian, Float64, order - 1)
    ref_pd = ReferenceFE(lagrangian, Float64, order)

    # Test spaces
    δVS = TestFESpace(ΩS, ref_us, conformity=:H1, dirichlet_tags=ΓS_tags)
    δQS = TestFESpace(ΩS, ref_ps, conformity=:H1) 
    δQD = TestFESpace(ΩD, ref_pd, conformity=:H1, dirichlet_tags=ΓD_tags)

    # Trial spaces
    VS = TrialFESpace(δVS, us0)
    QS = TrialFESpace(δQS)
    QD = TrialFESpace(δQD, pd0)

    δW = MultiFieldFESpace([δVS, δQS, δQD])
    W = MultiFieldFESpace([VS, QS, QD])    

    # Integration 
    degree = 2 * order
    dΩS = Measure(ΩS, degree)
    dΩD = Measure(ΩD, degree)
    dΓ = Measure(Γ, degree)
    dΓS = Measure(ΓS, degree)
    dΓD = Measure(ΓD, degree)

    # Normal and tangential vectors
    n̂Γ = get_normal_vector(Γ) # Should point from stokes into Darcy
    n̂ΓS = get_normal_vector(ΓS)
    t̂Γ = TensorValue(0, -1, 1, 0) ⋅ n̂Γ.⁺ 

    # --- Weak formulation --- #
    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))

    aΩs((us, ps), (vs, qs)) = ∫( 2*μ * ϵ(us) ⊙ ϵ(vs) - (∇⋅us)*qs - ps*(∇⋅vs) )dΩS
    aΩD((pd, qd)) = ∫( Κ*(∇(pd)⋅∇(qd)) )dΩD
    aΓ((us, pd), (vs, qd)) = ∫( α*(us.⁺⋅t̂Γ)*(vs.⁺⋅t̂Γ) - (us.⁺⋅n̂Γ.⁺)*qd.⁻ + pd.⁻*(n̂Γ.⁺⋅vs.⁺) )dΓ
    a((us, ps, pd), (vs, qs, qd)) =  aΩs((us, ps), (vs, qs)) + aΩD((pd, qd)) + aΓ((us, pd), (vs, qd))
    b((vs, qs, qd)) = ∫( fs0⋅vs )dΩS + ∫( fd0 * qd )dΩD - ∫(CellField(gΓ, Γ) * qd.⁻)dΓ # + neumann 


    # --- Solve --- #
    op = AffineFEOperator(a, b, W, δW)
    ush, psh, pdh = solve(op)
   
    # --- Check with manufactured solution --- #
    us_error = us0 - ush
    ps_error = ps0 - psh
    pd_error = pd0 - pdh

    us_l2norm = sqrt(sum(∫(us_error ⋅ us_error) * dΩS))
    ps_l2norm = sqrt(sum(∫(ps_error ⋅ ps_error) * dΩS))
    pd_l2norm = sqrt(sum(∫(pd_error ⋅ pd_error) * dΩD))


    @printf("us: l2 norm = %e \n", us_l2norm)
    @printf("ps: l2 norm = %e \n", ps_l2norm)
    @printf("pd: l2 norm = %e \n", pd_l2norm)

    write && writevtk(Ω, path * "vtu_files/" * "stokes_darcy_results", cellfields=["us" => ush, "ps" => psh, "pd" => pdh, "us_error" => us_error, "ps_error" => ps_error, "pd_error" => pd_error])
    
    return us_l2norm, ps_l2norm, pd_l2norm
    

end


function error_conv(solver, f0, g0, params; lc_start=2, num_points=5, show_plot=false)
    l2norm = zeros(Float64, num_points, 3)
    lc = zeros(Float64, num_points)
    slope = zeros(Float64, 3)


    # --- Decrease mesh size and get l2norm --- #
    write = false
    for p in 1:num_points
        lc[p] = lc_start * (1 / 2)^(p - 1)
        p == num_points && (write = true)
        model, pgs_dict = create_coupled_box(lc[p])
        l2norm[p, :] .= solver(model, pgs_dict, f0, g0, params; write)

        # Get the actual mesh size
        Ω = Triangulation(model)
        lc[p] = minimum(sqrt.(collect(get_cell_measure(Ω))))
    end

    println("mesh size: ", lc)
    println("us: l2 norm = ", l2norm[:, 1])
    println("ps: l2 norm = ", l2norm[:, 2])
    println("pd: l2 norm = ", l2norm[:, 3])


    # ---  Convergence rate --- #
    X = ones(num_points, 2)
    X[:, 2] = log.(lc)
    for i in 1:3
        y = log.(l2norm[:, i])
        b = inv(transpose(X) * X) * transpose(X) * y
        slope[i] = b[2]
      
    end

    println("us: Convergence rate: ", slope[1])
    println("ps: Convergence rate: ", slope[2])
    println("pd: Convergence rate: ", slope[3])

    if show_plot
        fig = plot(lc, l2norm[:, 1], xaxis=:log, yaxis=:log, label="us", marker=:x)
        fig = plot!(lc, l2norm[:, 2], xaxis=:log, yaxis=:log, label="ps", marker=:o)
        fig = plot!(lc, l2norm[:, 3], xaxis=:log, yaxis=:log, label="pd", marker=:o)

        display(fig)
    end

end



#MS1
μ = 1.0
Κ = 1.0 
α(x) = μ*π * (cos(π * x[2]) - sin(π * x[1]))/sin(π * x[2])
us0(x) = VectorValue(sin(π * x[2]), cos(π * x[1]))
ps0(x) = -sin(π * x[2])
fs0(x) = VectorValue(μ*π^2 * sin(π * x[2]), μ*π^2 * cos(π * x[1]) - π*cos(π * x[2]) )
σ(x) = TensorValue(-p0(x), μ*π * (cos(π * x[2]) - sin(π * x[1])), μ*π * (cos(π * x[2]) - sin(π * x[1])), -p0(x))
pd0(x) = ps0(x)
gΓ(x) = -cos(π * x[1]) + Κ*π*cos(π*x[2])
fd0(x) = -Κ * π^2 * sin(π * x[2])


# # MS4
# μ = 1.0
# Κ = 1.0 
# α(x) = μ*π*(cos(π*x[1]) + cos(π*x[2]))/sin(π*x[2])
# us0(x) = VectorValue(sin(π*x[2]), sin(π*x[1]))
# ps0(x) = sin(π*x[2])
# fs0(x) = VectorValue(μ*π^2*sin(π*x[2]), μ*π^2*sin(π*x[1]) + π*cos(π*x[2]))
# pd0(x) = 2*x[2] 
# gΓ(x) = -sin(π*x[1]) - 2*Κ
# fd0(x) = 0.0





params = Dict(:μ => μ, :Κ => Κ, :α => α)
f0 = (fs0, fd0)
g0 = (us0, ps0, pd0, gΓ)
# model, pgs_dict = create_coupled_box(2, false)
# stokes_darcy_solver(model, pgs_dict, f0, g0, params; write = true)
error_conv(stokes_darcy_solver, f0, g0, params; show_plot=false)





