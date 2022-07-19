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
Domain notation:

      ---ΛS---
ΓS → |   ΩS   | ← ΓS
      ----Γ---
ΓD → |   ΩD   | ← ΓD
      ---ΓD---

Physical tags:

   (12)--7--(13)
  5→ |  200  | ←6
   (10)--4--(11)
  2→ |  100  | ←3
    (8)--1--(9)

"""
function stokes_darcy_solver(model, pgs_dict, f0, g0, h0, params; write = false)
    # Unpack input
    fs0, fd0 = f0
    us0, ps0, pd0, gΓ = g0
    σ0, nab_pd0 = h0
    μ, Κ, α = [params[key] for key in [:μ, :Κ, :α]]



    # Brain
    ΩS_tags =  pgs_tags(pgs_dict, [2])
    ΩD_tags =  pgs_tags(pgs_dict, [1])
    ΛS_tags = pgs_tags(pgs_dict, [5, 12, 15]) 
    ΓD_tags = pgs_tags(pgs_dict, [3, 6, 8, 10, 11, 13, 14])
    ΓS_tags = pgs_tags(pgs_dict, [7, 9]) 

    # Boundary conditions
    ΛS_neutags = pgs_tags(pgs_dict, []) 
    ΓD_neutags = pgs_tags(pgs_dict, [3, 6, 8, 10, 11, 13, 14]) 
    ΛS_diritags = filter(x -> x ∉ ΛS_neutags, ΛS_tags)
    ΓD_diritags = filter(x -> x ∉ ΓD_neutags, ΓD_tags)
    ####

    # # --- Boundary tags --- #
    # ΩS_tags =  pgs_tags(pgs_dict, [200])
    # ΩD_tags =  pgs_tags(pgs_dict, [100])
    # ΛS_tags = pgs_tags(pgs_dict, [7, 12, 13]) 
    # ΓD_tags = pgs_tags(pgs_dict, [1, 2, 3, 8, 9, 10, 11])
    # ΓS_tags = pgs_tags(pgs_dict, [5, 6]) 

    # # Boundary conditions
    # ΛS_neutags = pgs_tags(pgs_dict, []) 
    # ΓD_neutags = pgs_tags(pgs_dict, [1, 2, 3, 8, 9, 10, 11]) 
    # ΛS_diritags = filter(x -> x ∉ ΛS_neutags, ΛS_tags)
    # ΓD_diritags = filter(x -> x ∉ ΓD_neutags, ΓD_tags)

    
    # --- Triangulation and spaces --- #
    # Surface domains
    Ω = Triangulation(model)
    ΩS = Triangulation(model, tags=ΩS_tags)
    ΩD = Triangulation(model, tags=ΩD_tags)

    # Boundary and interface
    Γ = InterfaceTriangulation(ΩS, ΩD)   
    ΛS = BoundaryTriangulation(model, tags=ΛS_tags)
    ΓS = BoundaryTriangulation(model, tags=ΓS_tags)
    ΓD = BoundaryTriangulation(model, tags=ΓD_tags)
    ΛS_neu = BoundaryTriangulation(model, tags=ΛS_neutags)
    ΓD_neu = BoundaryTriangulation(model, tags=ΓD_neutags)
    
    # Reference elementes
    order = 2
    ref_us = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    ref_ps = ReferenceFE(lagrangian, Float64, order - 1)
    ref_pd = ReferenceFE(lagrangian, Float64, order)
    
    # Test spaces
    δVS = TestFESpace(ΩS, ref_us, conformity=:H1, dirichlet_tags=ΛS_diritags)
    δQS = TestFESpace(ΩS, ref_ps, conformity=:H1) 
    δQD = TestFESpace(ΩD, ref_pd, conformity=:H1, dirichlet_tags=ΓD_diritags)
    
    # Trial spaces
    VS = TrialFESpace(δVS, us0)
    QS = TrialFESpace(δQS)
    QD = TrialFESpace(δQD, pd0)
    
    # Multispace
    δW = MultiFieldFESpace([δVS, δQS, δQD]) # Test
    W = MultiFieldFESpace([VS, QS, QD]) # Train
    
    # --- Integration ---- #
    degree = 2 * order
    dΩS = Measure(ΩS, degree)
    dΩD = Measure(ΩD, degree)
    dΓ = Measure(Γ, degree)
    dΛS = Measure(ΛS, degree)
    dΓS = Measure(ΓS, degree)
    dΓD = Measure(ΓD, degree)    
    dΛSneu = Measure(ΛS_neu, degree) 
    dΓDneu = Measure(ΓD_neu, degree) 
    
    # Normal and tangential vectors
    n̂Γ = get_normal_vector(Γ) 
    n̂ΛS = get_normal_vector(ΛS)
    n̂ΓS = get_normal_vector(ΓS)
    n̂D = get_normal_vector(ΓD)
    t̂Γ = TensorValue(0, -1, 1, 0) ⋅ n̂Γ.⁺ 
    t̂ΓS = TensorValue(0, -1, 1, 0) ⋅ n̂ΓS 

    
    # Nitsche 
    γ = 10^order # Nitsche penalty parameter
    h_ΓS = get_array(∫(1) * dΓS)
    h = CellField(lazy_map(h -> h, h_ΓS), ΓS) # Element size
    
    # --- Weak formulation --- #
    ϵ(u) = 1 / 2 * (∇(u) + transpose(∇(u)))
    σ(u,p) = 2 * μ * ε(u) - p * TensorValue(1, 0, 0, 1) # Stress matrix
    gΓS = us0 × n̂ΓS # = 0 on brain simulations
    
    # Nitsche on ΓS
    aNΓS((us, ps), (vs, qs)) = ∫(γ/h * (us⋅t̂ΓS) * (vs ⋅ t̂ΓS))dΓS  - ∫( ((n̂ΓS ⋅ σ(us,ps)) ⋅ t̂ΓS) * (vs⋅t̂ΓS) )dΓS - ∫( ((n̂ΓS ⋅ σ(vs,qs)) ⋅ t̂ΓS) * (us⋅t̂ΓS) )dΓS
    bNΓS((vs, qs)) = ∫(ps0 * (-vs ⋅ n̂ΓS))dΓS + ∫(gΓS * (γ/h * (vs ⋅ t̂ΓS) - ((n̂ΓS ⋅ σ(vs,qs)) ⋅ t̂ΓS)))dΓS
    
    # Stokes domain (left hand side)
    aΩs((us, ps), (vs, qs)) = ∫( 2*μ * ϵ(us) ⊙ ϵ(vs) - (∇⋅us)*qs - ps*(∇⋅vs) )dΩS

    # Darcy Domain (left hand side)
    aΩD((pd, qd)) = ∫( Κ*(∇(pd)⋅∇(qd)) )dΩD

    # Interface coupling (left hand side)
    aΓ((us, pd), (vs, qd)) = ∫( α*(us.⁺⋅t̂Γ)*(vs.⁺⋅t̂Γ) - (us.⁺⋅n̂Γ.⁺)*qd.⁻ + pd.⁻*(n̂Γ.⁺⋅vs.⁺) )dΓ

    # Neumann conditions (right hand side)
    b_neumann((vs, qd)) = ∫( (σ0 ⋅ get_normal_vector(ΛS_neu) ) ⋅ vs)dΛSneu + ∫( Κ*(get_normal_vector(ΓD_neu)⋅nab_pd0)*qd )dΓDneu

    # Gathering terms
    a((us, ps, pd), (vs, qs, qd)) =  aΩs((us, ps), (vs, qs)) + aΩD((pd, qd)) + aΓ((us, pd), (vs, qd)) + aNΓS((us, ps), (vs, qs))
    b((vs, qs, qd)) = ∫( fs0⋅vs )dΩS + ∫( fd0 * qd )dΩD - ∫(CellField(gΓ, Γ) * qd.⁻)dΓ + bNΓS((vs, qs)) + b_neumann((vs, qd))
        
    
    # --- Solve --- #
    println("Setting up matrix and vector")
    op = AffineFEOperator(a, b, W, δW)
    println("Solving")
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

    write && writevtk(Ω, path * "vtu_files/" * "nitsche_stokes_darcy_results", cellfields=["us" => ush, "ps" => psh, "pd" => pdh, "us_error" => us_error, "ps_error" => ps_error, "pd_error" => pd_error])
    
    return us_l2norm, ps_l2norm, pd_l2norm
    

end


function error_conv(solver, f0, g0, h0, params; lc_start=2, num_points=5, show_plot=false)
    l2norm = zeros(Float64, num_points, 3)
    lc = zeros(Float64, num_points)
    slope = zeros(Float64, 3)

    # --- Decrease mesh size and get l2norm --- #
    write = false
    for p in 1:num_points
        lc[p] = lc_start * (1 / 2)^(p - 1)
        p == num_points && (write = true)
        model, pgs_dict = create_brain(brain_params; view=false, write=false)
        # model, pgs_dict = create_coupled_box(lc[p])
        l2norm[p, :] .= solver(model, pgs_dict, f0, g0, h0, params; write)

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


# # MS1 (uₛ × n̂ₛ = 0 on ΓS)
# μ = 1.0
# Κ = 1.0 
# α(x) = μ*π * (cos(π*x[1]) + cos(π*x[2])) / sin(π*x[2])
# us0(x) = VectorValue(sin(π * x[2]), sin(π * x[1]))
# ps0(x) = sin(π * x[2])
# fs0(x) = VectorValue(μ*π^2 * sin(π * x[2]), μ*π^2 * sin(π * x[1]) + π*cos(π * x[2]) )
# pd0(x) = ps0(x)
# gΓ(x) = -sin(π*x[1]) - Κ*π*cos(π*x[2])
# fd0(x) = Κ * π^2 * sin(π * x[2])

# MS2 (uₛ × n̂ₛ ≠ 0 on ΓS)
μ = 1.0
Κ = 1.0 
us0(x) = VectorValue(sin(π * x[2]), cos(π * x[1])) 
ps0(x) = -sin(π * x[2])
pd0(x) = ps0(x)
fs0(x) = VectorValue(μ*π^2 * sin(π * x[2]), μ*π^2 * cos(π * x[1]) - π*cos(π * x[2]) )
fd0(x) = -Κ * π^2 * sin(π * x[2])

σ0(x) = TensorValue(-ps0(x), μ*π * (cos(π * x[2]) - sin(π * x[1])), μ*π * (cos(π * x[2]) - sin(π * x[1])), -ps0(x))
nab_pd0(x) = VectorValue(0, -π*cos(π * x[2]))

α(x) = μ*π * (cos(π * x[2]) - sin(π * x[1]))/sin(π * x[2])
gΓ(x) = -cos(π * x[1]) + Κ*π*cos(π*x[2])



params = Dict(:μ => μ, :Κ => Κ, :α => α)
f0 = (fs0, fd0)
g0 = (us0, ps0, pd0, gΓ)
h0 = (σ0, nab_pd0)



lc = 1.0
arcLen = (5, 0)
r_brain = 2
d_ratio = 0.5
r_curv = 50
inner_perturb(x, y) = 0
outer_perturb(x, y) = 0
BS_points = (arcLen[1]*20, arcLen[2]*10)
field_Lc_lim = [1 / 2, 1]
field_Dist_lim = [0.1, 0.5]

brain_params = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
# model, pgs_dict = create_brain(brain_params; view=false, write=false)



# model, pgs_dict = create_coupled_box(2, false)
# stokes_darcy_solver(model, pgs_dict, f0, g0, h0, params; write = false)
error_conv(stokes_darcy_solver, f0, g0, h0, params; show_plot=false)





