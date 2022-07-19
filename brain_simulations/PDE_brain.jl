using Gridap
using GridapGmsh
using Printf
using Plots


include("../mesh_generators/create_brain.jl")
include("../mesh_generators/unit_box_direct.jl")



path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/brain_simulations/"
end


"""
Domain notation:

      ---ΛS---
ΓS → |   ΩS   | ← ΓS
      ----Γ---
ΓD → |   ΩD   | ← ΓD
      ---ΓD---

Physical tags:

    (12)--5--(15)
  7→ |   <2>  | ←9
    (11)--4--(14)
  6→ |   <1>  | ←8
    (10)--3--(13)

"""
function brain_PDE(model, pgs_dict, data; write = false)
    # --- Boundary tags --- #
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




    # # Unit box 
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
    # #####

    
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
    
    # writevtk(ΩS, path*"model_S")
    # writevtk(ΩD, path*"model_D")
    # writevtk(Γ, path*"model_Gamma")
    # return


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
    VS = TrialFESpace(δVS, data[:us0])
    QS = TrialFESpace(δQS)
    QD = TrialFESpace(δQD, data[:pd0])
    
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
    tangent = TensorValue(0, -1, 1, 0) ⋅ n̂Γ.⁺ 
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
    σ(u,p) = 2 * data[:μ] * ε(u) - p * TensorValue(1, 0, 0, 1) # Stress matrix
    gΓS = us0 × n̂ΓS # REMEBER TO REMOVE THIS 
    
    # Nitsche on ΓS
    aNΓS((us, ps), (vs, qs)) = ∫(γ/h * (us⋅t̂ΓS) * (vs ⋅ t̂ΓS))dΓS  - ∫( ((n̂ΓS ⋅ σ(us,ps)) ⋅ t̂ΓS) * (vs⋅t̂ΓS) )dΓS - ∫( ((n̂ΓS ⋅ σ(vs,qs)) ⋅ t̂ΓS) * (us⋅t̂ΓS) )dΓS
    bNΓS((vs, qs)) = ∫(data[:ps0] * (-vs ⋅ n̂ΓS))dΓS + ∫(gΓS * (γ/h * (vs ⋅ t̂ΓS) - ((n̂ΓS ⋅ σ(vs,qs)) ⋅ t̂ΓS)))dΓS
    
    # Stokes domain (left hand side)
    aΩs((us, ps), (vs, qs)) = ∫( 2*data[:μ] * ϵ(us) ⊙ ϵ(vs) - (∇⋅us)*qs - ps*(∇⋅vs) )dΩS

    # Darcy Domain (left hand side)
    aΩD((pd, qd)) = ∫( data[:Κ]*(∇(pd)⋅∇(qd)) )dΩD

    # Interface coupling (left hand side)
    aΓ((us, pd), (vs, qd)) = ∫( data[:α]*(us.⁺⋅t̂Γ)*(vs.⁺⋅t̂Γ) - (us.⁺⋅n̂Γ.⁺)*qd.⁻ + pd.⁻*(n̂Γ.⁺⋅vs.⁺) )dΓ

    # Neumann conditions (right hand side)
    b_neumann((vs, qd)) = ∫( (data[:σ0] ⋅ get_normal_vector(ΛS_neu) ) ⋅ vs)dΛSneu + ∫( data[:Κ]*(get_normal_vector(ΓD_neu)⋅data[:nab_pd0])*qd )dΓDneu

    # Gathering terms
    a((us, ps, pd), (vs, qs, qd)) =  aΩs((us, ps), (vs, qs)) + aΩD((pd, qd)) + aΓ((us, pd), (vs, qd)) + aNΓS((us, ps), (vs, qs))
    b((vs, qs, qd)) = ∫( data[:fs0]⋅vs )dΩS + ∫( data[:fd0] * qd )dΩD - ∫(CellField(data[:gΓ], Γ) * qd.⁻)dΓ + bNΓS((vs, qs)) + b_neumann((vs, qd))
        
    
    # --- Solve --- #
    println("Setting up matrix and vector")
    op = AffineFEOperator(a, b, W, δW)
    println("Solving")
    ush, psh, pdh = solve(op)

    ph = psh + pdh 
    
    # --- Check with manufactured solution --- #
    write && writevtk(Ω, path * "vtu_files/" * "brain_sim_results", cellfields=["us" => ush, "ph" => ph])
    
    return  ush, psh, pdh
    

end



# --- Brain Model --- # 
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
model, pgs_dict = create_brain(brain_params; view=false, write=false)


# --- PDE parameters --- #
# μ = 1.0
# Κ = 1.0 
# α(x) = 0.0
# gΓS(x) = 0.0
# us0(x) = VectorValue(0.0, 0.0) # us = 0 on ΛS
# ps0(x) = -x # use angle instead?
# pd0(x) = 0.0
# fs0(x) = VectorValue(0.0, 0.0) #?
# fd0(x) = 0.0 #?

# σ0(x) = TensorValue(0.0, 0.0, 0.0, 0.0) # not used right now
# nab_pd0(x) = VectorValue(0.0, 0.0) # zero flux?
# gΓ(x) = 0 


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

PDE_params = Dict(:μ => μ, :Κ => Κ, :α => α, :gΓS => gΓS, :fs0 => fs0, :fd0 => fd0, :us0 => us0, :ps0 => ps0, :pd0 => pd0, :gΓ => gΓ, :σ0 => σ0, :nab_pd0 => nab_pd0) 


# --- Run simulation --- #
# model, pgs_dict  = create_coupled_box(1, false)
brain_PDE(model, pgs_dict, PDE_params; write = true)





