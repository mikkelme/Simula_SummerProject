using Gridap
using GridapGmsh
using Printf
using Plots


include("../mesh_generators/create_brain.jl")


path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/PDE_Solvers/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/PDE_Solvers/"
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
   
    
    # Nitsche on ΓS
    aNΓS((us, ps), (vs, qs)) = ∫(γ/h * (us⋅t̂ΓS) * (vs ⋅ t̂ΓS))dΓS  - ∫( ((n̂ΓS ⋅ σ(us,ps)) ⋅ t̂ΓS) * (vs⋅t̂ΓS) )dΓS - ∫( ((n̂ΓS ⋅ σ(vs,qs)) ⋅ t̂ΓS) * (us⋅t̂ΓS) )dΓS
    bNΓS((vs, qs)) = ∫(data[:ps0] * (-vs ⋅ n̂ΓS))dΓS + ∫(data[:gΓS] * (γ/h * (vs ⋅ t̂ΓS) - ((n̂ΓS ⋅ σ(vs,qs)) ⋅ t̂ΓS)))dΓS
    
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

    #try ph = psh + pdh better visualization
    
    # --- Check with manufactured solution --- #
    write && writevtk(Ω, path * "vtu_files/" * "brain_sim_results", cellfields=["us" => ush, "ps" => psh, "pd" => pdh])
    
    return us_l2norm, ps_l2norm, pd_l2norm
    

end



# --- Brain Model --- # 
lc = 0.5
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
μ = 1.0
Κ = 1.0 
α(x) = 0.0
gΓS(x) = 0.0
us0(x) = VectorValue(0.0, 0.0) # us = 0 on ΛS
ps0(x) = 0.0
pd0(x) = 0.0
fs0(x) = VectorValue(0.0, 0.0) #?
fd0(x) = 0.0 #?

σ0(x) = TensorValue(0.0, 0.0, 0.0, 0.0)
nab_pd0(x) = VectorValue(0.0, 0.0)
gΓ(x) = 0 

PDE_params = Dict(:μ => μ, :Κ => Κ, :α => α, :gΓS => gΓS, :fs0 => fs0, :fd0 => fd0, :us0 => us0, :ps0 => ps0, :pd0 => pd0, :gΓ => gΓ, :σ0 => σ0, :nab_pd0 => nab_pd0) 


# PDE_params = Dict(:fs0 => fs0, :gΓ => gΓ)

# f0 = (fs0, fd0)
# g0 = (us0, ps0, pd0, gΓ)
# h0 = (σ0, nab_pd0)



# fs0, fd0 = f0
# us0, ps0, pd0, gΓ = g0
# σ0, nab_pd0 = h0
# μ, Κ, α = [params[key] for key in [:μ, :Κ, :α]]

PDE_keys = collect(keys(PDE_params))
key = PDE_keys[1]
test = PDE_params[key]
println(test)

# brain_PDE(model, pgs_dict, params; write = false)





