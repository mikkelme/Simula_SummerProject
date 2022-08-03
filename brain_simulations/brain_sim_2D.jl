using Gridap
using GridapGmsh
using Gridap.CellData
using Printf
using Plots


include("../mesh_generators/create_brain.jl")



path = "/Users/mikkelme/Documents/Github/Simula_SummerProject/brain_simulations/"
if !ispath(path)
    path = "/home/mirok/Documents/MkSoftware/Simula_SummerProject/brain_simulations/"
end


mutable struct PDE_params
  μ::Float64  
  Κ::Float64   
  α::Function
  α_body::String
  ps0::Function
  ps0_body::String
  ∇pd0::Function
  ∇pd0_body::String
  function PDE_params(μ, Κ, α_body, ps0_body, ∇pd0_body)
      new(μ, Κ, (x -> Base.invokelatest(eval(Meta.parse(α_body)), x)), α_body, (x -> Base.invokelatest(eval(Meta.parse(ps0_body)), x)), ps0_body, (x -> Base.invokelatest(eval(Meta.parse(∇pd0_body)), x)), ∇pd0_body)
  end
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
    ΓD_neutags =  pgs_tags(pgs_dict, [3, 6, 8, 10, 11, 13, 14]) 
    ΛS_diritags = pgs_tags(pgs_dict, [5, 12, 15]) 

    
    # --- Triangulation and spaces --- #
    # Surface domains
    Ω = Triangulation(model)
    ΩS = Triangulation(model, tags=ΩS_tags)
    ΩD = Triangulation(model, tags=ΩD_tags)

    # Boundary & interface
    Γ = InterfaceTriangulation(ΩS, ΩD)   
    ΛS = BoundaryTriangulation(model, tags=ΛS_tags)
    ΓS = BoundaryTriangulation(model, tags=ΓS_tags)
    ΓD = BoundaryTriangulation(model, tags=ΓD_tags)
    ΓD_neu = BoundaryTriangulation(model, tags=ΓD_neutags)
    
    # Reference elementes
    order = 2
    ref_us = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    ref_ps = ReferenceFE(lagrangian, Float64, order - 1)
    ref_pd = ReferenceFE(lagrangian, Float64, order)
    
    # Test spaces
    δVS = TestFESpace(ΩS, ref_us, conformity=:H1, dirichlet_tags=ΛS_diritags)
    δQS = TestFESpace(ΩS, ref_ps, conformity=:H1) 
    δQD = TestFESpace(ΩD, ref_pd, conformity=:H1)
    
    # Trial spaces
    VS = TrialFESpace(δVS, VectorValue(0.0, 0.0)) # No slip on ΛS
    QS = TrialFESpace(δQS)
    QD = TrialFESpace(δQD)
    
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
    dΓDneu = Measure(ΓD_neu, degree) 
    
    # Normal & tangential vectors
    n̂Γ = get_normal_vector(Γ) 
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
    σ(u,p) = 2 * data.μ * ε(u) - p * TensorValue(1, 0, 0, 1) # Stress matrix

    # Nitsche on ΓS
    aNΓS((us, ps), (vs, qs)) = ∫(γ/h * (us⋅t̂ΓS) * (vs ⋅ t̂ΓS))dΓS  - ∫( ((n̂ΓS ⋅ σ(us,ps)) ⋅ t̂ΓS) * (vs⋅t̂ΓS) )dΓS - ∫( ((n̂ΓS ⋅ σ(vs,qs)) ⋅ t̂ΓS) * (us⋅t̂ΓS) )dΓS
    bNΓS((vs, qs)) = ∫(data.ps0 * (-vs ⋅ n̂ΓS))dΓS 
    
    # Stokes domain (left hand side)
    aΩs((us, ps), (vs, qs)) = ∫( 2*data.μ * ϵ(us) ⊙ ϵ(vs) - (∇⋅us)*qs - ps*(∇⋅vs) )dΩS

    # Darcy Domain (left hand side)
    aΩD((pd, qd)) = ∫( data.Κ/data.μ*(∇(pd)⋅∇(qd)) )dΩD

    # Interface coupling (left hand side)
    aΓ((us, pd), (vs, qd)) = ∫( data.α*(us.⁺⋅t̂Γ)*(vs.⁺⋅t̂Γ) - (us.⁺⋅n̂Γ.⁺)*qd.⁻ + pd.⁻*(n̂Γ.⁺⋅vs.⁺) )dΓ

    # Neumann conditions on Darcy (right hand side)
    b_neumann((vs, qd)) = ∫( data.Κ/data.μ*(get_normal_vector(ΓD_neu)⋅data.∇pd0)*qd )dΓDneu

    # Gathering terms
    a((us, ps, pd), (vs, qs, qd)) =  aΩs((us, ps), (vs, qs)) + aΩD((pd, qd)) + aΓ((us, pd), (vs, qd)) + aNΓS((us, ps), (vs, qs))
    b((vs, qs, qd)) = bNΓS((vs, qs)) + b_neumann((vs, qd)) # +∫( data[:fs0]⋅vs )dΩS + ∫( data[:fd0] * qd )dΩD
        
    
    # --- Solve --- #
    println("#--- Assembling ---#")
    op = AffineFEOperator(a, b, W, δW)
    println("#--- Solving ---#")
    @show size(op.op.matrix)
    ush, psh, pdh = solve(op)

    
    # --- Write & return results --- #
    write!=false && writevtk(ΩS, write[1] * "stokes_" * write[2], cellfields=["us" => ush, "psh" => psh])
    write!=false && writevtk(ΩD, write[1] * "darcy_" * write[2], cellfields=["pdh" => pdh])

    return  ush, psh, pdh, ΩS, ΩD, Γ
    

end


# function evaluate_along_centerline()


#   # f(x) = x[1] + x[2]
#   domain = (-0.5, 0.5, 9.5, 9.6)
#   partition = (5,5)
#   box = CartesianDiscreteModel(domain, partition)
#   B = BoundaryTriangulation(box)
#   # reffe₁ = ReferenceFE(lagrangian, Float64, 1)
#   # V₁ = FESpace(𝒯₁, reffe₁)
  

#   # fₕ = interpolate_everywhere(f,V₁)



#   iu = Interpolable(ush)
#   ip = Interpolable(psh; searchmethod=KDTreeSearch(num_nearest_vertices=4))


#   ip(VectorValue(-1.93, 9.3))
#   ph(VectorValue(-1.93, 9.3))



#   cl_model, cl_pgs_dict = create_centerline(brain_params; view = true)
#   # model = GmshDiscreteModel("./foo.msh") # Test the mesh
#   cl_tags =  pgs_tags(cl_pgs_dict, [1, 2, 3])
#   L = Triangulation(cl_model, tags=pgs_tags(cl_pgs_dict, [2, 5]))
#   # L = B
#   #L = Triangulation(box)
#   dL = Measure(L, 1)

#   Vi = TestFESpace(L, ReferenceFE(lagrangian, Float64, 1))
#   #phL = interpolate_everywhere(ip, Vi)
#   val = sum(∫(1)*dL)


#   # test = sqrt(sum(∫(ush ⋅ ush) * dΩS))
#   test = sum(∫(ip)dΩS)


#   # gp = Gridap.FESpaces.interpolate_everywhere(ip, V2)
#   # gp = interpolate(ip, V2)


#   # integral = sqrt(sum(∫(ush ⋅ ush)dL))
#   # integral = sqrt(sum(∫(iu ⋅ iu) * dL))
#   # val  = ∫( px )dL

# end




# Brain Model [length unit: meter] 
# lc = 2e-4 # size(op.op.matrix) = (425840, 425840)
lc = 1e-3
arcLen = (100e-3, 0)
r_brain = 10e-3  
d_ratio = 1.5e-3/r_brain
r_curv = 50e-3 
A = 1e-3
λ = 10*1e-3
ω(λ) = 2*pi/λ      
inner_perturb = @sprintf("(x,z) -> %f * sin(abs(x) * %f - pi/2) * fld(mod2pi(abs(x) * %f - pi/2),pi) ", A , ω(λ), ω(λ))
outer_perturb = "(x,z) -> 0.0"  
BS_points = (3000, 3000) 
field_Lc_lim = [1 / 2, 1]
field_Dist_lim = [1e-3, 5e-3] 

# PDE parameters 
μ = 0.8e-3  # Cerobrospinal fluid viscosity [Pa * s]
Κ = 1e-16   # Permeability in brain parenchyma [m^2] 
α = "(x) -> 1*μ/sqrt(Κ)" # Slip factor on Γ [Pa * s / m]
ps0 = "(x) -> x[1] < 0 ? 1*133.3224 : 0." # 1*mmHg [Pa]
∇pd0 = "(x) -> VectorValue(0.0, 0.0)" # Zero flux




# --- Run simulation --- #
brain_param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)
PDE_param = PDE_params(μ, Κ, α, ps0, ∇pd0)
model, pgs_dict = create_brain(brain_param; view=false, write=false)
ush, psh, pdh = brain_PDE(model, pgs_dict, PDE_param; write = (path * "vtu_files/", "2D_default"))


include("../mesh_generators/generate_radial_lines.jl")
include("../mesh_generators/generate_center_line.jl")
num_radial_lines = 100
rad_model, _ =  create_radial_lines(brain_param, num_radial_lines; view = true)
writevtk(rad_model, path * "vtu_files/radial_lines" )

S_center_model, _ = create_centerline(brain_param, "S", view = true)
D_center_model, _ = create_centerline(brain_param, "D")

SL = Triangulation(S_center_model)
DL = Triangulation(D_center_model)



writevtk(SL, path * "vtu_files/S_center_model" )
writevtk(DL, path * "vtu_files/D_center_model" )


