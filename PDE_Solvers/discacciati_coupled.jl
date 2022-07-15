



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
      ---Λf---
Γf → |        | ← Γf
      ----Γ---
Γp → |        | ← Γp
      ---Λp---

"""
function(model, pgs_dict)

    # --- Boundary tags --- #
    Λf_tags = pgs_tags(pgs_dict, [7])
    Λp_tags = pgs_tags(pgs_dict, [1])
    Γf_tags = pgs_tags(pgs_dict, [5, 6])
    Γp_tags = pgs_tags(pgs_dict, [2, 3])
    Γ_tags = pgs_tags(pgs_dict, [4])

    
    # --- Triangulation and spaces --- #
    Ω = Triangulation(model)
    Ωf = Triangulation(?)
    Ωp = Triangulation(?)

    Λf = BoundaryTriangulation(model, tags=Λf_tags)
    Λp = BoundaryTriangulation(model, tags=Λp_tags)
    Γf = BoundaryTriangulation(model, tags=Γf_tags)
    Γp = BoundaryTriangulation(model, tags=Γp_tags)
    Γ = BoundaryTriangulation(model, tags=Γ_tags)

    # assign velocity Vf = 0 on Γf
    # assign ϕ = ϕₚ on Λp 
    # assign Vₚ ⋅ n̂ₚ on Λp


    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order - 1)

    δW = TestFESpace(Ω, reffeᵤ, conformity=:H1, dirichlet_tags = Γf_tags + Γp_tags)
    δQ = TestFESpace(Ωf, reffeₚ, conformity=:L2)

    W = TrialFESpace(δW, _?_)
    Q = TrialFESpace(δQ)


    degree = 2 * order
    dΩ = Measure(Ω, degree)
    dΩf = Measure(Ωf, degree)
    dΩp = Measure(Ωp, degree)

    dΓD = Measure(ΓD, degree)
    dΓN = Measure(ΓN, degree)


  
    # --- Weak formulation --- #  
    ν = 1.0 # positive constant (don't know what it is describing )
    n = 1.0 # volumetric porosity
    

    ϵ(u) = (∇(u) + transpose(∇(u))) 

    ã(v, w) = ∫( n*ν*ϵ(v) ⋅ ∇(w) )dΩf + ∫( ρf*g*∇(ψ) ⋅ Κ * ∇(ϕ) )dΩp + ∫( n*ρf*g*ϕ* w⋅n̂f )dΓ - ∫( n*ρf*g*ψ* v⋅n̂f )dΓ
    b̃(w, q) = -∫( n*q* ∇⋅w )dΩf
    # F̃ = ∫( n* g⋅w )dΩf + ∫( n* l⋅w )dΛf - ∫( ρf*g* ∇ψ ⋅ Κ *∇(Ep)*ϕp )dΩp - ∫( n*ρf*g*Ep*ϕp* w ⋅ n̂f )dΓ

    LHS((u, p), (v, q)) = ã(u, v) + b̃(v, p) + b̃(u, q)
    RHS((v, q)) = ∫( n* g⋅v )dΩf + ∫( n* l⋅v )dΛf - ∫( ρf*g* ∇ψ ⋅ Κ *∇(Ep)*ϕp )dΩp - ∫( n*ρf*g*Ep*ϕp* v ⋅ n̂f )dΓ



     # --- Solve --- #
     op = AffineFEOperator(LHS, RHS, Multi_Trial, Multi_Test)
     uh, ph, ... = solve(op)


end




model, pgs_dict = create_coupled_box(2, true)





