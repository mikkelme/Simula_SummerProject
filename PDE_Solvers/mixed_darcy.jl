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




function mixed_darcy_solver(model, pgs_dict; write = false)


    ΓD_tags = pgs_tags(pgs_dict, [1, 3])

    order = 1
    reffeᵤ = ReferenceFE(raviart_thomas, Float64, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order)


    # FESpace = TestFESpace <---- ???
    δV = FESpace(model, reffeᵤ, conformity=:HDiv, dirichlet_tags=ΓD_tags)
    δQ = FESpace(model, reffeₚ, conformity=:L2)

    uD = VectorValue(0.0,0.0)
    V = TrialFESpace(δV, uD)
    Q = TrialFESpace(δQ)
    
    δW = MultiFieldFESpace([δV, δQ])
    W = MultiFieldFESpace([V, Q])

    degree = 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    ΓN_tags = pgs_tags(pgs_dict, [2])
    ΓN = BoundaryTriangulation(model, tags=ΓN_tags)
    dΓN = Measure(ΓN, degree) 



    # --- Weak formulation --- #
    κinv₁ = TensorValue(1.0, 0.0, 0.0, 1.0)
    κinv₂ = TensorValue(100.0, 90.0, 90.0, 100.0)

    function σ(x,u)
        if ((abs(x[1]-0.5) <= 0.1) && (abs(x[2]-0.5) <= 0.1))
           return κinv₁⋅u
        else
           return κinv₂⋅u
        end
     end


    px = get_physical_coordinate(Ω)

    a((u,p), (v,q)) = ∫(v⋅(σ ∘ (px,u)) - (∇⋅v) * p + q * (∇⋅u)) * dΩ


    nb = get_normal_vector(ΓN)
    h = -1.0

    b((v,q)) = ∫((v ⋅ nb) * h) * dΓN


    op = AffineFEOperator(a,b, W, δW)
    uh, ph = solve(op)


    write && writevtk(Ω, path * "vtu_files/" * "darcy_results", cellfields=["uh"=>uh,"ph"=>ph])

end




model, pgs_dict = create_unit_box(0.1, false)
mixed_darcy_solver(model, pgs_dict; write = true)




a((u,p), (v,q)) = ∫( κ⁻¹⋅u⋅v - p * (∇⋅v) - (∇⋅u) * q) * dΩ
b((v,q)) = -∫(gₚ * (v ⋅ n̂N)) * dΓN - ∫(f0 ⋅ v) * dΩ