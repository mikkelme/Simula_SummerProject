using Gridap
using GridapGmsh


function stokes(b_vals)
    model = GmshDiscreteModel(meshfile)
    b_tags = collect(keys(b_vals))

    # Define reference FE (Q2/P1(disc) pair)
    order = 2
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order - 1; space=:P)


    # Define test FESpaces
    V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=b_tags)
    Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)
    Y = MultiFieldFESpace([V, Q])

    # Define trial FESpaces from Dirichlet values
    b_conditions = [VectorValue(b_vals[tag]...) for tag in b_tags]
    U = TrialFESpace(V, b_conditions)
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U, P])

    # Define triangulation and integration measure
    degree = order
    Ωₕ = Triangulation(model)
    dΩ = Measure(Ωₕ, degree)

    # Define bilinear and linear form
    f = VectorValue(0.0, 0.0)
    # f(x) = x
    a((u, p), (v, q)) = ∫(∇(v) ⊙ ∇(u) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ
    l((v, q)) = ∫(v ⋅ f)dΩ

    # Build affine FE operator
    op = AffineFEOperator(a, l, X, Y)

    # Solve
    uh, ph = solve(op)

    # Export results to vtk
    writevtk(Ωₕ, "general_results", order=2, cellfields=["uh" => uh, "ph" => ph])
end


meshfile = "dummy.msh"
b_vals = Dict(1 => [0, 0], 2 => [0, 0], 3 => [10, 0], 4 => [0, 0]) # boundary values


stokes(b_vals)

# Diri = Dict("a" => 1, "b" => 2, "c" => 3) 
# Diri["d"] = 20
# println(Diri[Dtags])



#-------> Rad about: Method of Manufactured Solutions
#https://www.comsol.com/blogs/verify-simulations-with-the-method-of-manufactured-solutions/