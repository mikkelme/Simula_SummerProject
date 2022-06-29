using Gridap
using GridapGmsh

# Equation to solve
# Find scalar field u such that 
# -Δu = -∇²u = f, f ∈ Ω, 
# u = g, on boundary
# ∇u⋅n = h, 
# u g 

model = GmshDiscreteModel("dummy.msh")
# writevtk(model, "model")

order = 1
reffe = ReferenceFE(lagrangian, Float64, order) # type of FE interpolation
V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags = [1]) # Test space

g(x) = 3.0 # Boundary func
U = TrialFESpace(V, g) # Trial space

# Integration 
degree = 2
Ω = Triangulation(model) # Integration mesh of the domain Omega
dΩ = Measure(Ω, degree) # Gauss-like quadrature (how to connect point in numerical integration)

# Neumann conditions
# Γ = BoundaryTriangulation(model, tags = 1)
# dΓ = Measure(Γ, degree)

f(x) = 1.0
h(x) = 2.0

# --- Weak form --- #
a(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ
b(v) = ∫(v * f) * dΩ #+ ∫(v*h)*dΓ

LUSolver(A, b)
# --- Solve --- #

op = AffineFEOperator(a, b, U, V)
ls = LUSolver()
solver = LinearFESolver(ls)
uh = solve(solver, op)
writevtk(Ω, "results", cellfields=["uh" => uh])