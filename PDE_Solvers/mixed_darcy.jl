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




function mixed_darcy_solver(model, pgs_dict, f0, g0; write = false)

    # --- Boundary tags --- # (Corners?)
    ΓD_tags = pgs_tags(pgs_dict, [1, 2])#5, 6, 7, 8])
    ΓN_tags = pgs_tags(pgs_dict, [3, 4])
    u0 = g0[1]; p0 = g0[2]

    
    # --- Triangulation and spaces --- #
    Ω = Triangulation(model)
    ΓN = BoundaryTriangulation(model, tags=ΓN_tags)
    
    order = 0
    reffeᵤ = ReferenceFE(raviart_thomas, Float64, order) # the 0th order actually referes to raviart_thomas 1st order
    reffeₚ = ReferenceFE(lagrangian, Float64, 0)

    δV = TestFESpace(model, reffeᵤ, conformity=:Hdiv, dirichlet_tags=ΓD_tags)
    δQ = TestFESpace(model, reffeₚ, conformity=:L2)

    V = TrialFESpace(δV, u0) # effectively consider u ⋅ n̂ 
    Q = TrialFESpace(δQ)
    
    δW = MultiFieldFESpace([δV, δQ])
    W = MultiFieldFESpace([V, Q])

    degree = 2 * order
    dΩ = Measure(Ω, degree)
    dΓN = Measure(ΓN, degree) 
    n̂N = get_normal_vector(ΓN)



    # --- Weak formulation --- #
    κ⁻¹ = 1/κ
 
    a((u,p), (v,q)) = ∫( κ⁻¹ * (u ⋅ v) - p * (∇⋅v) - (∇⋅u) * q) * dΩ
    b((v,q)) = ∫(p0 * (-v ⋅ n̂N)) * dΓN - ∫(f0 * q) * dΩ


    # --- Solve --- #
    op = AffineFEOperator(a,b, W, δW)
    uh, ph = solve(op)


    # --- Check with manufactured solution --- #
    u_error = u0 - uh
    p_error = p0 - ph
    
    u_l2norm = sqrt(sum(∫(u_error ⋅ u_error) * dΩ))
    p_l2norm = sqrt(sum(∫(p_error ⋅ p_error) * dΩ))

    @printf("u: l2 norm = %e \n", u_l2norm)
    @printf("p: l2 norm = %e \n", p_l2norm)


    write && writevtk(Ω, path * "vtu_files/" * "mixed_darcy_results", cellfields=["uh" => uh, "ph" => ph, "u_error" => u_error, "p_error" => p_error])


    return u_l2norm, p_l2norm

end



function error_conv(solver, f0, g0; lc_start=2, num_points=5, show_plot=false)
    l2norm = zeros(Float64, num_points, 2)
    lc = zeros(Float64, num_points)
    slope = zeros(Float64, 2)

    # --- Decrease mesh size and get l2norm --- #
    write = false
    for p in 1:num_points
        lc[p] = lc_start * (1 / 2)^(p - 1)
        p == num_points && (write = true)
        model, pgs_dict = create_unit_box(lc[p])
        l2norm[p, :] .= solver(model, pgs_dict, f0, g0; write)

        # Get the actual mesh size
        Ω = Triangulation(model)
        lc[p] = minimum(sqrt.(collect(get_cell_measure(Ω))))
    end

    println("mesh size: ", lc)
    println("u: l2 norm = ", l2norm[:, 1])
    println("p: l2 norm =  ", l2norm[:, 2])


    # ---  Convergence rate --- #
    X = ones(num_points, 2)
    X[:, 2] = log.(lc)
    for i in 1:2
        y = log.(l2norm[:, i])
        b = inv(transpose(X) * X) * transpose(X) * y
        slope[i] = b[2]

    end

    println("u: Convergence rate: ", slope[1])
    println("p: Convergence rate: ", slope[2])

    if show_plot
        fig = plot(lc, l2norm[:, 1], xaxis=:log, yaxis=:log, label="Velocity", marker=:x)
        fig = plot(lc, l2norm[:, 2], xaxis=:log, yaxis=:log, label="Pressure", marker=:x)
        display(fig)
    end
end



κ = 1

# # MS 1
# p0(x) = cos(π * x[1] * x[2])
# u0(x) = VectorValue(κ * π * x[2] * sin(π * x[1] * x[2]), κ * π * x[1] * sin(π * x[1] * x[2]) )
# f0(x) = κ * π^2 * (x[1]^2 + x[2]^2) * cos(π * x[1] * x[2])

# MS 2
p0(x) = cos(π * (x[1] + x[2]))
u0(x) = VectorValue(κ*π * sin(π * (x[1] + x[2])), κ*π * sin(π * (x[1] + x[2])))
f0(x) = κ * π^2 * 2 * cos(π * (x[1] + x[2]))



g0 = (u0, p0)
model, pgs_dict = create_unit_box(0.1, false)
# mixed_darcy_solver(model, pgs_dict, f0, g0; write = true)
error_conv(mixed_darcy_solver, f0, g0)






