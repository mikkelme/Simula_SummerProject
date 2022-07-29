using LightGraphs
using Plots 


"""
Let V be a TestFESpace on a simple open curve. We compute a function 
with values of archlength coordinate.
"""
function arclength_coordinate(V; tol=1E-10)
    # FIXME:
    # Check that we are in 2d and have some scalar elements with dofs being 
    # the point evaluations
    Γ = get_triangulation(V)
    # The idea is to build a graph and walk it
    cell_idx = Γ.plus.trian.tface_to_mface
    parent_topology = Γ.plus.trian.model.grid_topology.n_m_to_nface_to_mfaces[2, 1]
    # Encode the cells of \Gamma in terms of parent vertices 
    cell_vtx = parent_topology[cell_idx, :]
    l2g = unique(hcat(cell_vtx...))  # In parent numbering
    # We want to build the graph in terms of local numbering
    g2l = Dict(map(reverse, enumerate(l2g)))
    G = SimpleGraph(length(l2g))
    for (g0, g1) ∈ cell_vtx
        add_edge!(G, g2l[g0], g2l[g1])
    end

    # What's the degree of each vertex?
    degrees = degree(G)
    # Check that this is a simple open curve
    @assert all(1 .<= degrees .<= 2)
    # Our path will march between two vertices connected only to single cell each
    start, stop = findall(isequal(1), degrees)
    path = enumerate_paths(dijkstra_shortest_paths(G, start), stop)
    
    cell_vtx_l = Dict{Tuple{Int, Int}, Int}()
    # Now we want to walk in terms of cells so we build a lookup
    for (ci, (g0, g1)) ∈ enumerate(cell_vtx)
        l0, l1 = g2l[g0], g2l[g1]
        # Sort the key
        if l0 < l1
            cell_vtx_l[(l0, l1)] = ci
        else
            cell_vtx_l[(l1, l0)] = ci
        end
    end

    # While we walk we want to build the arclength
    node_x = Γ.plus.trian.grid.parent.node_coordinates[l2g]
    # The idea is to insert nodes based on their distance; so let's get their position
    # Here we assume 2D
    dofs_x = get_free_dof_values(interpolate_everywhere(x -> x[1], V))
    dofs_y = get_free_dof_values(interpolate_everywhere(x -> x[2], V))
    dm = get_cell_dof_ids(V)
    # The vector of dof values we are building is 
    dist = similar(dofs_x)
    r_cord = similar(dofs_x)
 

    distance = 0
    for i ∈ 1:(length(path)-1)
        l0, l1 = (path[i], path[i+1])
        key = l0 < l1 ? (l0, l1) : (l1, l0)
        # Fing the cells with these two vertices
        ci = cell_vtx_l[key]

        # Start and end coord
        x0, y0 = node_x[l0]
        x1, y1 = node_x[l1]
        edge_length = sqrt((x0 - x1)^2 + (y0 - y1)^2)
        
        cell_dofs = dm[ci]
        # @show x_cord
        # For the dofs to be set get their distance from l0
        dofs_dist = sqrt.((dofs_x[cell_dofs] .- x0).^2 .+ (dofs_y[cell_dofs] .- y0).^2)
        @assert all(-tol .< dofs_dist .< edge_length+tol)
        # The arclength is based on the cumsum
        dist[cell_dofs] = dofs_dist .+ distance
        x_cord = [dofs_x[i] for i in cell_dofs]
        y_cord = [dofs_y[i] for i in cell_dofs]
        r_cord[cell_dofs] = sqrt.(x_cord.^2 .+ y_cord.^2)
        
        # For the next round we start with the l1 vertex
        distance += edge_length
    end

    
    return FEFunction(V, dist), r_cord
end

"""L^2 project expression into V"""
function project(expr, δV; qdegree=4)
    Ω = get_triangulation(δV)

    V = TrialFESpace(δV)
    dΩ = Measure(Ω, qdegree)

    a(u, v) = ∫(u⋅v)*dΩ
    L(v) = ∫(v⋅expr)*dΩ

    A = assemble_matrix(a, δV, δV)
    b = assemble_vector(L, δV)
    x = A\b 

    uh = FEFunction(V, x)
end



function plot_nflow_profile(brain_param, u, Γ, savepath = false, info = 0.0)
        # For representing flux
        Qelm = ReferenceFE(lagrangian, Float64, 1)
        δQ = TestFESpace(Γ, Qelm)

        ν = get_normal_vector(Γ)

        uh = project(u.⁺⋅ν.⁺, δQ)


        # Compute the arclength coordinate
        al, r_cord = arclength_coordinate(δQ)

        # NOTE: now we want to plot stuff against arc length. The dofs are not 
        # ordered in a way that x_val below is monotone...
        x_val = get_free_dof_values(al)
        # ... That's why we need to reorder
        idx = sortperm(x_val)
        y_val = get_free_dof_values(uh)

    
        red_idx = idx[5:length(idx)-5]
        if brain_param.inner_perturb_body != "(x,z) -> 0.0" 
            fig = plot(x_val[red_idx]*1e3, r_cord[red_idx]*1e3, label = "Interface radius", ylabel = "Interface radius [mm]", legend=:topleft, color = :black, alpha = 0.3, left_margin = 5Plots.mm, right_margin = 18Plots.mm)
            plot!(twinx(), x_val[red_idx]*1e3, y_val[red_idx]*1e3, label = "Interface normal velocity", ylabel = "Interface normal velocity [mm/s]", legend=:topright, left_margin = 5Plots.mm, right_margin = 18Plots.mm)
        else 
            fig = plot(x_val[red_idx]*1e3, y_val[red_idx]*1e3, label = "Interface normal velocity", ylabel = "Interface normal velocity [mm/s]", legend=:topleft, left_margin = 5Plots.mm, right_margin = 18Plots.mm)
        end
        
        xlabel!("Curve length coordinate [mm]")
        savepath != false && savefig(fig, savepath * "nflow_profile" * info * ".png")







end




# model_path, _ = split_square_mesh(0.2; offset=0.2, distance=2)

# model = GmshDiscreteModel(model_path)

# Ω0 = Triangulation(model, tags=["top_surface"])
# Ω1 = Triangulation(model, tags=["bottom_surface"])
# Γ = InterfaceTriangulation(Ω0, Ω1)
# # Let's have some vector space
# Velm = ReferenceFE(lagrangian, VectorValue{2, Float64}, 1)
# V = TestFESpace(Ω0, Velm)
# u = interpolate_everywhere(identity, V)

# # For representing flux
# Qelm = ReferenceFE(lagrangian, Float64, 1)
# δQ = TestFESpace(Γ, Qelm)

# ν = get_normal_vector(Γ)

# uh = project(u.⁺⋅ν.⁺, δQ)

# # Compute the arclength coordinate
# al = arclength_coordinate(δQ)


# using Plots 
# # NOTE: now we want to plot stuff against arc length. The dofs are not 
# # ordered in a way that x_val below is monotone...
# x_val = get_free_dof_values(al)
# # ... That's why we need to reorder
# idx = sortperm(x_val)

# y_val = get_free_dof_values(uh)

# plot(x_val[idx], y_val[idx])