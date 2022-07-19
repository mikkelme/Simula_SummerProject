
# --- Custom Function --- # 

function direct_wiring(gmsh; renumber=true)
    renumber && gmsh.model.mesh.renumberNodes()
    renumber && gmsh.model.mesh.renumberElements()

    Dc = GridapGmsh._setup_cell_dim(gmsh)
    Dp = GridapGmsh._setup_point_dim(gmsh, Dc)
    node_to_coords = GridapGmsh._setup_node_coords(gmsh, Dp)
    nnodes = length(node_to_coords)
    vertex_to_node, node_to_vertex = _setup_nodes_and_vertices(gmsh, node_to_coords)
    grid, cell_to_entity = _setup_grid(gmsh, Dc, Dp, node_to_coords, node_to_vertex)
    cell_to_vertices = _setup_cell_to_vertices(Gridap.Geometry.get_cell_node_ids(grid), node_to_vertex, nnodes)
    grid_topology = Gridap.Geometry.UnstructuredGridTopology(grid, cell_to_vertices, vertex_to_node)
    labeling = GridapGmsh._setup_labeling(gmsh, grid, grid_topology, cell_to_entity, vertex_to_node, node_to_vertex)

    pgs = gmsh.model.getPhysicalGroups()
 
    model = Gridap.Geometry.UnstructuredDiscreteModel(grid, grid_topology, labeling)
    pgs_dict = Dict(Int64(pgs[i][2]) => Int64(i) for i in 1:length(pgs)) # Physical group dictionary
  
    return model, pgs_dict
end


# --- SUPPORT FUNCTIONS: COPY-PASTED FROM https://github.com/gridap/GridapGmsh.jl/blob/master/src/GmshDiscreteModels.jl --- #


function  _setup_nodes_and_vertices(gmsh,node_to_coords)
  nnodes = length(node_to_coords)
  dimTags = gmsh.model.getEntities()
  if GridapGmsh._has_periodic_bcs(gmsh,dimTags)
    dimTags = gmsh.model.getEntities()
    vertex_to_node, node_to_vertex = _setup_nodes_and_vertices_periodic(gmsh,dimTags,nnodes)
  else
    vertex_to_node = 1:nnodes
    node_to_vertex = vertex_to_node
  end
  vertex_to_node, node_to_vertex
end


function _setup_nodes_and_vertices_periodic(gmsh,dimTags,nnodes)
    # Assumes linear grid
    node_to_node_master = fill(GridapGmsh.UNSET,nnodes)
    GridapGmsh._node_to_node_master!(node_to_node_master,gmsh,dimTags)
    slave_to_node_slave = findall(node_to_node_master .!= GridapGmsh.UNSET)
    slave_to_node_master = node_to_node_master[slave_to_node_slave]
    node_to_vertex = fill(GridapGmsh.UNSET,nnodes)
    vertex_to_node = findall(node_to_node_master .== GridapGmsh.UNSET)
    node_to_vertex[vertex_to_node] = 1:length(vertex_to_node)
    nmax = 20
    for i in 1:nmax
      node_to_vertex[slave_to_node_slave] = node_to_vertex[slave_to_node_master]
      if all(j->j!=0,node_to_vertex)
        break
      end
      if i == nmax
        Gridap.Helpers.@unreachable
      end
    end
    vertex_to_node, node_to_vertex
  end


function _setup_grid(gmsh, Dc, Dp, node_to_coords, node_to_vertex)

    if Dp == 3 && Dc == 2
        orient_if_simplex = false
    else
        orient_if_simplex = true
    end

    cell_to_nodes, nminD = _setup_connectivity(gmsh, Dc, node_to_vertex, orient_if_simplex)
    cell_to_type, reffes, orientation = GridapGmsh._setup_reffes(gmsh, Dc, orient_if_simplex)
    cell_to_entity = GridapGmsh._setup_cell_to_entity(
        gmsh, Dc, length(cell_to_nodes), nminD)

    if Dp == 3 && Dc == 2
        cell_coords = lazy_map(Broadcasting(Reindex(node_to_coords)), cell_to_nodes)
        ctype_shapefuns = map(Gridap.ReferenceFEs.get_shapefuns, reffes)
        cell_shapefuns = Gridap.ReferenceFEs.expand_cell_data(ctype_shapefuns, cell_to_type)
        cell_map = lazy_map(Gridap.ReferenceFEs.linear_combination, cell_coords, cell_shapefuns)
        ctype_x = fill(zero(VectorValue{Dc,Float64}), length(ctype_shapefuns))
        cell_x = Gridap.ReferenceFEs.expand_cell_data(ctype_x, cell_to_type)
        cell_Jt = lazy_map(∇, cell_map)
        cell_n = lazy_map(Operation(GridapGmsh._unit_outward_normal), cell_Jt)
        cell_nx = lazy_map(evaluate, cell_n, cell_x) |> collect
        facet_normal = lazy_map(Gridap.ReferenceFEs.constant_field, cell_nx)
    else
        facet_normal = nothing
    end

    grid = Gridap.Geometry.UnstructuredGrid(
        node_to_coords,
        cell_to_nodes,
        reffes,
        cell_to_type,
        orientation,
        facet_normal)

    (grid, cell_to_entity)

end


function _setup_connectivity(gmsh, d, node_to_vertex, orient_if_simplex)

    elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(d)

    if length(elemTypes) == 0
        ncells = 0
        ndata = 0
        nmin = 1
        cell_to_nodes_prts = zeros(Int, ncells + 1)
        cell_to_nodes_data = zeros(Int32, ndata)
        cell_to_nodes = Table(cell_to_nodes_data, cell_to_nodes_prts)
        return (cell_to_nodes, nmin)
    end

    ncells, nmin, nmax = GridapGmsh._check_cell_tags(elemTags)

    etype_to_nlnodes = GridapGmsh._setup_etype_to_nlnodes(elemTypes, gmsh)

    ndata = sum([length(t) for t in nodeTags])

    cell_to_nodes_data = zeros(Int, ndata)
    cell_to_nodes_prts = zeros(Int32, ncells + 1)

    _fill_connectivity!(
        cell_to_nodes_data,
        cell_to_nodes_prts,
        nmin - 1,
        etype_to_nlnodes,
        elemTypes,
        elemTags,
        nodeTags,
        d,
        node_to_vertex,
        orient_if_simplex)

    cell_to_nodes = Gridap.Arrays.Table(cell_to_nodes_data, cell_to_nodes_prts)

    (cell_to_nodes, nmin)

end


function _fill_connectivity!(
    cell_to_nodes_data,
    cell_to_nodes_prts,
    o,
    etype_to_nlnodes,
    elemTypes,
    elemTags,
    nodeTags,
    d,
    node_to_vertex,
    orient_if_simplex)

    for (j, etype) in enumerate(elemTypes)
        nlnodes = etype_to_nlnodes[etype]
        i_to_cell = elemTags[j]
        for cell in i_to_cell
            cell_to_nodes_prts[cell+1-o] = nlnodes
        end
    end

    Gridap.Geometry.length_to_ptrs!(cell_to_nodes_prts)

    for (j, etype) in enumerate(elemTypes)
        nlnodes = etype_to_nlnodes[etype]
        i_to_cell = elemTags[j]
        i_lnode_to_node = nodeTags[j]
        if (nlnodes == d + 1) && orient_if_simplex
            # what we do here has to match with the OrientationStyle we
            # use when building the UnstructuredGrid
            _orient_simplex_connectivities!(nlnodes, i_lnode_to_node, node_to_vertex)
        elseif (nlnodes == 4)
            GridapGmsh._sort_quad_connectivites!(nlnodes, i_lnode_to_node)
        elseif (nlnodes == 8)
            _sort_hex_connectivites!(nlnodes, i_lnode_to_node)
        end
        for (i, cell) in enumerate(i_to_cell)
            a = cell_to_nodes_prts[cell-o] - 1
            for lnode in 1:nlnodes
                node = i_lnode_to_node[(i-1)*nlnodes+lnode]
                cell_to_nodes_data[a+lnode] = node
            end
        end
    end

end


function _setup_cell_to_vertices(cell_to_nodes, node_to_vertex, nnodes)
    if isa(node_to_vertex, AbstractVector)
        cell_to_vertices = Gridap.Arrays.Table(lazy_map(Broadcasting(Reindex(node_to_vertex)), cell_to_nodes))
    else
        @assert node_to_vertex == 1:nnodes
        cell_to_vertices = cell_to_nodes
    end
    cell_to_vertices
end


function _orient_simplex_connectivities!(nlnodes, i_lnode_to_node, node_to_vertex)
    aux = zeros(eltype(i_lnode_to_node), nlnodes)
    offset = nlnodes - 1
    for i in 1:nlnodes:length(i_lnode_to_node)
        nodes = i_lnode_to_node[i:i+offset]
        vertices = view(node_to_vertex, nodes)
        perm = sortperm(vertices)
        i_lnode_to_node[i:i+offset] = nodes[perm]
    end
end