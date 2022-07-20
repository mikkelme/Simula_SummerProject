include("./brain_mesh_utils.jl")

mutable struct geo2D
    origo::Int32             # Center of coordinate system
    vertex::Array{Int32,2}   # Vertices of all corners
    arc::Array{Int32,1}      # arcs
    vline::Array{Int32,2}    # Vertical lines
    surf::Array{Int32,1}     # Surfaces
    angle::Float64           # Angle span (y → x-axis)

    function geo2D() # Constructor
        new(gmsh.model.occ.addPoint(0.0, 0.0, 0.0),
            Matrix{Int32}(undef, 3, 2),
            Vector{Int32}(undef, 3),
            Matrix{Int32}(undef, 2, 2),
            Vector{Int32}(undef, 2),
             0)
    end
end




function create_arc(brain::geo2D, r)
    vertex = Array{Int32,1}(undef, 2)
    vertex[1] = gmsh.model.occ.addPoint(spherical_to_cartesian(r, pi, brain.angle / 2)...)
    vertex[2] = gmsh.model.occ.addPoint(spherical_to_cartesian(r, 0, brain.angle / 2)...)
    arc = gmsh.model.occ.addCircleArc(vertex[1], brain.origo, vertex[2])
    gmsh.model.occ.synchronize()
    return vertex, arc
end

function create_perturbed_arc(brain::geo2D, r, perturbation_func, BS_points)
    pointTags = []

    for i in LinRange(0, 1, BS_points)
        phi = - brain.angle / 2 + brain.angle * i
        x_arcLen = r * phi
        r_pert = r + perturbation_func(x_arcLen, 0)
        append!(pointTags, gmsh.model.occ.addPoint(spherical_to_cartesian(r_pert, 0, phi)...))
    end

    vertex = [pointTags[1], last(pointTags)]
    BSpline = gmsh.model.occ.addBSpline(pointTags)
    gmsh.model.occ.synchronize()
    return vertex, BSpline
end

function connect_and_surfize(brain::geo2D)
    for i in 1:2
        brain.vline[i,:] = [gmsh.model.occ.addLine(brain.vertex[i, j], brain.vertex[i+1, j]) for j in 1:2]
        Loop = gmsh.model.occ.addCurveLoop([brain.vline[i,1], brain.arc[i+1], -brain.vline[i,2], -brain.arc[i]])
        brain.surf[i] = gmsh.model.occ.addPlaneSurface([Loop])
        gmsh.model.occ.synchronize()
    end
end

function add_mesh_field(brain::geo2D, param::model_params)
    # Linearly change mesh size from LcMin → LcMax ∈ [DistMin, DistMax]

    # LcMin -                       /------------------
    #                              /
    #                             /
    #                            /
    # LcMax   -o----------------/
    #          |                |    |
    #       Surface        DistMin  DistMax

    # Get parameters
    LcMin, LcMax = param.lc * param.field_Lc_lim
    DistMin, DistMax = param.field_Dist_lim

    # Define distance field
    F_distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumber(F_distance, "NNodesByEdge", 100)
    gmsh.model.mesh.field.setNumbers(F_distance, "EdgesList", [brain.arc[2]])

    

    # Define trheshold field
    F_threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(F_threshold, "IField", F_distance)
    gmsh.model.mesh.field.setNumber(F_threshold, "LcMin", LcMin)
    gmsh.model.mesh.field.setNumber(F_threshold, "LcMax", LcMax)
    gmsh.model.mesh.field.setNumber(F_threshold, "DistMin", DistMin)
    gmsh.model.mesh.field.setNumber(F_threshold, "DistMax", DistMax)
    gmsh.model.mesh.field.setAsBackgroundMesh(F_threshold)

    # Prevent over-refinement due to small mesh sizes on the boundary.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.occ.synchronize()
end


function apply_periodic_meshing(brain::geo2D)
    ϕ = brain.angle

    Rz = [  cos(ϕ), -sin(ϕ), 0, 0,
            sin(ϕ), cos(ϕ), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1 ]

    # Copy mesh from right (brain.vline[i, 2]) to left (brain.vline[i, 1])  
    println(brain.vline[1, 1], " ", brain.vline[1, 2])
    gmsh.model.mesh.setPeriodic(1, [brain.vline[1, 1]], [brain.vline[1, 2]], Rz)  # Inner part
    gmsh.model.mesh.setPeriodic(1, [brain.vline[2, 1]], [brain.vline[2, 2]], Rz)  # Outer part
    # dim, tags, tagsMaster, affine
    # rotate tags tagsMaster to tags to copy-paste the mesh onto tags
end


function create_brain_2D_xy(param::model_params)
    brain = geo2D() # Struct for holding tags and angle

    # Calculate derived parameters
    rI = param.r_curv - param.r_brain                   # Inner radius  
    rD = param.r_curv - param.d_ratio * param.r_brain   # Radius for dividing line
    arcLen = param.arcLen[1]                            # Length of outer arc
    brain.angle = arcLen / param.r_curv                 # Angle span y -> x -axis
    
    
    # Add arcs
    brain.vertex[1, :], brain.arc[1] = create_arc(brain, rI)
    brain.vertex[2, :], brain.arc[2] = create_perturbed_arc(brain, rD, param.inner_perturb, param.BS_points[1])
    brain.vertex[3, :], brain.arc[3] = create_perturbed_arc(brain, param.r_curv, param.outer_perturb, param.BS_points[1])
    
    # Connect and refine mesh
    connect_and_surfize(brain)
    add_mesh_field(brain, param)
    apply_periodic_meshing(brain)
    
    
    # Add physical groups 
    obj = [(brain.surf, 2), (brain.arc, 1), (brain.vline, 1), (brain.vertex, 0)] # (tag, dim)
    for k in 1:length(obj)
        dim = obj[k][2]
        for tag in obj[k][1]
            gmsh.model.addPhysicalGroup(dim, [tag])
        end
    end
    

    # Generate mesh
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)


end




