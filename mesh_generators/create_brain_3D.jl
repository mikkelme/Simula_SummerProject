include("./brain_mesh_utils.jl")

mutable struct geo3D
    origo::Int32                      # Center of coordinate system
    vertex::Array{Int32,2}            # Vertices of all corners
    arc::Array{Int32,2}               # Tangential arcs
    tan_surf::Array{Int32,1}          # Tangential surfaces
    rad_surf::Array{Int32,2}          # Radial surfaces
    angle::Tuple{Float64,Float64}     # Angle span (z → x-axis, z → y-axis)

    function geo3D() # Constructor
        new(gmsh.model.occ.addPoint(0.0, 0.0, 0.0),
            Matrix{Int32}(undef, 3, 4),
            Matrix{Int32}(undef, 3, 4),
            Vector{Int32}(undef, 3),
            Matrix{Int32}(undef, 2, 4),
            (0, 0))
    end
end



function sphere_patch_corners(r, angle, pert_func=f(x, y) = 0)
    # Calculate corners of sphere patch, 
    # using interception of inclined circles: 
    # x^2 + (z + α_y * y)^2 = y^2 + (z + α_x * x)^2 = r^2

    yx_angle, yz_angle = angle

    x, y_x, _ = spherical_to_cartesian(r, 0, yx_angle / 2)
    _, y_z, z = spherical_to_cartesian(r, pi / 2, yz_angle / 2)



    alpha_x = (r - y_x) / x
    alpha_z = (r - y_z) / z

    a = 2 * alpha_x / (1 - alpha_x^2)
    b = 2 * alpha_z / (1 - alpha_z^2)
    c = r / sqrt(a^2 + b^2 + 1)

    
    A = [a * c, c, b * c]
    B = [-a * c, c, b * c]
    C = [-a * c, c, -b * c]
    D = [a * c, c, -b * c]
    
    NonPerturbed = [A, B, C, D]

    # Add perturbation
    vertex = Array{Int32,1}(undef, 4)
    for i in 1:4
        P = NonPerturbed[i] * (r + pert_func(cartesian_to_surface_cord(NonPerturbed[i]...)...)) / vecNorm(NonPerturbed[i])
        vertex[i] = gmsh.model.occ.addPoint(P...)
    end


    return vertex, NonPerturbed
end



function create_surface(brain::geo3D, r)
    vertex, _ = sphere_patch_corners(r, brain.angle) # Get sphere patch corners
    arc = [gmsh.model.occ.addCircleArc(vertex[i], brain.origo, vertex[i%4 + 1]) for i in 1:4]

    CurveLoop = gmsh.model.occ.addCurveLoop([arc...])
    surf = gmsh.model.occ.addSurfaceFilling(CurveLoop)

    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [surf])
    return vertex, arc, surf

end


function cartesian_to_surface_cord(x, y, z)
    # Go from 3D cartesian coordinates to 
    # 2D cartesian surface coordinates. 
    # Direction to travel along circle arcs in 
    # x and y direction respectively

    # Is all calculations safe from trigonometric traps?

    r = sqrt(x^2 + y^2 + z^2)

    c = y
    a = x / c
    b = z / c


    alpha_x = -(1 - sqrt(1 + a^2)) / a
    alpha_z = -(1 - sqrt(1 + b^2)) / b

    angle_x = 2 * atan(alpha_x)
    angle_z = 2 * atan(alpha_z)

    x_arcLen = r * angle_x
    z_arcLen = r * angle_z

   return x_arcLen, z_arcLen

end


function add_perturbed_arc(start_tag, center_tag, end_tag, NonPerturbed_start, NonPerturbed_end, pert_func, BS_points)
    pointTags = [start_tag]

    start_point = NonPerturbed_start
    end_point = NonPerturbed_end
    origo = gmsh.model.getValue(0, center_tag, [])
    r = vecNorm(start_point)


    # Draw arc
    direction = end_point - start_point

    # Calculate spline points in given direction
    BS_points = Int(sum(abs.(direction / vecNorm(direction))[[1, 3]] .* BS_points))
  

    range = LinRange(0, 1, BS_points)
    for i in 2:BS_points-1
        # Step toward end_point in straight line
        point = start_point + range[i] * direction

        # Normalize to be on sphere 
        point *= r / vecNorm(point)

        # Add perturbation to radius
        perturbed_radius = r + pert_func(cartesian_to_surface_cord(point...)...)

        # Normalize to be on perturbed sphere
        point *= perturbed_radius / vecNorm(point)
        append!(pointTags, gmsh.model.occ.addPoint(point...))
    end

    append!(pointTags, end_tag)
    BSpline = gmsh.model.occ.addBSpline(pointTags)
    return BSpline

end


function create_perturbed_surface(brain::geo3D, r, pert_func, BS_points)

    vertex, NonPerturbed = sphere_patch_corners(r, brain.angle, pert_func) # Get sphere patch corners
    arc = [add_perturbed_arc(vertex[i], brain.origo, vertex[i%4+1], NonPerturbed[i], NonPerturbed[i%4+1], pert_func, BS_points) for i in 1:4]

    CurveLoop = gmsh.model.occ.addCurveLoop([arc...])
    surf = gmsh.model.occ.addBSplineFilling(CurveLoop, -1, "Stretch") # Alternatively use "Coons"
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [surf])

    return vertex, arc, surf

end


function connect_and_volumize(brain::geo3D)
    # Connect surfaces and create volumes
    for i in 1:2
        vline = [gmsh.model.occ.addLine(brain.vertex[i, 1], brain.vertex[i+1, 1])] # Vertical lines
        Loop = [] # curve loops
        for j in 2:4
            append!(vline, gmsh.model.occ.addLine(brain.vertex[i, j], brain.vertex[i+1, j]))
            append!(Loop, gmsh.model.occ.addCurveLoop([vline[j-1], brain.arc[i+1, j-1], -vline[j], -brain.arc[i, j-1]]))
        end
        append!(Loop, gmsh.model.occ.addCurveLoop([vline[4], brain.arc[i+1, 4], -vline[1], -brain.arc[i, 4]]))
    
        brain.rad_surf[i, :] = [gmsh.model.occ.addSurfaceFilling(l) for l in Loop]
    

        surfLoop = gmsh.model.occ.addSurfaceLoop([brain.tan_surf[i], brain.rad_surf[i,:]..., brain.tan_surf[i+1]])  
        vol = gmsh.model.occ.addVolume([surfLoop])
        gmsh.model.occ.synchronize()

        # Add physical gorups
        [gmsh.model.addPhysicalGroup(2, [s]) for s in brain.rad_surf[i]] # surface
        gmsh.model.addPhysicalGroup(3, [vol])  # volume
    
    end
end


function add_mesh_field(brain::geo3D, param::model_params)
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
    gmsh.model.mesh.field.setNumbers(F_distance, "SurfacesList", [brain.tan_surf[2]])



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


function apply_periodic_meshing(brain::geo3D)
    yx_angle, yz_angle = brain.angle

    

    # Rotation around x-axis (do not understand why the minues is needed)
    Rx = [1, 0, 0, 0,   
        0, cos(-yz_angle), sin(-yz_angle), 0,
        0, -sin(-yz_angle), cos(-yz_angle), 0,
        0, 0, 0, 1]


    # Rotation around z-axis 
    Rz = [  cos(yx_angle), -sin(yx_angle), 0, 0,
            sin(yx_angle), cos(yx_angle), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1 ]

 
    gmsh.model.mesh.setPeriodic(2, [brain.rad_surf[1, 1]], [brain.rad_surf[1, 3]], Rx)
    gmsh.model.mesh.setPeriodic(2, [brain.rad_surf[1, 2]], [brain.rad_surf[1, 4]], Rz)


    for i in 1:2
        gmsh.model.mesh.setPeriodic(2, [brain.rad_surf[i, 1]], [brain.rad_surf[i, 3]], Rx) # From -z → +z surf
        gmsh.model.mesh.setPeriodic(2, [brain.rad_surf[i, 2]], [brain.rad_surf[i, 4]], Rz) # From +x → -x surf        
    end
end


function create_brain_3D(param::model_params)
    #?--> Safety check of parameters? # Should enforce  0 < angle < pi

    brain = geo3D() # Struct for holding tags and angle

    # Calculate derived parameters
    rI = param.r_curv - param.r_brain                   # Inner radius  
    rD = param.r_curv - param.d_ratio * param.r_brain   # Radius for dividing line
    x_arcLen, z_arcLen = param.arcLen
    yx_angle = x_arcLen / param.r_curv                  # Angle span y -> x -axis
    yz_angle = z_arcLen / param.r_curv                  # Angle span y -> z -axis


    brain.angle = (yx_angle, yz_angle)

    # Add radial surfaces
    brain.vertex[1, :], brain.arc[1, :], brain.tan_surf[1] = create_surface(brain, rI)
    brain.vertex[2, :], brain.arc[2, :], brain.tan_surf[2] = create_perturbed_surface(brain, rD, param.inner_perturb, param.BS_points)
    brain.vertex[3, :], brain.arc[3, :], brain.tan_surf[3] = create_perturbed_surface(brain, param.r_curv, param.outer_perturb, param.BS_points)


    connect_and_volumize(brain)
    add_mesh_field(brain, param)
    apply_periodic_meshing(brain)

    # To do: set physical groups in a more organized way here (in a loop: se 2D case)


    # Generate mesh
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)



end # End of create_brain_3D



