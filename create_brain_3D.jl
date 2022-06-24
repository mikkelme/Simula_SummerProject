using GridapGmsh: gmsh
using Printf


struct model_params
    # Parameters defining the brain model
    lc::Float64                       # Mesh size
    arcLen::Tuple{Float64,Float64}    # Outer arc length (x, y)-direction
    r_brain::Float64                  # Radial total length
    d_ratio::Float64                  # Radial relative length of fluid section
    r_curv::Float64                   # Radius of curvature
    inner_perturb::Function           # Radial perturbation on inner surface f(x,y,z)
    outer_perturb::Function           # Radial perturbation on inner surface f(x,y,z)
    BS_points::Int64                  # Number of points in BSpline curves
    field_Lc_lim::Array{Float64,1}    # Field strengh as (LcMin, LcMax) multiplied by Lc when used
    field_Dist_lim::Array{Float64,1}  # field distance (DistMin, DistMax)
end


mutable struct geo3D
    # For holding gmsh tags
    origo::Int32                      # Center of coordinate system
    vertex::Array{Int32,2}            # Vertices of all corners
    arc::Array{Int32,2}               # Tangential arcs
    tan_surf::Array{Int32,1}          # Tangential surfaces
    rad_surf::Array{Int32,2}          # Radial surfaces
    angle::Tuple{Float64,Float64}     # Angle span (z → x-axis, z → x-axis)

    function geo3D() # Constructor
        new(gmsh.model.occ.addPoint(0.0, 0.0, 0.0),
            Matrix{Int32}(undef, 3, 4),
            Matrix{Int32}(undef, 3, 4),
            Vector{Int32}(undef, 3),
            Matrix{Int32}(undef, 2, 4),
            (0, 0))
    end
end


function spherical_to_cartesian(r, theta, phi)
    # theta: angle from x-axis in x-y-plane (0 → 2π)
    # phi: angle from z-axis towards x-y-plane (0 → π)
    x = r * cos(theta) * sin(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(phi)
    return x, y, z
end


function sphere_patch_corners(r, angle, pert_func=f(x, y, z) = 0)
    # Calculate corners of sphere patch, 
    # using interception of inclined circles: 
    # x^2 + (z + α_y * y)^2 = y^2 + (z + α_x * x)^2 = r^2

    zx_angle, zy_angle = angle
    vertex = Array{Int32,1}(undef, 4)



    x, _, z_x = spherical_to_cartesian(r, 0, zx_angle / 2)
    _, y, z_y = spherical_to_cartesian(r, pi / 2, zy_angle / 2)

    alpha_x = (r - z_x) / x
    alpha_y = (r - z_y) / y

    a = 2 * alpha_x / (1 - alpha_x^2)
    b = 2 * alpha_y / (1 - alpha_y^2)
    c = r / sqrt(a^2 + b^2 + 1)

    # Add perturbation
    A = [a * c, b * c, c]
    B = [-a * c, b * c, c]
    C = [-a * c, -b * c, c]
    D = [a * c, -b * c, c]

    A_pert = A * (r + pert_func(A...)) / sqrt(A[1]^2 + A[2]^2 + A[3]^2)
    B_pert = B * (r + pert_func(B...)) / sqrt(B[1]^2 + B[2]^2 + B[3]^2)
    C_pert = C * (r + pert_func(C...)) / sqrt(C[1]^2 + C[2]^2 + C[3]^2)
    D_pert = D * (r + pert_func(D...)) / sqrt(D[1]^2 + D[2]^2 + D[3]^2)

    # Add points
    vertex[1] = gmsh.model.occ.addPoint(A_pert...)
    vertex[2] = gmsh.model.occ.addPoint(B_pert...)
    vertex[3] = gmsh.model.occ.addPoint(C_pert...)
    vertex[4] = gmsh.model.occ.addPoint(D_pert...)

    #-----> Working here perturbating the cornes, not sure if it works completely yet
    # After doing that go back and finish implementation of boundary conditions

    return vertex
end



function create_surface(brain::geo3D, r)
    vertex = sphere_patch_corners(r, brain.angle) # Get sphere patch corners
    arc = Array{Int32,1}(undef, 4)

    arc[1] = gmsh.model.occ.addCircleArc(vertex[1], brain.origo, vertex[2])
    arc[2] = gmsh.model.occ.addCircleArc(vertex[2], brain.origo, vertex[3])
    arc[3] = gmsh.model.occ.addCircleArc(vertex[3], brain.origo, vertex[4])
    arc[4] = gmsh.model.occ.addCircleArc(vertex[4], brain.origo, vertex[1])

    CurveLoop = gmsh.model.occ.addCurveLoop([arc[1], arc[2], arc[3], arc[4]])
    surf = gmsh.model.occ.addSurfaceFilling(CurveLoop)

    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [surf])
    return vertex, arc, surf

end


function add_perturbed_arc(start_tag, center_tag, end_tag, pert_func, BS_points)
    pointTags = [start_tag]

    # Convert: tag -> coordinates 
    start_point = gmsh.model.getValue(0, start_tag, []) #dim, tag, parametrization
    end_point = gmsh.model.getValue(0, end_tag, [])
    origo = gmsh.model.getValue(0, center_tag, [])
    r = sqrt(start_point[1]^2 + start_point[2]^2 + start_point[3]^2)

    # Draw arc
    direction = end_point - start_point
    range = LinRange(0, 1, BS_points)
    for i in 2:BS_points-1
        # Step toward end_point in straight line
        point = start_point + range[i] * direction
    
        # Normalize to be on normal sphere arc
        point *= r / sqrt((point[1] - origo[1])^2 + (point[2] - origo[2])^2 + (point[3] - origo[3])^2)
    
        # Add perturbation to radius
        # perturbed_radius = r + pert_func((point - start_point)...)
        perturbed_radius = r + pert_func((point)...)
    
    
        # Normalize to be on perturbed sphere
        point *= perturbed_radius / sqrt((point[1] - origo[1])^2 + (point[2] - origo[2])^2 + (point[3] - origo[3])^2)
        append!(pointTags, gmsh.model.occ.addPoint(point...))
    end

    append!(pointTags, end_tag)
    BSpline = gmsh.model.occ.addBSpline(pointTags)
    return BSpline

end


function create_perturbed_surface(brain::geo3D, r, pert_func, BS_points)
    # 2D func porjected to sphere surface
    vertex = sphere_patch_corners(r, brain.angle, pert_func) # Get sphere patch corners
    arc = Array{Int32,1}(undef, 4)
    gmsh.model.occ.synchronize()

    arc[1] = add_perturbed_arc(vertex[1], brain.origo, vertex[2], pert_func, BS_points)
    arc[2] = add_perturbed_arc(vertex[2], brain.origo, vertex[3], pert_func, BS_points)
    arc[3] = add_perturbed_arc(vertex[3], brain.origo, vertex[4], pert_func, BS_points)
    arc[4] = add_perturbed_arc(vertex[4], brain.origo, vertex[1], pert_func, BS_points)

    CurveLoop = gmsh.model.occ.addCurveLoop([arc[1], arc[2], arc[3], arc[4]])
    surf = gmsh.model.occ.addBSplineFilling(CurveLoop, -1, "Stretch") # Alternatively use "Coons"
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [surf])

    return vertex, arc, surf

end


function connect_and_volumize(brain::geo3D)
    # Connect surfaces and create volumes
    for i in 1:2
        vline = [gmsh.model.occ.addLine(brain.vertex[i, 1], brain.vertex[i+1, 1])] # Vertical line
        Loop = [] # curve loops
        for j in 2:4
            append!(vline, gmsh.model.occ.addLine(brain.vertex[i, j], brain.vertex[i+1, j]))
            append!(Loop, gmsh.model.occ.addCurveLoop([vline[j-1], brain.arc[i+1, j-1], -vline[j], -brain.arc[i, j-1]]))
        end
        append!(Loop, gmsh.model.occ.addCurveLoop([vline[4], brain.arc[i+1, 4], -vline[1], -brain.arc[i, 4]]))

        brain.rad_surf[i, :] = [gmsh.model.occ.addSurfaceFilling(l) for l in Loop]

        surfLoop = gmsh.model.occ.addSurfaceLoop([brain.tan_surf[i], brain.rad_surf[i]..., brain.tan_surf[i+1]])
        Vol = gmsh.model.occ.addVolume([surfLoop])

        gmsh.model.occ.synchronize()
        [gmsh.model.addPhysicalGroup(2, [s]) for s in brain.rad_surf[i]]
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
end


function apply_periodic_meshing(brain::geo3D)
    θy, θx = brain.angle # Notice θx and θy is purposely switched 
    zx_angle, zy_angle = brain.angle
    println(θy)
    println(θx)



    # Rotation around x-axis 
    Rx = [1, 0, 0, 0,
        0, cos(θy), sin(θy), 0,
        0, -sin(θy), cos(θy), 0,
        0, 0, 0, 1]

    # Rotation around y-axis 
    Ry = [cos(θx), 0, -sin(θx), 0,
        0, 1, 0, 0,
        sin(θx), 0, cos(θx), 0,
        0, 0, 0, 1]




    # gmsh.model.mesh.setPeriodic(2, [brain.rad_surf[1, 1]], [brain.rad_surf[1, 3]], Rx) 





    # # Manual test
    # point3 = gmsh.model.getValue(0, 3, []) #dim, tag, parametrization
    # cord = [point3..., 1]
    # θ = -θy
    # A = [1 0 0 0; 0 cos(θ) sin(θ) 0; 0 -sin(θ) cos(θ) 0; 0 0 0 1]
    # B = [cos(θ) 0 -sin(θ) 0; 0 1 0 0; sin(θ) 0 cos(θ) 0; 0 0 0 1]
    # res = A * cord
    # println(cord)
    # println(res)

    

    # for i in 1:2
    #     gmsh.model.mesh.setPeriodic(2, [brain.rad_surf[i, 1]], [brain.rad_surf[i, 3]], Ry) # From -y → y surf
    #     gmsh.model.mesh.setPeriodic(2, [brain.rad_surf[i, 2]], [brain.rad_surf[i, 4]], Rx) # From -x → x surf        
    # end



end


function create_brain_3D(param::model_params, view=true)
    #?--> Safety check of parameters? # Should enforce  0 < angle < pi

    gmsh.initialize(["", "-clmax", string(param.lc)])
    brain = geo3D() # Struct for holding tags and angle

    # Calculate derived parameters
    rI = param.r_curv - param.r_brain                   # Inner radius  
    rD = param.r_curv - param.d_ratio * param.r_brain   # Radius for dividing line
    x_arcLen, y_arcLen = param.arcLen
    zx_angle = x_arcLen / param.r_curv                  # Angle span z -> x -axis
    zy_angle = y_arcLen / param.r_curv                  # Angle span z -> x -axis


    brain.angle = (zx_angle, zy_angle)

    # Add radial surfaces
    brain.vertex[1, :], brain.arc[1, :], brain.tan_surf[1] = create_surface(brain, rI)
    brain.vertex[2, :], brain.arc[2, :], brain.tan_surf[2] = create_perturbed_surface(brain, rD, param.inner_perturb, param.BS_points)
    # brain.vertex[2, :], brain.arc[2, :], brain.tan_surf[2] = create_surface(brain, rD)

    brain.vertex[3, :], brain.arc[3, :], brain.tan_surf[3] = create_perturbed_surface(brain, param.r_curv, param.outer_perturb, param.BS_points)


    connect_and_volumize(brain)
    add_mesh_field(brain, param)
    apply_periodic_meshing(brain)



    # Generate mesh
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(1)


    # View and finalize
    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    gmsh.finalize()




end # End of create_brain_3D




if abspath(PROGRAM_FILE) == @__FILE__

    lc = 0.5
    arcLen = (5, 5)
    r_brain = 5
    d_ratio = 0.5
    r_curv = 20
    # inner_perturb(x, y, z) = 0.2 * sin(pi * x / 2) + 0.2 * sin(pi * y / 2)
    inner_perturb(x, y, z) = 0.05*x + 0.05*y ### Why isn't the lines continous with this? FIGURE THIS OUT MONDAY :D

    outer_perturb(x, y, z) = 0.2 * sin(pi * x / 2) + 0.2 * sin(pi * y / 1)
    BS_points = 50 # Make direction depending x,y
    field_Lc_lim = [1 / 2, 1]
    field_Dist_lim = [0.3, 0.5]

    param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, field_Lc_lim, field_Dist_lim)

    A = Array{Int32,1}(undef, 3)
    B = Array{Int32,2}(undef, 2, 4)



    create_brain_3D(param)
end



