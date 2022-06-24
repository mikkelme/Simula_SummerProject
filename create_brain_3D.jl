using GridapGmsh: gmsh
using Printf

# # struct that can be changed along the way: 
# mutable struct Person
#     name::String # Fields
#     age::Float64
#     isActive

#     function Person(name, age) # Think this is a constructor
#         new(name, age, true) # sets isActive=true bu default        
#     end
# end

# logan = Person("Logan", 44)
# logan.age += 1

# function birtday(person::Person)
#     person.age += 1
# end


struct model_params
    lc::Float64             # Mesh size
    arcLen::Tuple{Float64,Float64}  # Outer arc length (x, y)-direction
    r_brain::Float64        # Radial total length
    d_ratio::Float64        # Radial relative length of fluid section
    r_curv::Float64         # Radius of curvature
    inner_perturb::Function # Radial perturbation on inner surface f(x,y,z)
    outer_perturb::Function # Radial perturbation on inner surface f(x,y,z)
    BS_points::Int64        # Number of points in BSpline curves
end



mutable struct geo3D
    origo::Int32
    vertex::Array{Int32,2}
    arc::Array{Int32,2}
    rad_surf::Array{Int32, 1} # Radial surfaces

    function geo3D() # Constructor
        new(0, Matrix{Int32}(undef, 3, 4), Matrix{Int32}(undef, 3, 4), Vector{Int32}(undef, 3))
    end
end


function spherical_to_cartesian(r, theta, phi)
    # theta: angle from x-axis (0 to pi)
    # phi: angle from z-axis (0 to 2pi)
    x = r * cos(theta) * sin(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(phi)
    return x, y, z
end




function create_surface(brain::geo3D, r, angle)
    # Make boundary of surface
    ###################
    # Use inclined circle (for y incline)
    # eq.: x^2 + (z + αy)^2 = r^2
    ###################s
    

    zx_angle, zy_angle = angle


    vertex = Array{Int32,1}(undef, 4)
    arc = Array{Int32,1}(undef, 4)


    x, _, zx = spherical_to_cartesian(r, 0, zx_angle / 2)
    _, y, zy = spherical_to_cartesian(r, pi / 2, zy_angle / 2)

    alpha_x = (1 - zx) / x
    alpha_y = (1 - zy) / y


    a = 2 * alpha_x / (1 - alpha_x^2)
    b = 2 * alpha_y / (1 - alpha_y^2)
    c = r / sqrt(a^2 + b^2 + 1)


    vertex[1] = gmsh.model.occ.addPoint(a * c, b * c, c)
    vertex[2] = gmsh.model.occ.addPoint(-a * c, b * c, c)
    vertex[3] = gmsh.model.occ.addPoint(-a * c, -b * c, c)
    vertex[4] = gmsh.model.occ.addPoint(a * c, -b * c, c)


    arc[1] = gmsh.model.occ.addCircleArc(vertex[1], brain.origo, vertex[2])
    arc[2] = gmsh.model.occ.addCircleArc(vertex[2], brain.origo, vertex[3])
    arc[3] = gmsh.model.occ.addCircleArc(vertex[3], brain.origo, vertex[4])
    arc[4] = gmsh.model.occ.addCircleArc(vertex[4], brain.origo, vertex[1])

    CurveLoop = gmsh.model.occ.addCurveLoop([arc[1], arc[2], arc[3], arc[4]])
    surf = gmsh.model.occ.addSurfaceFilling(CurveLoop)

    return vertex, arc, surf

end


function add_perturbed_arc(start_tag, center_tag, end_tag, pert_func, BS_points)
    


    pointTags = [start_tag]
    start_point = gmsh.model.getValue(0, start_tag, []) #dim, tag, parametrization
    end_point = gmsh.model.getValue(0, end_tag, [])
    origo = gmsh.model.getValue(0, center_tag, [])
    r = sqrt(start_point[1]^2 + start_point[2]^2 + start_point[3]^2)


    direction = end_point - start_point
    range = LinRange(0, 1, BS_points)
    for i in 2:BS_points-1
        # Step toward end_point in straight line
        point = start_point + range[i] * direction

        # Add perturbation
        perturbed_radius = r + pert_func((point - start_point)...)

        # Normalize to be on perturbed sphere
        point *= perturbed_radius / sqrt((point[1] - origo[1])^2 + (point[2] - origo[2])^2 + (point[3] - origo[3])^2)

        append!(pointTags, gmsh.model.occ.addPoint(point...))
    end

    append!(pointTags, end_tag)
    BSpline = gmsh.model.occ.addBSpline(pointTags)
    return BSpline

end





function create_perturbed_surface(brain::geo3D, r, angle, pert_func, BS_points)
    
    zx_angle, zy_angle = angle



    vertex = Array{Int32,1}(undef, 4)
    arc = Array{Int32,1}(undef, 4)


    x, _, zx = spherical_to_cartesian(r, 0, zx_angle / 2)
    _, y, zy = spherical_to_cartesian(r, pi / 2, zy_angle / 2)

    alpha_x = (1 - zx) / x
    alpha_y = (1 - zy) / y



    a = 2 * alpha_x / (1 - alpha_x^2)
    b = 2 * alpha_y / (1 - alpha_y^2)
    c = r / sqrt(a^2 + b^2 + 1)


    vertex[1] = gmsh.model.occ.addPoint(a * c, b * c, c)
    vertex[2] = gmsh.model.occ.addPoint(-a * c, b * c, c)
    vertex[3] = gmsh.model.occ.addPoint(-a * c, -b * c, c)
    vertex[4] = gmsh.model.occ.addPoint(a * c, -b * c, c)
    gmsh.model.occ.synchronize()



    arc[1] = add_perturbed_arc(vertex[1], brain.origo, vertex[2], pert_func, BS_points)
    arc[2] = add_perturbed_arc(vertex[2], brain.origo, vertex[3], pert_func, BS_points)
    arc[3] = add_perturbed_arc(vertex[3], brain.origo, vertex[4], pert_func, BS_points)
    arc[4] = add_perturbed_arc(vertex[4], brain.origo, vertex[1], pert_func, BS_points)

    CurveLoop = gmsh.model.occ.addCurveLoop([arc[1], arc[2], arc[3], arc[4]])
    surf = gmsh.model.occ.addBSplineFilling(CurveLoop, -1, "Stretch") # Alternatively use "Coons"


    gmsh.model.occ.synchronize()



    return vertex, arc, surf

end


function connect_and_fill(brain::geo3D)
    

    for i in 1:2    
        vline = [gmsh.model.occ.addLine(brain.vertex[i, 1], brain.vertex[i+1, 1])] # Vertical line
        Loop = [] # curve loops
        for j in 2:4
            append!(vline, gmsh.model.occ.addLine(brain.vertex[i, j], brain.vertex[i+1, j]))
            append!(Loop, gmsh.model.occ.addCurveLoop([vline[j-1], brain.arc[i+1, j-1], -vline[j], -brain.arc[i, j-1]]))
        end
        append!(Loop, gmsh.model.occ.addCurveLoop([vline[4], brain.arc[i+1, 4], -vline[1], -brain.arc[i, 4]]))
    
        Surf = [gmsh.model.occ.addSurfaceFilling(l) for l in Loop]
        SurfLoop = gmsh.model.occ.addSurfaceLoop([brain.rad_surf[i], Surf..., brain.rad_surf[i+1]])
        Vol = gmsh.model.occ.addVolume([SurfLoop])
    end





end

function create_brain_3D(param::model_params, view=true)
    # Consider mutable struct for stuff in here
    # Check parameters

    brain = geo3D()

    rI = param.r_curv - param.r_brain                   # Inner radius  
    rD = param.r_curv - param.d_ratio * param.r_brain   # Radius for dividing line

    # Should enforce  0 < angle < pi
    x_arcLen, y_arcLen = param.arcLen
    zx_angle = x_arcLen / param.r_curv                  # Angle span z -> x -axis
    zy_angle = y_arcLen / param.r_curv                  # Angle span z -> x -axis
    angle = (zx_angle, zy_angle)


    #--- occmetry ---# (Only curved slab for now)
    gmsh.initialize(["", "-clmax", string(0.1)])
    

    brain.origo = gmsh.model.occ.addPoint(0.0, 0.0, 0.0)



    brain.vertex[1, :], brain.arc[1, :], brain.rad_surf[1] = create_surface(brain, rI, angle)
    brain.vertex[2, :], brain.arc[2, :], brain.rad_surf[2] = create_perturbed_surface(brain, rD, angle, param.inner_perturb, param.BS_points)
    brain.vertex[3, :], brain.arc[3, :], brain.rad_surf[3] = create_perturbed_surface(brain, param.r_curv, angle, param.outer_perturb, param.BS_points)


    connect_and_fill(brain)





    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)





    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    gmsh.finalize()








end # End of create_brain_3D




if abspath(PROGRAM_FILE) == @__FILE__

    lc = 0.5
    arcLen = (5.0, 2.0)
    r_brain = 5
    d_ratio = 0.5
    r_curv = 15
    inner_perturb(x, y, z) = 0.2 * sin(pi * x / 0.1) + 0.2 * sin(pi * y / 0.5)
    outer_perturb(x, y, z) = 0.2 * sin(pi * x / 2) + 0.2 * sin(pi * y / 1)
    BS_points = 50 # Make direction depending x,y
    param = model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points)

    A = Array{Int32,1}(undef, 3)
    B = Array{Int32,2}(undef, 2, 4)



    create_brain_3D(param)
end





##### Leftovers #####



#Draw cross for reference
# X = gmsh.model.occ.addPoint(spherical_to_cartesian(r, 0, zx_angle / 2)...)
# Y = gmsh.model.occ.addPoint(spherical_to_cartesian(r, pi / 2, zy_angle / 2)...)
# Xneg = gmsh.model.occ.addPoint(spherical_to_cartesian(r, pi, zx_angle / 2)...)
# Yneg = gmsh.model.occ.addPoint(spherical_to_cartesian(r, 3 * pi / 2, zy_angle / 2)...)
# z_sphere = gmsh.model.occ.addPoint(0, 0, 1)

# # gmsh.model.occ.synchronize()

# # # Helpful naming
# # X_G = model.addPhysicalGroup(0, [X])
# # model.setPhysicalName(0, X_G, "X_G")

# # Y_G = model.addPhysicalGroup(0, [Y])
# # model.setPhysicalName(0, Y_G, "Y_G")

# # Xneg_G = model.addPhysicalGroup(0, [Xneg])
# # model.setPhysicalName(0, Xneg_G, "Xneg_G")

# # Yneg_G = model.addPhysicalGroup(0, [Yneg])
# # model.setPhysicalName(0, Yneg_G, "Yneg_G")


# gmsh.model.occ.addCircleArc(X, origo, Xneg)
# gmsh.model.occ.addCircleArc(Y, origo, Yneg)

# gmsh.model.occ.addLine(origo, X)
# gmsh.model.occ.addLine(origo, Y)
# gmsh.model.occ.addLine(origo, Xneg)
# gmsh.model.occ.addLine(origo, Yneg)
# gmsh.model.occ.addLine(origo, z_sphere)

