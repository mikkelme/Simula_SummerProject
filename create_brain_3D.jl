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


mutable struct model_params
    lc::Float64             # Mesh size
    arcLen::Tuple{Float64,Float64}  # Outer arc length (x, y)-direction
    r_brain::Float64        # Radial total length
    d_ratio::Float64        # Radial relative length of fluid section
    r_curv::Float64         # Radius of curvature
    inner_perturb::Function # Radial perturbation on inner surface f(x,y,z)
    outer_perturb::Function # Radial perturbation on inner surface f(x,y,z)
    BS_points::Int64        # Number of points in BSpline curves
    vertex::Array{Int64, 2}
    # origo::Int64

    function model_params(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points) # Constructor
        new(lc, arcLen, r_brain, d_ratio, r_curv, inner_perturb, outer_perturb, BS_points, Matrix{Int64}(undef, 2, 4) )
    end
    # ::Array{Float64,3}
end


function spherical_to_cartesian(r, theta, phi)
    # theta: angle from x-axis (0 to pi)
    # phi: angle from z-axis (0 to 2pi)
    x = r * cos(theta) * sin(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(phi)
    return x, y, z
end




function create_surface(param, r, angle)
    model = gmsh.model

    zx_angle, zy_angle = angle

    # origo = model.occ.addPoint(0.0, 0.0, 0.0) # Center ### Get rid of this later
    origo = param.vertex[1,1]

    #-----> Working form here with using param struct <------

    # Make boundary of surface
    ###################
    # alpha = (1 - z) / y
    # Use inclined circle (for y incline)
    # eq.: x^2 + (z + αy)^2 = r^2
    ###################

    vertex = Array{Float64,1}(undef, 4)
    arc = Array{Float64,1}(undef, 4)




    x, _, zx = spherical_to_cartesian(r, 0, zx_angle / 2)
    _, y, zy = spherical_to_cartesian(r, pi / 2, zy_angle / 2)

    alpha_x = (1 - zx) / x
    alpha_y = (1 - zy) / y


    a = 2 * alpha_x / (1 - alpha_x^2)
    b = 2 * alpha_y / (1 - alpha_y^2)
    c = r / sqrt(a^2 + b^2 + 1)


    vertex[1] = model.occ.addPoint(a * c, b * c, c)
    vertex[2] = model.occ.addPoint(-a * c, b * c, c)
    vertex[3] = model.occ.addPoint(-a * c, -b * c, c)
    vertex[4] = model.occ.addPoint(a * c, -b * c, c)


    arc[1] = model.occ.addCircleArc(vertex[1], origo, vertex[2])
    arc[2] = model.occ.addCircleArc(vertex[2], origo, vertex[3])
    arc[3] = model.occ.addCircleArc(vertex[3], origo, vertex[4])
    arc[4] = model.occ.addCircleArc(vertex[4], origo, vertex[1])

    CurveLoop = model.occ.addCurveLoop([arc[1], arc[2], arc[3], arc[4]])
    surf = model.occ.addSurfaceFilling(CurveLoop)


    return vertex, arc, surf

end


function add_perturbed_arc(param, start_tag, center_tag, end_tag, pert_func, BS_points)
    model = gmsh.model


    pointTags = [start_tag]
    start_point = model.getValue(0, start_tag, []) #dim, tag, parametrization
    end_point = model.getValue(0, end_tag, [])
    origo = model.getValue(0, center_tag, [])
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

        append!(pointTags, model.occ.addPoint(point...))
    end

    append!(pointTags, end_tag)
    BSpline = model.occ.addBSpline(pointTags)
    return BSpline

end





function create_perturbed_surface(param, r, angle, pert_func, BS_points)
    model = gmsh.model
    zx_angle, zy_angle = angle

    origo = model.occ.addPoint(0.0, 0.0, 0.0) # Center

    vertex = Array{Float64,1}(undef, 4)
    arc = Array{Float64,1}(undef, 4)


    x, _, zx = spherical_to_cartesian(r, 0, zx_angle / 2)
    _, y, zy = spherical_to_cartesian(r, pi / 2, zy_angle / 2)

    alpha_x = (1 - zx) / x
    alpha_y = (1 - zy) / y



    a = 2 * alpha_x / (1 - alpha_x^2)
    b = 2 * alpha_y / (1 - alpha_y^2)
    c = r / sqrt(a^2 + b^2 + 1)


    vertex[1] = model.occ.addPoint(a * c, b * c, c)
    vertex[2] = model.occ.addPoint(-a * c, b * c, c)
    vertex[3] = model.occ.addPoint(-a * c, -b * c, c)
    vertex[4] = model.occ.addPoint(a * c, -b * c, c)
    model.occ.synchronize()



    arc[1] = add_perturbed_arc(model, vertex[1], origo, vertex[2], pert_func, BS_points)
    arc[2] = add_perturbed_arc(model, vertex[2], origo, vertex[3], pert_func, BS_points)
    arc[3] = add_perturbed_arc(model, vertex[3], origo, vertex[4], pert_func, BS_points)
    arc[4] = add_perturbed_arc(model, vertex[4], origo, vertex[1], pert_func, BS_points)

    CurveLoop = model.occ.addCurveLoop([arc[1], arc[2], arc[3], arc[4]])
    surf = model.occ.addBSplineFilling(CurveLoop, 2, "Stretch") # Alternatively use "Coons" is also nice


    model.occ.synchronize()



    return vertex, arc, surf

end




function create_brain_3D(param::model_params, view=true)
    # Consider mutable struct for stuff in here
    # Check parameters



    rI = param.r_curv - param.r_brain                         # Inner radius  
    rD = param.r_curv - param.d_ratio * param.r_brain        # Radius for dividing line

    # Should enforce  0 < angle < pi
    x_arcLen, y_arcLen = param.arcLen
    zx_angle = x_arcLen / param.r_curv # Angle span from z-axis towards x-axis
    zy_angle = y_arcLen / param.r_curv # Angle span from z-axis towards x-axis
    angle = (zx_angle, zy_angle)


    #--- occmetry ---# (Only curved slab for now)
    gmsh.initialize(["", "-clmax", string(0.1)])
    model = gmsh.model

    origo = model.occ.addPoint(0.0, 0.0, 0.0)
    param.vertex[1,1] = origo

    # println(typeof(model))
    # exit()
    I_vertex, I_arc, I_surf = create_surface(param, rI, angle)

    D_vertex, D_arc, D_surf = create_perturbed_surface(param, rD, angle, param.inner_perturb, param.BS_points)
    # O_vertex, O_arc, O_surf = create_perturbed_surface(model, param.r_curv, angle, param.outer_perturb, param.BS_points)


    #--- Connect surfaces ---#
    # Lines
    # A = model.occ.addLine(I_vertex[1], D_vertex[1])
    # B = model.occ.addLine(I_vertex[2], D_vertex[2])
    # C = model.occ.addLine(I_vertex[3], D_vertex[3])
    # D = model.occ.addLine(I_vertex[4], D_vertex[4])


    # loop1 = model.occ.addCurveLoop([A, D_arc[1], -B, -I_arc[1]])
    # Y_surf = model.occ.addPlaneSurface([loop1])

    # loop2 = model.occ.addCurveLoop([B, D_arc[2], -C, -I_arc[2]])
    # Xneg_surf = model.occ.addPlaneSurface([loop2])

    # loop3 = model.occ.addCurveLoop([C, D_arc[3], -D, -I_arc[3]])
    # Xneg_surf = model.occ.addPlaneSurface([loop3])

    # loop4 = model.occ.addCurveLoop([D, D_arc[4], -A, -I_arc[4]])
    # Yneg_surf = model.occ.addPlaneSurface([loop4])






    model.occ.synchronize()
    model.mesh.generate(2)





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

    create_brain_3D(param)
end





##### Leftovers #####



#Draw cross for reference
# X = model.occ.addPoint(spherical_to_cartesian(r, 0, zx_angle / 2)...)
# Y = model.occ.addPoint(spherical_to_cartesian(r, pi / 2, zy_angle / 2)...)
# Xneg = model.occ.addPoint(spherical_to_cartesian(r, pi, zx_angle / 2)...)
# Yneg = model.occ.addPoint(spherical_to_cartesian(r, 3 * pi / 2, zy_angle / 2)...)
# z_sphere = model.occ.addPoint(0, 0, 1)

# # model.occ.synchronize()

# # # Helpful naming
# # X_G = model.addPhysicalGroup(0, [X])
# # model.setPhysicalName(0, X_G, "X_G")

# # Y_G = model.addPhysicalGroup(0, [Y])
# # model.setPhysicalName(0, Y_G, "Y_G")

# # Xneg_G = model.addPhysicalGroup(0, [Xneg])
# # model.setPhysicalName(0, Xneg_G, "Xneg_G")

# # Yneg_G = model.addPhysicalGroup(0, [Yneg])
# # model.setPhysicalName(0, Yneg_G, "Yneg_G")


# model.occ.addCircleArc(X, origo, Xneg)
# model.occ.addCircleArc(Y, origo, Yneg)

# model.occ.addLine(origo, X)
# model.occ.addLine(origo, Y)
# model.occ.addLine(origo, Xneg)
# model.occ.addLine(origo, Yneg)
# model.occ.addLine(origo, z_sphere)

