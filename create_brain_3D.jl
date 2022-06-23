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
    lc::Float64         # Mesh size
    x_arcLen::Float64   # Outer arc length (x-direction)
    y_arcLen::Float64   # Outer arc length (y-direction)
    r_brain::Float64    # Radial total length
    d_ratio::Float64    # Radial relative length of fluid section
    r_curv::Float64     # Radius of curvature
    BS_points::Int64    # Number of points in BSpline curves
end


function spherical_to_cartesian(r, theta, phi)
    # theta: angle from x-axis (0 to pi)
    # phi: angle from z-axis (0 to 2pi)
    x = r * cos(theta) * sin(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(phi)
    return x, y, z
end

function cartesian_to_spherical(x, y, z)


end


function create_surface(model, r, zx_angle, zy_angle)

    origo = model.geo.addPoint(0.0, 0.0, 0.0) # Center

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


    vertex[1] = model.geo.addPoint(a * c, b * c, c)
    vertex[2] = model.geo.addPoint(-a * c, b * c, c)
    vertex[3] = model.geo.addPoint(-a * c, -b * c, c)
    vertex[4] = model.geo.addPoint(a * c, -b * c, c)


    arc[1] = model.geo.addCircleArc(vertex[1], origo, vertex[2])
    arc[2] = model.geo.addCircleArc(vertex[2], origo, vertex[3])
    arc[3] = model.geo.addCircleArc(vertex[3], origo, vertex[4])
    arc[4] = model.geo.addCircleArc(vertex[4], origo, vertex[1])

    CurveLoop = model.geo.addCurveLoop([arc[1], arc[2], arc[3], arc[4]])
    surf = model.geo.addSurfaceFilling([CurveLoop])


    return vertex, arc, surf



# function create_perturbed_surface(model, r, zx_angle, zy_angle)

   
#     vertex = Array{Float64,1}(undef, 4)
#     arc = Array{Float64,1}(undef, 4)




#     x, _, zx = spherical_to_cartesian(r, 0, zx_angle / 2)
#     _, y, zy = spherical_to_cartesian(r, pi / 2, zy_angle / 2)

#     alpha_x = (1 - zx) / x
#     alpha_y = (1 - zy) / y


#     a = 2 * alpha_x / (1 - alpha_x^2)
#     b = 2 * alpha_y / (1 - alpha_y^2)
#     c = r / sqrt(a^2 + b^2 + 1)


#     vertex[1] = model.geo.addPoint(a * c, b * c, c)
#     vertex[2] = model.geo.addPoint(-a * c, b * c, c)
#     vertex[3] = model.geo.addPoint(-a * c, -b * c, c)
#     vertex[4] = model.geo.addPoint(a * c, -b * c, c)


#     arc[1] = model.geo.addCircleArc(vertex[1], origo, vertex[2])
#     arc[2] = model.geo.addCircleArc(vertex[2], origo, vertex[3])
#     arc[3] = model.geo.addCircleArc(vertex[3], origo, vertex[4])
#     arc[4] = model.geo.addCircleArc(vertex[4], origo, vertex[1])

#     CurveLoop = model.geo.addCurveLoop([arc[1], arc[2], arc[3], arc[4]])
#     surf = model.geo.addSurfaceFilling([CurveLoop])


#     return vertex, arc, surf





end

function create_brain_3D(params::model_params, view=true)
    # Consider mutable struct for stuff in here
    # Check parameters


    rI = params.r_curv - params.r_brain                         # Inner radius  
    rD = params.r_curv - params.d_ratio * params.r_brain        # Radius for dividing line

    # Should enforce  0 < angle < pi
    zx_angle = params.x_arcLen / params.r_curv # Angle span from z-axis towards x-axis
    zy_angle = params.y_arcLen / params.r_curv # Angle span from z-axis towards x-axis


    #--- Geometry ---# (Only curved slab for now)
    gmsh.initialize(["", "-clmax", string(0.1)])
    model = gmsh.model



    I_vertex, I_arc, I_surf = create_surface(model, rI, zx_angle, zy_angle)
    D_vertex, D_arc, D_surf = create_surface(model, rD, zx_angle, zy_angle)
    # create_surface(model, params.r_curv, zx_angle, zy_angle)


    #--- Connect surfaces ---#
    # Lines
    A = model.geo.addLine(I_vertex[1], D_vertex[1])
    B = model.geo.addLine(I_vertex[2], D_vertex[2])
    C = model.geo.addLine(I_vertex[3], D_vertex[3])
    D = model.geo.addLine(I_vertex[4], D_vertex[4])


    loop1 = model.geo.addCurveLoop([A, D_arc[1], -B, -I_arc[1]])
    Y_surf = model.geo.addPlaneSurface([loop1])

    loop2 = model.geo.addCurveLoop([B, D_arc[2], -C, -I_arc[2]])
    Xneg_surf = model.geo.addPlaneSurface([loop2])

    loop3 = model.geo.addCurveLoop([C, D_arc[3], -D, -I_arc[3]])
    Xneg_surf = model.geo.addPlaneSurface([loop3])

    loop4 = model.geo.addCurveLoop([D, D_arc[4], -A, -I_arc[4]])
    Yneg_surf = model.geo.addPlaneSurface([loop4])






    model.geo.synchronize()
    model.mesh.generate(2)





    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    gmsh.finalize()








end # End of create_brain_3D




if abspath(PROGRAM_FILE) == @__FILE__

    lc = 0.5
    x_arcLen = 5
    y_arcLen = 2
    r_brain = 5
    d_ratio = 0.5
    r_curv = 15
    BS_points = 50
    params = model_params(lc, x_arcLen, y_arcLen, r_brain, d_ratio, r_curv, BS_points)

    create_brain_3D(params)
end





##### Leftvers #####



#Draw cross for reference
# X = model.geo.addPoint(spherical_to_cartesian(r, 0, zx_angle / 2)...)
# Y = model.geo.addPoint(spherical_to_cartesian(r, pi / 2, zy_angle / 2)...)
# Xneg = model.geo.addPoint(spherical_to_cartesian(r, pi, zx_angle / 2)...)
# Yneg = model.geo.addPoint(spherical_to_cartesian(r, 3 * pi / 2, zy_angle / 2)...)
# z_sphere = model.geo.addPoint(0, 0, 1)

# # model.geo.synchronize()

# # # Helpful naming
# # X_G = model.addPhysicalGroup(0, [X])
# # model.setPhysicalName(0, X_G, "X_G")

# # Y_G = model.addPhysicalGroup(0, [Y])
# # model.setPhysicalName(0, Y_G, "Y_G")

# # Xneg_G = model.addPhysicalGroup(0, [Xneg])
# # model.setPhysicalName(0, Xneg_G, "Xneg_G")

# # Yneg_G = model.addPhysicalGroup(0, [Yneg])
# # model.setPhysicalName(0, Yneg_G, "Yneg_G")


# model.geo.addCircleArc(X, origo, Xneg)
# model.geo.addCircleArc(Y, origo, Yneg)

# model.geo.addLine(origo, X)
# model.geo.addLine(origo, Y)
# model.geo.addLine(origo, Xneg)
# model.geo.addLine(origo, Yneg)
# model.geo.addLine(origo, z_sphere)

