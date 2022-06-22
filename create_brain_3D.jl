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
    arcLen::Float64     # Outer arc length 
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


function create_surface(model, r, zx, zy)

    # x, y, z = spherical_to_cartesian(r, theta / 2, phi / 2)
    origo = model.geo.addPoint(0.0, 0.0, 0.0) # Center

    zx_angle = pi / 3
    zy_angle = pi / 4


    # Make boundary of surface

    x, _, rx = spherical_to_cartesian(r, 0, zx_angle / 2)
    _, y, ry = spherical_to_cartesian(r, pi / 2, zy_angle / 2)

    a = sqrt(r^2 - rx^2)
    b = sqrt(r^2 - ry^2)
    c = sqrt(rx^2 + ry^2 - r^2)


    ###################
    alpha = (1 - z) / y
    # Use inclined circle (for y incline)
    # eq.: x^2 + (z + αy)^2 = r^2
    ###################
    


    # x -= 0.04
    # y -= 0.04
    # theta = atan(y / x)
    # phi = asin(x / (r * cos(theta)))

    A = model.geo.addPoint(a, b, c)
    B = model.geo.addPoint(-a, b, c)
    C = model.geo.addPoint(-a, -b, c)
    D = model.geo.addPoint(a, -b, c)


    # A = model.geo.addPoint(spherical_to_cartesian(r, theta, phi)...)
    # B = model.geo.addPoint(spherical_to_cartesian(r, -theta + pi, phi)...)
    # C = model.geo.addPoint(spherical_to_cartesian(r, theta + pi, phi)...)
    # D = model.geo.addPoint(spherical_to_cartesian(r, -theta, phi)...)

    AB_arc = model.geo.addCircleArc(A, origo, B)
    BC_arc = model.geo.addCircleArc(B, origo, C)
    CD_arc = model.geo.addCircleArc(C, origo, D)
    DA_arc = model.geo.addCircleArc(D, origo, A)

    #AB midpoint
    # x1, y1, z1 = spherical_to_cartesian(r, theta, phi)
    # x2, y2, z2 = spherical_to_cartesian(r, -theta + pi, phi)

    # AB_mid = [x1 + x2, y1 + y2, z1 + z2]
    # theta_mid = atan((y1 + y2) / (x1 + x2))
    # phi_mid = asin((x1 + x2) / (2*r * cos(theta_mid)))
    # theta_mid = (theta -theta + pi/2)
    # AB_mid = model.geo.addPoint(spherical_to_cartesian(r, theta_mid, phi)...)



    # B = model.geo.addPoint(spherical_to_cartesian(r, theta + pi/2, phi)...)
    # C = model.geo.addPoint(spherical_to_cartesian(r, theta + pi, phi)...)
    # D = model.geo.addPoint(spherical_to_cartesian(r, theta + 3*pi/2, phi)...)


    # phi1 = asin(x/(r*cos(theta)))
    # phi2 = asin(y / (r * sin(theta)))
    # println(phi1*pi)
    # println(phi2*pi)


    # Circle around y-axis
    #



    #Draw cross for reference
    X = model.geo.addPoint(spherical_to_cartesian(r, 0, zx_angle / 2)...)
    Y = model.geo.addPoint(spherical_to_cartesian(r, pi / 2, zy_angle / 2)...)
    Xneg = model.geo.addPoint(spherical_to_cartesian(r, pi, zx_angle / 2)...)
    Yneg = model.geo.addPoint(spherical_to_cartesian(r, 3 * pi / 2, zy_angle / 2)...)
    z_sphere = model.geo.addPoint(0, 0, 1)

    model.geo.synchronize()

    # Helpful naming
    X_G = model.addPhysicalGroup(0, [X])
    model.setPhysicalName(0, X_G, "X_G")

    Y_G = model.addPhysicalGroup(0, [Y])
    model.setPhysicalName(0, Y_G, "Y_G")

    Xneg_G = model.addPhysicalGroup(0, [Xneg])
    model.setPhysicalName(0, Xneg_G, "Xneg_G")

    Yneg_G = model.addPhysicalGroup(0, [Yneg])
    model.setPhysicalName(0, Yneg_G, "Yneg_G")



    model.geo.addCircleArc(X, origo, Xneg)
    model.geo.addCircleArc(Y, origo, Yneg)

    model.geo.addLine(origo, X)
    model.geo.addLine(origo, Y)
    model.geo.addLine(origo, Xneg)
    model.geo.addLine(origo, Yneg)
    model.geo.addLine(origo, z_sphere)




    #     arc1 = model.geo.addCircleArc(IL, C, IR) # Inner arc






    # point2 = model.geo.addPoint(spherical_to_cartesian(r, pi / 2 + theta / 2, phi / 2)...)
    # point3 = model.geo.addPoint(spherical_to_cartesian(r, pi / 2 - theta / 2, phi / 2)...)
    # point4 = model.geo.addPoint(spherical_to_cartesian(r, pi / 2 + pi + theta / 2, phi / 2)...)
    # point5 = model.geo.addPoint(spherical_to_cartesian(r, pi / 2 + pi - theta / 2, phi / 2)...)

    # line1 = model.geo.addLine(point2, point3)
    # line2 = model.geo.addLine(point3, point4)
    # line3 = model.geo.addLine(point4, point5)
    # line4 = model.geo.addLine(point5, point2)









    # sph = model.occ.addSphere(0, 0, 0, r, -1, -2*pi, pi/2, 2*pi)



    # model.occ.synchronize()
    model.geo.synchronize()

    # R = 1
    # R1 = 0.95
    # sph = model.occ.addSphere(0, 0, 0, R, -1, 0, pi / 2, pi / 2)
    # b1 = model.occ.addBox(R1, 0, 0, R, R, R)
    # b2 = model.occ.addBox(0, R1, 0, R, R, R)
    # b3 = model.occ.addBox(0, 0, R1, R, R, R)
    # model.occ.cut([(3, sph)], [(3, b1), (3, b2), (3, b3)])
    # model.occ.synchronize()
    # model.removeEntities([(3, sph)])
    # model.removeEntities([(2, 2), (2, 4), (2, 6)], true)






end

function create_brain_3D(params::model_params, view=true)
    # Consider mutable struct for stuff in here
    # Check parameters

    #--- Geometry ---# (Only curved slab for now)
    gmsh.initialize(["", "-clmax", string(0.1)])
    model = gmsh.model


    # Bottom surface
    r = 1
    theta = pi/4
    phi = pi/8
    create_surface(model, r, theta, phi)




    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    gmsh.finalize()








end # End of create_brain_3D




if abspath(PROGRAM_FILE) == @__FILE__

    lc = 0.1       
    arcLen = 10    
    r_brain = 5    
    d_ratio = 0.5  
    r_curv = 20    
    BS_points = 50 
    params = model_params(lc, arcLen, r_brain, d_ratio, r_curv, BS_points)

    create_brain_3D(params)
end