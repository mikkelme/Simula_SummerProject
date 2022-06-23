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


function create_surface(model, r, zx_angle,  zy_angle)

    origo = model.geo.addPoint(0.0, 0.0, 0.0) # Center

    # Make boundary of surface
    ###################
    # alpha = (1 - z) / y
    # Use inclined circle (for y incline)
    # eq.: x^2 + (z + αy)^2 = r^2
    ###################


    x, _, zx = spherical_to_cartesian(r, 0, zx_angle / 2)
    _, y, zy = spherical_to_cartesian(r, pi / 2, zy_angle / 2)

    alpha_x = (1 - zx) / x
    alpha_y = (1 - zy) / y


    a = 2 * alpha_x / (1 - alpha_x^2)
    b = 2 * alpha_y / (1 - alpha_y^2)
    c = r / sqrt(a^2 + b^2 + 1)


    A = model.geo.addPoint(a * c, b * c, c)
    B = model.geo.addPoint(-a * c, b * c, c)
    C = model.geo.addPoint(-a * c, -b * c, c)
    D = model.geo.addPoint(a * c, -b * c, c)


    AB_arc = model.geo.addCircleArc(A, origo, B)
    BC_arc = model.geo.addCircleArc(B, origo, C)
    CD_arc = model.geo.addCircleArc(C, origo, D)
    DA_arc = model.geo.addCircleArc(D, origo, A)



 


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




    # model.occ.synchronize()
    model.geo.synchronize()

 




end

function create_brain_3D(params::model_params, view=true)
    # Consider mutable struct for stuff in here
    # Check parameters

    #--- Geometry ---# (Only curved slab for now)
    gmsh.initialize(["", "-clmax", string(0.1)])
    model = gmsh.model


    # Bottom surface
    r = 1
    # 0 < angle < pi
    zx_angle = pi / 2
    zy_angle = pi / 3
    create_surface(model, r, zx_angle, zy_angle)




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