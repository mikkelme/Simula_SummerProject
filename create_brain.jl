using GridapGmsh: gmsh
using Printf
    


function polar_to_cartesian(r, phi)
    x = r * cos(phi)
    y = r * sin(phi)
    return x, y
end

function add_manual_curve(model, r, theta)
    lc = 0.15 # Remember to get rid of HARDCODING HERE


    num_points = 5
    angle_space =  LinRange(0, 1, num_points)
    pointTags = [] # rethink how to tag them without risking duplicates 

    for i in range(1, num_points)
        x, y = polar_to_cartesian(r0, pi/2 + theta)
        append!(pointTags, model.geo.addPoint(x, y, 0.0))

    end
    Bspline = model.geo.addBSpline(pointTags)

    # left_tag = pointTags[1]
    # right_tag = last(pointTags)
   
    println(pointTags)
    exit()
    # xL, yL = polar_to_cartesian(r1 - d, pi/2 + theta)
    # xR, yR = polar_to_cartesian(r1 - d, pi/2 - theta)
    # model.geo.addPoint(xL, yL, 0.0, lc, 3) # Left
    # model.geo.addPoint(xR, yR, 0.0, lc, 4) # Right
    # model.geo.addCircleArc(left_tag, 0, right_tag, curve_tag) 
    # exit()


    # xL, yL = polar_to_cartesian(r1 - d, pi/2 + theta)
    # xR, yR = polar_to_cartesian(r1 - d, pi/2 - theta)
    # L = model.geo.addPoint(xL, yL, 0.0, lc) # Left
    # R = model.geo.addPoint(xR, yR, 0.0, lc) # Right
    # arc = model.geo.addCircleArc(L, C, R) 
    # return left_tag, right_tag
    # model.geo.addCircleArc(3, 0, 4, curve_tag) 

    
 



end


function create_brain(r0, r1, d, theta, view = true)

    if r1-r0 <= d
        @printf("r1 - r0 = %g must be greater than d = %g\n", r1-r0, d)
        exit()
    end

    gmsh.initialize()
    model = gmsh.model
    # factory = model.occ

    lc = 0.15

    C = model.geo.addPoint(0.0, 0.0, 0.0, lc) # Center


    #Inner arc
    xL, yL = polar_to_cartesian(r0, pi/2 + theta)
    xR, yR = polar_to_cartesian(r0, pi/2 - theta)
    IL = model.geo.addPoint(xL, yL, 0.0, lc) # Inner Left
    IR = model.geo.addPoint(xR, yR, 0.0, lc) # Inner Right
    I_arc = model.geo.addCircleArc(IL, C, IR) # Inner arc

    

    # Dividing arc

    # xL, yL = polar_to_cartesian(r1 - d, pi/2 + theta)
    # xR, yR = polar_to_cartesian(r1 - d, pi/2 - theta)
    # model.geo.addPoint(xL, yL, 0.0, lc, 3) # Left
    # model.geo.addPoint(xR, yR, 0.0, lc, 4) # Right
    # model.geo.addCircleArc(3, 0, 4, curve_tag) 
   
    # xL, yL = polar_to_cartesian(r1 - d, pi/2 + theta)
    # xR, yR = polar_to_cartesian(r1 - d, pi/2 - theta)
    # L = model.geo.addPoint(xL, yL, 0.0, lc) # Left
    # R = model.geo.addPoint(xR, yR, 0.0, lc) # Right
    # arc = model.geo.addCircleArc(L, C, R) 
    add_manual_curve(model, r1-d, theta)


    # Outer arc
    xL, yL = polar_to_cartesian(r1, pi/2 + theta)
    xR, yR = polar_to_cartesian(r1, pi/2 - theta)
    OL = model.geo.addPoint(xL, yL, 0.0, lc) # Outer Left
    OR = model.geo.addPoint(xR, yR, 0.0, lc) # Outer Right
    O_arc = model.geo.addCircleArc(OL, C, OR)  # Outer arc



    # # Inner tisue
    # model.geo.addLine(1, left_tag, right_tag) 
    # model.geo.addLine(2, right_tag, 5) 
    # model.geo.addCurveLoop([4, 2, -5, -1], 1)
    # model.geo.addPlaneSurface([1], 1)
    
    
    
    # # Outer fluid
    # model.geo.addLine(left_tag, 5, 6) 
    # model.geo.addLine(right_tag, 6, 7) 
    # model.geo.addCurveLoop([6,3,-7, -2], 2)
    # model.geo.addPlaneSurface([2], 2)




    model.geo.synchronize()
    # model.mesh.generate(2) 
    
    
    
    # exit()

  

     if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end
    
    
    gmsh.finalize()



end




r0 = 7
r1 = 10
d = 1
angle = pi/10
create_brain(r0, r1, d, angle)