using GridapGmsh: gmsh
using Printf
    


function polar_to_cartesian(r, phi)
    x = r * cos(phi)
    y = r * sin(phi)
    return x, y
end


function add_manual_curve(model, r, theta)
    # lc = 0.15 # Remember to get rid of HARDCODING HERE


    num_points = 10
    angle_space =  LinRange(0, 1, num_points)
    pointTags = [] # rethink how to tag them without risking duplicates 

    for i in range(1, num_points)
        angle = pi/2 + (1 - 2*(i-1)/(num_points-1)) * theta/2 #angle from +theta/2 -> -theta/2
        # angle = pi/2 +  (num_points + 1 - 2*i)/(num_points - 1) * theta/2

        x, y = polar_to_cartesian(r, angle)
        append!(pointTags, model.geo.addPoint(x, y, 0.0))

    end
    BSpline = model.geo.addBSpline(pointTags)

    return pointTags[1], last(pointTags), BSpline


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
    xL, yL = polar_to_cartesian(r0, pi/2 + theta/2)
    xR, yR = polar_to_cartesian(r0, pi/2 - theta/2)
    IL = model.geo.addPoint(xL, yL, 0.0, lc) # Inner Left
    IR = model.geo.addPoint(xR, yR, 0.0, lc) # Inner Right
    I_arc = model.geo.addCircleArc(IL, C, IR) # Inner arc


    # Dividing arc
    DL, DR, D_arc = add_manual_curve(model, r1-d, theta)
 

    # Outer arc
    OL, OR, O_arc = add_manual_curve(model, r1, theta)


    # Inner tisue
    IL_line = model.geo.addLine(IL, DL) 
    IR_line = model.geo.addLine(IR, DR) 
    I_CurveLoop = model.geo.addCurveLoop([IL_line, D_arc, -IR_line, -I_arc])
    I_surf = model.geo.addPlaneSurface([I_CurveLoop])
    
    
    
    # Outer fluid
    OL_line = model.geo.addLine(OL, DL) 
    OR_line = model.geo.addLine(OR, DR) 
    O_CurveLoop = model.geo.addCurveLoop([OL_line, D_arc, -OR_line, -O_arc])
    O_surf = model.geo.addPlaneSurface([O_CurveLoop])
    
  

    model.geo.synchronize()
    model.mesh.generate(2) 
    
    
  
     if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end
    
    
    gmsh.finalize()



end




r0 = 7
r1 = 10
d = 1
angle = pi/6
create_brain(r0, r1, d, angle)