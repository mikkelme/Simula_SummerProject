using GridapGmsh: gmsh
using Printf
    


function polar_to_cartesian(r, phi)
    x = r * cos(phi)
    y = r * sin(phi)
    return x, y
end


function add_manual_curve(model, r, theta, perturbation_func, lc, num_points = 50)
    """
    Description

    Model: 
    r:
    theta:
    perturbation_func: Should be 0 at endpoints x = 0,1
    num_points

    return
    """

    # num_points = 2
    angle_space =  LinRange(0, 1, num_points)
    pointTags = [] 

    for i in range(1, num_points)
        angle = pi/2 + theta/2 - theta * angle_space[i]
        x, y = polar_to_cartesian(r + perturbation_func(angle_space[i]), angle)
        append!(pointTags, model.geo.addPoint(x, y, 0.0, lc))

    end
    BSpline = model.geo.addBSpline(pointTags)

    return pointTags[1], last(pointTags), BSpline


end


function create_brain(arcLen, r_brain, d_ratio, r_curv, D_func, O_func, view = true)

    # Compute main geometry parameters
    rI = r_curv - r_brain # > 0    
    d = d_ratio * r_brain
    rD = r_curv - d
    theta = arcLen/r_curv # > 0 & < pi
    
    # Check parameters
    illegal_param = false
    tol = 1e-10
 
    if d_ratio <= 0 || d_ratio >= 1 
        @printf("d_ratio = %g must be in the open interval (0,1). \n", d_ratio)
        illegal_param = true
    end
    
    if rI <= 0
        @printf("Inner radius = %g must be positive. Choose r_curv > r_brain. \n", rI)
        illegal_param = true
    end

    if theta <= 0 || theta >= pi 
        @printf("theta = %g must be in the open interval (0,pi). Choose arcLen and r_curv differently. \n", theta)
        illegal_param = true
    end

    if abs(D_func(0)) > tol || abs(D_func(1)) > tol
        @printf("D_func(0) = %g and D_func(1) = %g should be 0. \n", D_func(0), D_func(1))
        illegal_param = true
    end

    if abs(O_func(0)) > tol || abs(O_func(1)) > tol
        @printf("O_func(0) = %g and O_func(1) = %g should be 0. \n", O_func(0), O_func(1))
        illegal_param = true
    end

    illegal_param ? exit(0) :
    


    #--- Meshing ---#

    gmsh.initialize()
    model = gmsh.model
    # factory = model.occ

    lc = 0.15

    C = model.geo.addPoint(0.0, 0.0, 0.0, lc) # Center


    #Inner arc
    xL, yL = polar_to_cartesian(rI, pi/2 + theta/2)
    xR, yR = polar_to_cartesian(rI, pi/2 - theta/2)
    IL = model.geo.addPoint(xL, yL, 0.0, lc) # Inner Left
    IR = model.geo.addPoint(xR, yR, 0.0, lc) # Inner Right
    I_arc = model.geo.addCircleArc(IL, C, IR) # Inner arc


    # Dividing arc
    DL, DR, D_arc = add_manual_curve(model, rD, theta, D_func, lc)
 

    # Outer arc
    OL, OR, O_arc = add_manual_curve(model, r_curv, theta, O_func, lc)


    # Inner tisue
    IL_line = model.geo.addLine(IL, DL) 
    IR_line = model.geo.addLine(IR, DR) 
    I_CurveLoop = model.geo.addCurveLoop([IL_line, D_arc, -IR_line, -I_arc])
    I_surf = model.geo.addPlaneSurface([I_CurveLoop])

    I_PGroup = model.addPhysicalGroup(2, [I_surf])
    model.setPhysicalName(2, I_PGroup, "Brain tissue")
    
    
    
    # Outer fluid
    OL_line = model.geo.addLine(OL, DL) 
    OR_line = model.geo.addLine(OR, DR) 
    O_CurveLoop = model.geo.addCurveLoop([OL_line, D_arc, -OR_line, -O_arc])
    O_surf = model.geo.addPlaneSurface([O_CurveLoop])

    O_PGroup = model.addPhysicalGroup(2, [O_surf])
    model.setPhysicalName(2, O_PGroup, "Brain fluid")
    
  
    

    model.geo.synchronize()
    model.mesh.generate(2) 
    
     if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end
    
    gmsh.finalize()



end



arcLen = 10     # Outer arc length 
r_brain = 5    # radial length of computed area
d_ratio = 0.2  # thickness of fluid section relative to to r_brain
r_curv = 20    # Radius of curvature

D_func(x) = 0.2*sin(pi*x*20)
O_func(x) = 0.2*sin(pi*x*2)
create_brain(arcLen, r_brain, d_ratio, r_curv, D_func, O_func)