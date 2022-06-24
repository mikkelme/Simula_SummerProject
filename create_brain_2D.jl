using GridapGmsh: gmsh
using Printf
    


function polar_to_cartesian(r, phi)
    x = r * cos(phi)
    y = r * sin(phi)
    return x, y
end



"""
Description

Model: 
r:
theta:
perturbation_func: Should be 0 at endpoints x = 0,1
num_points

return
"""
function add_perturbed_arc(model, r, theta, perturbation_func, BS_points)
    pointTags = [] 

    for i in LinRange(0, 1, BS_points)
        angle = pi/2 + theta/2 - theta * i
        x, y = polar_to_cartesian(r + perturbation_func(i), angle)
        append!(pointTags, model.geo.addPoint(x, y, 0.0))
    end

    BSpline = model.geo.addBSpline(pointTags)
    return pointTags[1], last(pointTags), BSpline
end



function add_perturbed_line(model, L, Dh, perturbation_func, BS_points)
    input_space = LinRange(0, 1, BS_points)
    pointTags = []

    for i in LinRange(0, 1, BS_points)
        x = i * L
        y = Dh + perturbation_func(i)
        append!(pointTags, model.geo.addPoint(x, y, 0.0))

    end

    BSpline = model.geo.addBSpline(pointTags)
    return pointTags[1], last(pointTags), BSpline
end




function check_parameters(arcLen, r_brain, d_ratio, r_curv, D_func, O_func, BS_points)

    illegal_param = false
    tol = 1e-10

    if d_ratio <= 0 || d_ratio >= 1
        @printf("d_ratio = %g must be in the open interval (0,1). \n", d_ratio)
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

    if BS_points < 2
        @printf("BS_points = %g should be >= 2. \n", BS_points)
        illegal_param = true
    end


    if r_curv == 0 # Box
        L = arcLen              # Length
        h = r_brain             # Height
        Dh = (1 - d_ratio) * h    # Dividing height
    
    
        if L <= 0
            printf("L (arcLen) = %g must be > 0.\n", h)
            illegal_param = true
        end
    
        if h <= 0
            printf("h (r_brain) = %g must be > 0.\n", h)
            illegal_param = true
        end
    
    
        for i in LinRange(0, 1, BS_points)
            dis = ((h + O_func(i)) - (Dh + D_func(i)))
            if dis < 0
                @printf("BSpline points (outer and dividing line) overlap by a distance %g at x = %g\n", dis, i * L)
                illegal_param = true
                break
            end
        end
    
    
    
    else # Curved slab
        rI = r_curv - r_brain
        d = d_ratio * r_brain
        rD = r_curv - d
        theta = arcLen / r_curv
    
        if rI <= 0
            @printf("Inner radius = %g must be positive. Choose r_curv > r_brain. \n", rI)
            illegal_param = true
        end
    
        if theta <= 0 || theta >= pi
            @printf("theta = %g must be in the open interval (0,pi). Choose arcLen and r_curv differently. \n", theta)
            illegal_param = true
        end
    
 
        for i in LinRange(0, 1, BS_points)
            angle = pi / 2 + theta / 2 - theta * i
            dis = ((r_curv + O_func(i)) - (rD + D_func(i)))
            if dis < 0
                @printf("BSpline points (outer and dividing curve) overlap by a distance %g intersect at angle = %g\n", dis, angle)
                illegal_param = true
                break
            end
        end
    
    
    end
    
    if illegal_param exit() end
    

end # End of function




function create_brain_2D(arcLen, r_brain, d_ratio, r_curv, D_func, O_func, BS_points, lc, view=true)

    # Check for valid parameters
    check_parameters(arcLen, r_brain, d_ratio, r_curv, D_func, O_func, BS_points)


    #--- Geometry ---#
    gmsh.initialize(["", "-clmax", string(lc)])
    model = gmsh.model

    
    if r_curv == 0  # Straight Box
        box_mode = true
        L = arcLen              # Length
        h = r_brain             # Height
        Dh = (1 - d_ratio) * h  # Dividing height
    
        # Box corners
        BL = model.geo.addPoint(0.0, 0.0, 0.0)      # Bottom Left
        BR = model.geo.addPoint(arcLen, 0.0, 0.0)   # Bottom Right
    
        # Inner tisue
        DL, DR, D_line = add_perturbed_line(model, L, Dh, D_func, BS_points) #Dividing Left, Dividing Right, Inner Top
        IL_line = model.geo.addLine(BL, DL) # Left line 
        IR_line = model.geo.addLine(DR, BR) # Right line
        IB = model.geo.addLine(BR, BL) # Bottom line
        I_CurveLoop = model.geo.addCurveLoop([IL_line, D_line, IR_line, IB])
        I_surf = model.geo.addPlaneSurface([I_CurveLoop])
    
        # Outer fluid
        TL, TR, O_line = add_perturbed_line(model, L, h, O_func, BS_points) # Top left, Top right, Outer Top
        OL_line = model.geo.addLine(DL, TL) # Left line 
        OR_line = model.geo.addLine(TR, DR) # Right line
        O_CurveLoop = model.geo.addCurveLoop([OL_line, O_line, OR_line, -D_line])
        O_surf = model.geo.addPlaneSurface([O_CurveLoop])
    
    
        # For generalisation in the following
        I_arc = IB
        D_arc = D_line
        O_arc = O_line
    
    
    
    
    else # Curved slap 
        box_mode = false
        rI = r_curv - r_brain   # Inner radius  
        d = d_ratio * r_brain   # Thickness of fluid
        rD = r_curv - d         # Radius for dividing line
        theta = arcLen / r_curv # Spanning angle 
    
    
        C = model.geo.addPoint(0.0, 0.0, 0.0) # Center
    
        #Inner arc
        xL, yL = polar_to_cartesian(rI, pi / 2 + theta / 2)
        xR, yR = polar_to_cartesian(rI, pi / 2 - theta / 2)
        IL = model.geo.addPoint(xL, yL, 0.0) # Inner Left
        IR = model.geo.addPoint(xR, yR, 0.0) # Inner Right
        I_arc = model.geo.addCircleArc(IL, C, IR) # Inner arc
    
    
        # Dividing arc
        DL, DR, D_arc = add_perturbed_arc(model, rD, theta, D_func, BS_points)
    
        # Outer arc
        OL, OR, O_arc = add_perturbed_arc(model, r_curv, theta, O_func, BS_points)
    
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
    
    
    end





    #--- Mesh Field: Distance to dividing BSpline ---#
    LcMin = lc / 2  # Hardcoded for now
    LcMax = lc      # Hardcoded for now
    DistMin = 0.3
    DistMax = 0.5

    # LcMin -                       /------------------
    #                              /
    #                             /
    #                            /
    # LcMax   -o----------------/
    #          |                |    |
    #       BSpline        DistMin  DistMax

    F_distance = model.mesh.field.add("Distance")
    model.mesh.field.setNumber(F_distance, "NNodesByEdge", 100) # Try 2 if complaining
    model.mesh.field.setNumbers(F_distance, "EdgesList", [D_arc])
    # model.mesh.field.setNumbers(F_distance, "EdgesList", [IL_line]) # To verify periodic boundary mesh



    F_threshold = model.mesh.field.add("Threshold")
    model.mesh.field.setNumber(F_threshold, "IField", F_distance)
    model.mesh.field.setNumber(F_threshold, "LcMin", LcMin)
    model.mesh.field.setNumber(F_threshold, "LcMax", LcMax)
    model.mesh.field.setNumber(F_threshold, "DistMin", DistMin)
    model.mesh.field.setNumber(F_threshold, "DistMax", DistMax)

    gmsh.model.mesh.field.setAsBackgroundMesh(F_threshold)

    # The following should prevent over-refinement due to small mesh sizes on the boundary.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    model.geo.synchronize()


    #--- Periodic meshing ---#   
    if box_mode
        # Translation along x-direction
        translation = [ 1, 0, 0, arcLen,
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, 1 ]
         affineTransform = translation
    else 
        # Rotation around z-axis (in negative direction)
        rotation = [    cos(-theta) , -sin(-theta)  , 0, 0,
                        sin(-theta) , cos(-theta)   , 0, 0,
                        0           , 0             , 1, 0,
                        0           , 0             , 0, 1 ]
        affineTransform = rotation 
    end

    # Impose mesh from left boundary on right boundar 
    model.mesh.setPeriodic(1, [IR_line], [IL_line], affineTransform)  # Inner part
    model.mesh.setPeriodic(1, [OR_line], [OL_line], affineTransform)  # Outer part




    #--- Physics groups ---#


    I_arc_face = model.addPhysicalGroup(1, [I_arc])
    IL_line_face = model.addPhysicalGroup(1, [IL_line])
    IR_line_face = model.addPhysicalGroup(1, [IR_line])
    D_arc_face = model.addPhysicalGroup(1, [D_arc])
    OL_line_face = model.addPhysicalGroup(1, [OL_line])
    OR_line_face = model.addPhysicalGroup(1, [OR_line])
    O_arc_face = model.addPhysicalGroup(1, [O_arc])

    I_surf_face = model.addPhysicalGroup(2, [I_surf])
    model.setPhysicalName(2, I_surf_face, "Brain tissue")

    O_surf_face = model.addPhysicalGroup(2, [O_surf])
    model.setPhysicalName(2, O_surf_face, "Brain fluid")

    model.geo.synchronize()


    #--- Finish mesh and view ---#
    model.mesh.generate(2)

    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end

    gmsh.finalize()



end


lc = 0.1       # Mesh size (definition?)
arcLen = 10    # Outer arc length 
r_brain = 5    # radial length of computed area
d_ratio = 0.5  # thickness of fluid section relative to to r_brain
r_curv = 20    # Radius of curvature
BS_points = 50 # Number of points in BSpline curves

arcLen = 10
# r_curv = 0

D_func(x) = 0.2*sin(pi*x/0.1)
O_func(x) = 0.2*sin(pi*x/0.5)
create_brain_2D(arcLen, r_brain, d_ratio, r_curv, D_func, O_func, BS_points, lc)




# Things to do:
    # √ add num_points
    # √ option for flat section
# periodic mesh (equal in left and right)
    # √ take care of overlap between D_func and O_func
    # √ Work with fields
    # √ Define physics groups on all facets