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
            println(dis)
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




function create_brain_2D(arcLen, r_brain, d_ratio, r_curv, D_func, O_func, BS_points, view=true)

    # Check for valid parameters
    check_parameters(arcLen, r_brain, d_ratio, r_curv, D_func, O_func, BS_points)

    
    #--- Meshing ---#
    gmsh.initialize(["", "-clmax", "0.1"])
    model = gmsh.model

    if r_curv == 0  # Straight Box
        L = arcLen              # Length
        h = r_brain             # Height
        Dh = (1 - d_ratio) * h  # Dividing height

        # Box corners
        BL = model.geo.addPoint(0.0, 0.0, 0.0)      # Bottom Left
        BR = model.geo.addPoint(arcLen, 0.0, 0.0)   # Bottom Right

        # Inner tisue
        DL, DR, D_line = add_perturbed_line(model, L, Dh, D_func, BS_points) #Dividing Left, Dividing Right, Inner Top
        IL = model.geo.addLine(BL, DL) # Left line 
        IR = model.geo.addLine(DR, BR) # Right line
        IB = model.geo.addLine(BR, BL) # Bottom line
        I_CurveLoop = model.geo.addCurveLoop([IL, D_line, IR, IB])
        I_surf = model.geo.addPlaneSurface([I_CurveLoop])

        # Outer fluid
        TL, TR, O_line = add_perturbed_line(model, L, h, O_func, BS_points) # Top left, Top right, Outer Top
        OL = model.geo.addLine(DL, TL) # Left line 
        OR = model.geo.addLine(TR, DR) # Right line
        O_CurveLoop = model.geo.addCurveLoop([OL, O_line, OR, -D_line])
        O_surf = model.geo.addPlaneSurface([O_CurveLoop])



    else # Curved slap 
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

    model.geo.synchronize()

    # Physics gorup
    I_PGroup = model.addPhysicalGroup(2, [I_surf])
    model.setPhysicalName(2, I_PGroup, "Brain tissue")

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



arcLen = 10    # Outer arc length 
r_brain = 5    # radial length of computed area
d_ratio = 0.2  # thickness of fluid section relative to to r_brain
r_curv = 50    # Radius of curvature
BS_points = 50 # Number of points in BSpline curves


D_func(x) = 0.2*sin(pi*x/0.1)
O_func(x) = 0.2*sin(pi*x/0.5)
create_brain_2D(arcLen, r_brain, d_ratio, r_curv, D_func, O_func, BS_points)




# Things to do:
# √ add num_points
# √ option for flat section
# periodic mesh (equal in left and right)
# √ take care of overlap between D_func and O_func
