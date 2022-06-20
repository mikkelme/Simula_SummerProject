using GridapGmsh: gmsh



function create_geometry(view = true)
    """
    Goal: Make a rectangle on top of square box 
       ____
      [_  _]
        []
    """
 
    
    
    gmsh.initialize()
    model = gmsh.model
    factory = model.occ

    # upperbox = factory.addBox(0, 0, 0, 1, 1, 1) #x, y, z, dx, dy, dz, tag = -1
    # lowerbox = factory.addBox(0.25, 0.25, 0.25, 0.5, 0.5, 0.5)

    # factory.synchronize()


    lc = 0.15

    # Upper box
    model.geo.addPoint(0.0, 0.0, 0.0, lc, 1) 
    model.geo.addPoint(2.0, 0.0, 0.0, lc, 2)
    model.geo.addPoint(2.0, 1.0, 0.0, lc, 3)
    model.geo.addPoint(0, 1.0, 0.0, lc, 4)
    model.geo.addPoint(1, 0.5, 0.0, lc, 5) # For setting surface

    model.geo.addLine(1,2,1) 
    model.geo.addLine(2,3,2)
    model.geo.addLine(3,4,3)
    model.geo.addLine(4,1,4)

    model.geo.addCurveLoop([1,2,3,4], 5)
    model.geo.addPlaneSurface([5], 6)

    model.addPhysicalGroup(1, [5], 6)
    model.setPhysicalName(1, 6, "Upper box")
    # model.geo.synchronize()
    
    # Lower box
    model.geo.addPoint(0.5, 0.0, 0.0, lc, 11) 
    model.geo.addPoint(1.5, 0.0, 0.0, lc, 12)
    model.geo.addPoint(1.5, -1.0, 0.0, lc, 13)
    model.geo.addPoint(0.5, -1.0, 0.0, lc, 14)
    model.geo.addPoint(1, -0.5, 0.0, lc, 15)

    model.geo.addLine(11,12,11) 
    model.geo.addLine(12,13,12)
    model.geo.addLine(13,14,13)
    model.geo.addLine(14,11,14)

    model.geo.addCurveLoop([11,12,13,14], 15)
    model.geo.addPlaneSurface([15], 16)

    model.addPhysicalGroup(2, [15], 16)
    model.setPhysicalName(2, 16, "Lower box")
    model.geo.synchronize()

    model.mesh.generate(2) # dim
    
    # model.addPhysicalGroup(0, [1, 2], 1) # dim, tags, tag
    # model.addPhysicalGroup(1, [1, 2], 2)
    # model.addPhysicalGroup(2, [1], 6)
    # model.setPhysicalName(2, 6, "My surface")
    
    
    # upperbox = factory.addBox(0, 0, 0, 1, 1, 1)
    # lowerbox = factory.addBox(0.25, 0.25, 0.25, 0.5, 0.5, 0.5)
    # factory.synchronize()







    # lc = .15
    # gmsh.model.geo.addPoint(0.0, 0.0, 0, lc, 1) #x, y, z, meshsize = 0, tag = -1
    # gmsh.model.geo.addPoint(1, 0.0, 0, lc, 2)
    # gmsh.model.geo.addPoint(1, 1, 0, lc, 3)
    # gmsh.model.geo.addPoint(0, 1, 0, lc, 4)
    # gmsh.model.geo.addPoint(0.2, .5, 0, lc, 5)

    # gmsh.model.geo.addLine(1, 2, 1) #start tag, endtag, tag = -1
    # gmsh.model.geo.addLine(2, 3, 2)
    # gmsh.model.geo.addLine(3, 4, 3)
    # gmsh.model.geo.addLine(4, 1, 4)

    # gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 5) # Go around the lines (connect)
    # gmsh.model.geo.addPlaneSurface([5], 6)

    # gmsh.model.geo.synchronize() # to gather what you have done I think



    
    
    
    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end
    
    
    gmsh.finalize()



end





# create_geometry()










