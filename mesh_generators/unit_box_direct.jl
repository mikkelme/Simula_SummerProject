using Gridap
using GridapGmsh
using GridapGmsh: gmsh


include("./DiscreteModel_utils.jl")

function pgs_tags(pgs_dict, tags)
    return [pgs_dict[tag] for tag in tags]
end


"""
Unit square (sidelength = 1)
# Physical tags
   (8)--3--(7)
 4→ |      | ←2
   (5)--1--(6)
"""
function create_unit_box(lc, view=false)
    gmsh.initialize(["", "-clscale", string(lc)])

    A = gmsh.model.occ.addPoint(0, 0, 0)
    B = gmsh.model.occ.addPoint(1, 0, 0)
    C = gmsh.model.occ.addPoint(1, 1, 0)
    D = gmsh.model.occ.addPoint(0, 1, 0)

    P = [A, B, C, D]
    line = []

    for i in 1:4
        from, to = P[1+(i-1)%4], P[1+i%4]
        append!(line, gmsh.model.occ.addLine(from, to))

        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(1, [line[i]], i)
    end


    
    gmsh.model.addPhysicalGroup(0, [A], 5)
    gmsh.model.addPhysicalGroup(0, [B], 6)
    gmsh.model.addPhysicalGroup(0, [C], 7)
    gmsh.model.addPhysicalGroup(0, [D], 8)

    loop = gmsh.model.occ.addCurveLoop([line[1], line[2], line[3], line[4]])
    surf = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setTransfiniteSurface(surf)
    # gmsh.option.setNumber("Mesh.RecombineAll", 1) # For square mesh

    # surf_group = gmsh.model.addPhysicalGroup(2, [surf]) # This must be here for some reason when writing



    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end
    
    
    
    model, pgs_dict = direct_wiring(gmsh)
    gmsh.finalize()

    return model, pgs_dict
end



"""
Unit square (sidelength = 1) with dividing line from (0, 0.5) to (1, 0.5)
Physical tags
  (12)--7--(13)
 5→ |  200  | ←6
  (10)--4--(11)
 2→ |  100  | ←3
   (8)--1--(9)

"""
function create_coupled_box(lc, view=false)
    gmsh.initialize(["", "-clscale", string(lc)])




    A = gmsh.model.occ.addPoint(0, 0, 0)
    B = gmsh.model.occ.addPoint(1, 0, 0)
    C = gmsh.model.occ.addPoint(1, 1, 0)
    D = gmsh.model.occ.addPoint(0, 1, 0)
    ΓL = gmsh.model.occ.addPoint(0, 0.5, 0)
    ΓR = gmsh.model.occ.addPoint(1, 0.5, 0)

    P = Array{Int32,2}(undef, 2, 4)
    line = Array{Int32,2}(undef, 2, 4)
    surf = Array{Int32,1}(undef, 2)


    PD = [A, B, ΓR, ΓL] # Points Darcy
    PS = [ΓL, ΓR, C, D] # Points Stokes
    P[1, :] = PD
    P[2, :] = PS
    
    Γline = gmsh.model.occ.addLine(ΓL, ΓR)
   
    for i in 1:2
        for j in 1:4
            from, to = P[i, 1+(j-1)%4], P[i, 1+j%4]
            if (from, to) == (ΓL, ΓR)
                line[i,j] = Γline
            elseif (from, to) == (ΓR, ΓL)
                line[i,j] = -Γline

            else 
                line[i, j] = gmsh.model.occ.addLine(from, to)
            end   
        end
        loop = gmsh.model.occ.addCurveLoop([line[i, 1], line[i, 2], line[i, 3], line[i, 4]])     
        surf[i] = gmsh.model.occ.addPlaneSurface([loop])
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.setTransfiniteSurface(surf[i]) # Ordering to nice triangles og squares
    end

    # Physical group
    gmsh.model.addPhysicalGroup(1, [line[1, 1]], 1)
    gmsh.model.addPhysicalGroup(1, [line[1, 4]], 2)
    gmsh.model.addPhysicalGroup(1, [line[1, 2]], 3)
    gmsh.model.addPhysicalGroup(1, [Γline], 4)
    gmsh.model.addPhysicalGroup(1, [line[2, 4]], 5)
    gmsh.model.addPhysicalGroup(1, [line[2, 2]], 6)
    gmsh.model.addPhysicalGroup(1, [line[2, 3]], 7)

    gmsh.model.addPhysicalGroup(0, [A], 8)
    gmsh.model.addPhysicalGroup(0, [B], 9)
    gmsh.model.addPhysicalGroup(0, [ΓL],10)
    gmsh.model.addPhysicalGroup(0, [ΓR],11)
    gmsh.model.addPhysicalGroup(0, [D], 12)
    gmsh.model.addPhysicalGroup(0, [C], 13)

    gmsh.model.addPhysicalGroup(2, [surf[1]], 100)
    gmsh.model.addPhysicalGroup(2, [surf[2]], 200)





    gmsh.option.setNumber("Mesh.RecombineAll", 1) # For square mesh
    gmsh.model.mesh.generate(2)

    if view
        gmsh.fltk.initialize()
        gmsh.fltk.run()
    end


     
    model, pgs_dict = direct_wiring(gmsh)
    gmsh.finalize()

    return model, pgs_dict

end



# create_unit_box(2, true)
# create_coupled_box(2, true)

# if abspath(PROGRAM_FILE) == @__FILE__
#     lc = 0.1
#     lc = 2
# end



