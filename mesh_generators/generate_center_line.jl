using Gridap
using GridapGmsh
using GridapGmsh: gmsh
include("./brain_mesh_utils.jl")
include("./DiscreteModel_utils.jl")



# function get_points(arcLen, r_brain, d_ratio, r_curv, perturbation_func, BS_points)
#     points = []

   
#     r_cl = r_curv - (d_ratio * r_brain)/2 
#     angle = arcLen / r_curv   
  

#     angle = arcLen / r_curv  
#     for i in LinRange(0, 1, BS_points)
#         phi = - angle / 2 + angle * i
#         x_arcLen = r_cl * phi
#         r_pert = r_cl + perturbation_func(x_arcLen, 0)
#         append!(points, [spherical_to_cartesian(r_pert, 0, phi)])
#     end

#     return points
# end




function create_centerline(arcLen, r_brain, d_ratio, r_curv, perturbation_func, BS_points)
    pointTags = []
   
    r_cl = r_curv - (d_ratio * r_brain)/2 
    angle = arcLen / r_curv   
  

    angle = arcLen / r_curv  
    for i in LinRange(0, 1, BS_points)
        phi = - angle / 2 + angle * i
        x_arcLen = r_cl * phi
        r_pert = r_cl + perturbation_func(x_arcLen, 0)
        append!(pointTags, gmsh.model.occ.addPoint(spherical_to_cartesian(r_pert, 0, phi)...))
    end

    vertex = [pointTags[1], last(pointTags)]
    BSpline = gmsh.model.occ.addBSpline(pointTags)
    gmsh.model.occ.synchronize()

    return vertex, BSpline
end


pert_func(x,z) = 0.0
# points = get_points(5, 2, 0.5, 50, pert_func, 10)
lc = 0.1

gmsh.initialize(["", "-clmax", string(lc)])
gmsh.option.setNumber("Mesh.SaveAll", 1)  # For direct wiring

vertex, BSpline = create_centerline(5, 2, 0.5, 50, pert_func, 3)
gmsh.model.addPhysicalGroup(0, [vertex[1]]) # center line (left point)
gmsh.model.addPhysicalGroup(1, [BSpline])   # center line (arc)
gmsh.model.addPhysicalGroup(0, [vertex[2]]) # center line (right point)


gmsh.model.mesh.generate(1)


if true
    gmsh.fltk.initialize()
    gmsh.fltk.run()
end

model, pgs_dict = direct_wiring(gmsh)
gmsh.finalize()






# gmsh.initialize(["", "-clscale", string(0.5)])
# A = gmsh.model.occ.addPoint(points[1]...)
# B = gmsh.model.occ.addPoint(points[2]...)
# C = gmsh.model.occ.addPoint(points[3]...)

# AB = gmsh.model.occ.addLine(A, B)
# BC = gmsh.model.occ.addLine(B, C)

# gmsh.model.occ.synchronize()
# gmsh.model.mesh.generate(1)


# if false
#     gmsh.fltk.initialize()
#     gmsh.fltk.run()
# end


# model, pgs_dict = direct_wiring(gmsh)
# gmsh.finalize()












# model = CartesianDiscreteModel((0, 1), (3, ))





# desc = Gridap.Geometry.CartesianDescriptor()

# CartesianGrid(desc)


# node_cords = Gridap.Geometry.CartesianCoordinates{1, Float64, typeof(identity)
# node_cords[1] = VectorValue{1, Float64}



# println(node_cords)

# desc = Gridap.Geometry.CartesianDescriptor{D, T, F}



# Gridap.Geometry.CartesianCoordinates{1, Float64, typeof(identity)}:
#                 VectorValue{1, Float64}(0.0,)
#  VectorValue{1, Float64}(0.3333333333333333,)
#  VectorValue{1, Float64}(0.6666666666666666,)
#                 VectorValue{1, Float64}(1.0,)

# CartesianGrid(desc::CartesianDescriptor{node_coords,cell_node_ids,cell_type})



# model = CartesianDiscreteModel((0, 1), (3, ))



# --- Goal --- #
# grid = UnstructuredGrid(node_to_coords, cell_to_nodes, reffes, cell_to_type, orientation, facet_normal)

# UnstructuredDiscreteModel(grid,grid_topology,labeling)