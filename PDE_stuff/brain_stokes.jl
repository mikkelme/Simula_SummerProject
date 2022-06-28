using Gridap
using GridapGmsh



model = GmshDiscreteModel("brain.msh")
# writevtk(model, "model")

dividing_arc = 2
lower_boundary_tags = [1, 4, 5]
upper_boundary_tags = [3, 7, 8]



order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
# V = TestFESpace(model, reffe, dirichlet_tags=["side"])
V = TestFESpace(model, reffe, dirichlet_tags=[1])




exit()
# Define Dirichlet boundaries
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "diri1", [6,])
add_tag_from_tags!(labels, "diri0", [1, 2, 3, 4, 5, 7, 8])
