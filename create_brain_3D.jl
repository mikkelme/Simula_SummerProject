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

function create_brain_3D(params::model_params)
    # Consider mutable struct for stuff in here
    # Check parameters

    #--- Geometry ---# (Only curved slab for now)
    gmsh.initialize(["", "-clmax", string(0.1)])
    model = gmsh.model

    
    


end # End of create_brain_3D


lc = 0.1       # Mesh size (definition?)
arcLen = 10    # Outer arc length 
r_brain = 5    # radial length of computed area
d_ratio = 0.5  # thickness of fluid section relative to to r_brain
r_curv = 20    # Radius of curvature
BS_points = 50 # Number of points in BSpline curves

params = model_params(lc, arcLen, r_brain, d_ratio, r_curv, BS_points)
create_brain_3D(params)