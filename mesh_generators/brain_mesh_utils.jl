

mutable struct model_params
    # Parameters defining the brain model
    lc::Float64                       # Mesh size
    arcLen::Tuple{Float64,Float64}    # Outer arc length (x, y)-direction
    r_brain::Float64                  # Radial total length
    d_ratio::Float64                  # Radial relative length of fluid section
    r_curv::Float64                   # Radius of curvature
    inner_perturb::Function           # Radial perturbation on inner surface f(x,y,z)
    inner_perturb_body::String
    outer_perturb::Function           # Radial perturbation on inner surface f(x,y,z)
    outer_perturb_body::String
    BS_points::Tuple{Int64, Int64}    # Number of points in BSpline curves
    field_Lc_lim::Array{Float64,1}    # Field strengh as (LcMin, LcMax) multiplied by Lc when used
    field_Dist_lim::Array{Float64,1}  # field distance (DistMin, DistMax)

    function model_params(lc,arcLen, r_brain, d_ratio, r_curv, inner_perturb_body, outer_perturb_body, BS_points, field_Lc_lim, field_Dist_lim)
        new(lc,arcLen, r_brain, d_ratio, r_curv,((x,z) -> Base.invokelatest(eval(Meta.parse(inner_perturb_body)),x,z)), inner_perturb_body, ((x,z) -> Base.invokelatest(eval(Meta.parse(outer_perturb_body)),x,z)), outer_perturb_body, BS_points, field_Lc_lim, field_Dist_lim)
    end
end



function spherical_to_cartesian(r, theta, phi)
    # theta: angle from x-axis in x-z-plane (0 → 2π)
    # phi: angle from y-axis towards x-z-plane (0 → π)
    x = r * cos(theta) * sin(phi)
    z = r * sin(theta) * sin(phi)
    y = r * cos(phi)
    return x, y, z
end


function vecNorm(point)
    return sqrt(point[1]^2 + point[2]^2 + point[3]^2)
end



