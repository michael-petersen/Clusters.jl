

using DelimitedFiles # ability to read the tab-delimited text files
using Statistics     # access to mean

function return_density_centre(filename::String)
    # read a simple centering file
    allpos = readdlm(filename, Float32,skipstart=1)
    
    x  = allpos[:,1]
    y  = allpos[:,2]
    z  = allpos[:,3]    

    return x,y,z
end



function return_particles(filename::String)
    # function to read in ascii versions of nbody outputs (x y z vx vy vz)
    allpos = readdlm(filename, Float32,skipstart=1)

    x  = allpos[:,1]
    y  = allpos[:,2]
    z  = allpos[:,3]
    vx = allpos[:,4]
    vy = allpos[:,5]
    vz = allpos[:,6]

    return x,y,z,vx,vy,vz
end


function find_rbary(x::Vector{Float32},y::Vector{Float32},z::Vector{Float32},vx::Vector{Float32},vy::Vector{Float32},vz::Vector{Float32})

    # function to compute the barycentric radius of particles

    # these are the barycentre values
    xmean  = mean(x)
    ymean  = mean(y)
    zmean  = mean(z)
    vxmean = mean(vx)
    vymean = mean(vy)
    vzmean = mean(vz)

    # rbary = barycentre radius
    rbary = ((x.-xmean).^2 + (y.-ymean).^2 + (z.-zmean).^2).^0.5

    # sort these to look at the largest particle?
    return rbary
end





function perturbation_density(xvals::Vector{Float32},
                              yvals::Vector{Float32},
                              zvals::Vector{Float32},
                              radius::Float64,
                              mass::Float64=1.0)
    #
    # estimating the perturbation density
    # formally written as (eq 1)
    #   \rho_1(r,\theta,\phi) = a_1(r)\sin\theta\cos\phi + a_2(r)\sin\theta\sin\phi + a_3(r)\cos\theta
    # but we can do an estimate of the prefactors a_{1,2,3} using summation
    # choosing some radius to truncate the particles
    # (eq 3)
    
    phi   = atan.(yvals,xvals)
    theta = atan.(sqrt.(xvals.^2+yvals.^2)./zvals)
    
    volume = (4/3)*pi*radius^3
    
    
    a1 = mass.*(3/volume)*sum(sin.(theta).*cos.(phi))
    a2 = mass.*(3/volume)*sum(sin.(theta).*sin.(phi))
    a3 = mass.*(3/volume)*sum(cos.(theta))
    
    return a1,a2,a3
    
end



