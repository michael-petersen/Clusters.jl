

using DelimitedFiles # ability to read the tab-delimited text files
using Statistics     # access to mean

function return_density_centre(filename::String)
    # read a simple centering file
    allpos = readdlm(filename, Float32,skipstart=0)
    
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


function find_rbary(x::Vector{Float32},y::Vector{Float32},z::Vector{Float32})

    # function to compute the barycentric radius of particles
    # this could use functionality to also take a centre, not just
    # compute the mean (or use the Julia ability to have
    # multiply-defined names?)

    # these are the barycentre values
    xmean  = mean(x)
    ymean  = mean(y)
    zmean  = mean(z)

    # rbary = barycentre radius
    rbary = ((x.-xmean).^2 + (y.-ymean).^2 + (z.-zmean).^2).^0.5

    # sort these to look at the largest particle?
    return rbary
end

function crosscorrelation(x::Vector{Float32},y::Vector{Float32})
    # y should be optional, to default to autocorrelation?
    #
    # compute the crosscorrelation. will compute for all possible
    # differences in time, but may want to make this an option later.

    nsamples = length(x)

    # initialise a vector
    corr = Vector{Float32}(undef,nsamples)

    for t=1:nsamples
        corr[t] = 0.
        n = 0
        for dt=1:nsamples
            # if greater than the max number of samples, move on
            if t+dt > nsamples
                #print(t," ",dt,"\n")
                continue
            end
        
            n+=1
            corr[t] += x[dt] * y[t+dt]
        end

        # average by the number of valid samples
        print("number of samples:",n,"\n")
        corr[t] /= (n+1.e-60)
    end
    return corr

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



        
function discrete_fourier_k(xmean::Vector{Float32},k::Int)

    # 
    # function to compute the discrete fourier transform of the density centre run with time for a coordinate,
    #   here called 'x'.
    # this is equation (4) of Heggie+ (2020)
    #
    
    # initialise the sum
    sum = 0
    
    # how many numbers in the array?
    tmax = length(xmean)
    
    # loop through the times
    for t=1:tmax
        sum += xmean[t] * exp(-2*1im*pi*t*k/tmax)
    end

    # this will be complex, so need to treat appropriately!
    # e.g. we are interested in the real part of the frequency only
    # normalise and return
    return real(sum/tmax)

       
end



##################################################################################
##################################################################################
# escaper analysis block
##################################################################################
##################################################################################

function plummer_potential(r::Float32,a::Float32=1.0)
    # using G=M=1
    return  -1/sqrt(a*a + r*r)
end

function plummer_potential(r::Vector{Float64},a::Float64=1.0)
    # using G=M=1
    # vectorised.
    return  -1 ./ sqrt.(a^2 .+ (r).^2)
end

function escaper_energy(x::Vector{Float32} ,y::Vector{Float32} ,z::Vector{Float32},
                        vx::Vector{Float32},vy::Vector{Float32},vz::Vector{Float32})
    # assumes positions have already been corrected to the desired frame!
    r  = ( (x).^2 +  (y).^2 +  (z).^2).^0.5
    v2 = ((vx).^2 + (vy).^2 + (vz).^2)
    
    pot = plummer_potential(r)
    
    energy = 0.5*v2 + pot
    
    return energy
end




function spherical_coordinates(x::Vector{Float32} ,y::Vector{Float32} ,z::Vector{Float32},
                               vx::Vector{Float32},vy::Vector{Float32},vz::Vector{Float32})
    # assumes positions have already been corrected to the desired frame!
    r  = ( (x).^2 +  (y).^2 +  (z).^2).^0.5
    r2 = ( (x).^2 +  (y).^2).^0.5
    th = atan.(r2 ./ z)
    ph = atan.(y,x)
        
    vr = (x.*vx + y.*vy + z.*vz) ./ r
    vt = ((x.*vx + y.*vy).*z .- r2.*r2.*vz) ./ (r.*r2)
    vp = (y.*vx + x.*vy) ./ (r2)
    
    return r,th,ph,vr,vt,vp
end

