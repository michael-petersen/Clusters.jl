#=
Tools for the Henon anomaly analysis

=#

function drdu(u::Float64,a::Float64,ecc::Float64,TOL::Float64=0.01)
    dr = (3/2)*(a*ecc)*(1-(u^2))
    ddr = -3a*ecc*u

    if (1 - u^2)<TOL
        # expand from value
        # is there a point to this? Should always be well-defined.
        if u > 0 # u near 1
            h = u-1
            dr = h*ddr
        else # u near -1
            h = u+1
            dr = h*ddr
        end
    end

    return dr

end
