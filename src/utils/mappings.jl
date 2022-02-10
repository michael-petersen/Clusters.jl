


function find_w_min_max(n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.)
    #=
    function to find omega extrema for a given resonance pair

    Fouvry & Prunet B3
    =#

    omega_func(x) = n1*Omega1_circular(dpotential,ddpotential,x) + n2*Omega2_circular(dpotential,x)

    m = extremise_function(omega_func,24,0.,rmax,false) # return a min/max flag?

    # for the problem of a cored cluster with an infinite extent, the w minima and maxima are either
    # in the very centre
    # in the very outskirts
    # along the circular velocity
    w_min,w_max = extrema([omega_func(1.e-10),0.0,omega_func(m)])
    return w_min,w_max
end

function hu(u::Float64,wmin::Float64,wmax::Float64)
    #=
    convenience function, Fouvry & Prunet B8
    =#
    return 0.5*(wmax+wmin + u*(wmax-wmin))
end

function root_of_h_omega(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64,vbound::Float64)
    #=
    function to find the roots of hu - omega_func for a circular orbit; sets constraint on vmin/vmax

    =#
    hval = hu(u,wmin,wmax)
    rootequation(x) = hval - n1*x - n2*x*beta_c[x]
    # two roots to try: bounded by [0,v(u=1)] and [v(u=1),1]
    if rootequation(0+1.e-6)*rootequation(vbound) < 0
        r1 = fzero(rootequation, 0+1.e-6,vbound)
    else
        r1 = 0.0
    end

    if rootequation(1-1.e-6)*rootequation(vbound) < 0
        r2 = fzero(rootequation, vbound,1-1.e-6)
    else
        r2 = 1.0
    end

    # edge curing for vbounds in case of sloppy input
    # from upper left to lower right cure
    if (vbound > 0.999999) & (u==-1.) & (hu(-1.,wmin,wmax)<0)
        r1 = 0.999999
        r2 = 1.0
    end

    # from lower left to upper right cure
    if (vbound > 0.999999) & (u==1.) & (hu(1.,wmin,wmax)>0)
        r1 = 0.999999
        r2 = 1.
    end

    # greater than r1, less than r2
    return r1,r2
end



function constraint_three(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64)
    #=
    constraint to set vmin, vmax limits

    must be greater than this
    =#
    hval = hu(u,wmin,wmax)
    return hval/(n2/2 + n1)
end

function vmin_vmax(u::Float64,wmin::Float64,wmax::Float64,n1::Int64,n2::Int64,vbound::Float64)
    #=
    compute the vmin, vmax boundary for a given u value

    inputs
    -------
    u       : the u coordinate value
    wmin    : the minimum frequency for the resonance pair
    wmax    : the maximum frequency for the resonance pair
    n1      : the first resonance integer
    n2      : the second resonance integer
    vbound  : the value of v at u +\- 1; sets the root-finding limits

    =#

    if (n2==0)

        hval = hu(u,wmin,wmax)

        vmin = 0.5

        # put in guards for the very edges. SLOPPY
        vmax = beta_c(minimum([0.99999,maximum([hval/n1,0.00001])]))

    else

        r1,r2 = root_of_h_omega(u,wmin,wmax,n1,n2,vbound)

        r3 = constraint_three(u,wmin,wmax,n1,n2)

        if (abs(n2)<=abs(n1)) | (abs(n2)>=abs(2n1))
            vmin = maximum([0.0,r1])
            vmax = minimum([r2,1.0,r3])
        else
            vmin = maximum([0.0,r1,r3])
            vmax = minimum([r2,1.0])
        end

        if isnan(vmax)
            vmax = 1.
        end

        if isnan(vmin)
            vmin = 0.
        end
    end

    return vmin,vmax
end

function alphabeta_from_uv(u::Float64,v::Float64,
                           n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.)
    #=
    the inverse mapping from (u,v) -> (alpha,beta)

    weirdly imperfect: why?

    =#

    wmin,wmax = find_w_min_max(n1,n2,dpotential,ddpotential,rmax)

    if n2 == 0
        beta  = v
        alpha = (1/(2n1))*((wmax-wmin)*u + wmin + wmax)
    else
        alpha = v
        beta  = (1/(n2*v))*(0.5*((wmax-wmin)*u + wmin + wmax) - n1*v)
    end

    return alpha,beta
end

function uv_from_alphabeta(alpha::Float64,beta::Float64,
                           n1::Int64,n2::Int64,dpotential::Function,ddpotential::Function,rmax::Float64=1000.)
    #=
    the mapping from (alpha,beta) -> (u,v)

    =#
    wmin,wmax = find_w_min_max(n1,n2,dpotential,ddpotential,rmax)

    wval = n1*alpha + n2*beta*alpha

    u = (2*wval - wmax - wmin)/(wmax-wmin)

    if (n2==0)
        v = beta
    else
        v = alpha
    end

    return u,v

end
