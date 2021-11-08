# Find empirically-estimated frequencies given a potential and gradient.
##################################################
using Optim # for zero-finding

##################################################
# some basic potential options
##################################################
function isochrone_psi(r::Float64)
    M = 1.
    G = 1.
    bc= 1.
    rbc = r^2 + bc^2
    return -G*M*(bc+sqrt(rbc))^(-1)
end

function isochrone_dpsi_dr(r::Float64)
    M = 1.
    G = 1.
    bc= 1.
    rbc = r^2 + bc^2
    return G*M*r*(sqrt(rbc)*(sqrt(rbc)+bc)^2)^(-1)
end

function plummer_psi(r::Float64)
    M = 1.
    G = 1.
    bc= 1.
    rbc = r^2 + bc^2
    return -G*M*(sqrt(rbc))^(-1)
end

function plummer_dpsi_dr(r::Float64)
    M = 1.
    G = 1.
    bc= 1.
    rbc = r^2 + bc^2
    return G*M*r*((rbc)^(-3/2))
end



##################################################
# the potential wrappers: select your potential here
##################################################
function dpotential_numerical(r::Float64,eps::Float64=0.0000001)
    # rough numerical derivative, but cheap
    # making eps any smaller incurs numerical noise
    return (potential(r+eps)-potential(r))/eps
end

function potential(r::Float64)
    return isochrone_psi(r)
end    
    
function dpotential(r::Float64)
    return isochrone_dpsi_dr(r)
end


##################################################
# some basic orbit functions
##################################################

function find_j(r::Float64,kappa::Float64)
    # compute the angular momentum for a given orbit at radius r, with given kappa
    dudr = dpotential(r)
    jmax = sqrt(r*r*r*dudr);
    J = jmax*kappa;
    return J
end


function Ecirc(r::Float64,E::Float64)
    # compute the energy of a circular orbit at some radius
    # must define the potential and potential derivative a priori
    ur = potential(r)
    dudr = dpotential(r)
    return  abs(E - 0.5*r*dudr - ur)
end


function denom(r::Float64,E::Float64,J::Float64)
    # the main function to root-find, this is the effective potential
    ur = potential(r)
    return abs(2.0*(E-ur)*r*r - J*J)
end



##################################################
# the wrapper to convert E,K to r_peri,r_apo
##################################################

function make_orbit(E::Float64,K::Float64,rmax::Float64=1000.)
    # initialise an orbit in E,K space
    # will not check that E is valid for the model!
    # rmax must be defined in order to recover appropriate roots
    r_circ = optimize(r -> Ecirc(r,E), 0.    ,rmax  , Brent()).minimizer
    J      = find_j(r_circ,K)
    r_apo  = optimize(r -> denom(r,E,J), r_circ,rmax  , Brent()).minimizer
    r_peri = optimize(r -> denom(r,E,J), 0.    ,r_circ, Brent()).minimizer
    return r_peri,r_apo,r_circ,J
end



##################################################
# the workhorse: take orbit-defining parameters and return the frequencies
##################################################

function kepler_anomaly_frequencies(r_apo::Float64,r_peri::Float64,ee::Float64,jj::Float64)
    # set the integration width
    FRECS = 16
    # ,FRECS::Int32=16
    dt = pi/FRECS;
    
    # define some auxilliaries
    ap = 0.5*(r_apo + r_peri);
    am = 0.5*(r_apo - r_peri);
    sp = ap/(r_apo*r_peri);
    sm = am/(r_apo*r_peri);
    
    # set the accumulators to zero
    #accum0 = 0.0;
    accum1 = 0.0;
    accum2 = 0.0;
    
    # proceed as centred rectangle integration
    t = 0.5*(dt-pi)
  
    for i=1:FRECS

        r = ap + am*sin(t)
      
        ur = potential(r)
        cost = cos(t)
        
        tmp = sqrt(2.0*(ee-ur) - (jj*jj)/(r*r));
        
        #accum0 += cost * tmp;
      
        accum1 += cost / tmp;
      
        s = sp + sm*sin(t);
      
        ur = potential(1.0/s)
      
        accum2 += cost/sqrt(2.0*(ee-ur) - (jj*jj*s*s));
      
        # advance the counter
        t += dt
    end

  
    freq1 = pi/(am*accum1*dt);
    freq2 = freq1 * jj * sm * accum2 * dt/pi;
    
    # note that we could also compute the actions if we wanted:
    #action1 = am*accum0*dt/pi;
    #action2 = jj;
    
    return freq1,freq2
    
end


##################################################
# the wrappers for interfacing and getting frequencies
##################################################

function compute_frequencies(r_apo::Float64,r_peri::Float64,ee::Float64,jj::Float64)
    
    freq1,freq2 = kepler_anomaly_frequencies(r_apo,r_peri,ee,jj)
    
    return freq1,freq2
end
    
function compute_frequencies(r_apo::Float64,r_peri::Float64)
    
    # ee, jj are dangerous for circular orbits... may need the limited development here
    ee = (r_apo*r_apo*potential(r_apo) - r_peri*r_peri*potential(r_peri))/(r_apo^2 - r_peri^2)
    jj = sqrt(2*(potential(r_apo) - potential(r_peri))/(r_peri^(-2) - r_apo^(-2)))
    #print("ee/jj",ee," ",jj," ","\n")
    
    freq1,freq2 = kepler_anomaly_frequencies(r_apo,r_peri,ee,jj)
    
    return freq1,freq2
end

function compute_frequencies_EK(E::Float64,K::Float64)
    # put the computations together so we only need E,K
    
    r_peri,r_apo,r_circ,J = make_orbit(E,K)
    
    freq1,freq2 = kepler_anomaly_frequencies(r_apo,r_peri,E,J)
    
    return freq1,freq2
end
    

# https://stackoverflow.com/questions/46842510/how-to-pass-a-function-as-an-argument-for-another-function-in-julia
# https://docs.julialang.org/en/v1/manual/performance-tips/

