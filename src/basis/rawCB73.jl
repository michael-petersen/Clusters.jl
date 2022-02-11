##################################################
# Construction of the radial basis elements
# from Clutton-Brock (1973)
##################################################
# Precomputing the prefactors of the Gegenbauer functions
#####
const Rbasis = 1.
const nradial = 12
const lmax    = 6
const G = 1
const M = 1
const tabPrefCB73_Ulnp = zeros(Float64,lmax+1,nradial) # Table of the prefactors of the POTENTIAL basis functions, i.e. the Ulnp !! ATTENTION, size is lmax+1 as l=0 exists
const tabPrefCB73_Dlnp = zeros(Float64,lmax+1,nradial) # Table of the prefactors of the DENSITY   basis functions, i.e. the Dlnp !! ATTENTION, size is lmax+1 as l=0 exists
#####
# Reading the table of the pre-computed prefactors of the Gegenbauer functions
#####
const lmaxCB73 = 50 # Maximum ell for the pre-computed prefactors
const nmaxCB73 = 200 # Maximum np  for the pre-computed prefactors
filename = "data_Basis/data_CB73_lmax_"*string(lmaxCB73)*"_nmax_"*string(nmaxCB73)*".h5" # Name of the file where the prefactors were dumped
const namefileCB73 = joinpath(@__DIR__, filename)
const tabalphaCB73 = h5read(namefileCB73,"tab_alphalnp") # Reading the prefactors alpha_lnp
const tabbetaCB73  = h5read(namefileCB73,"tab_betalnp")  # Reading the prefactors beta_lnp
##################################################
# Function that fills in the prefactors
# of the Gegenbauer functions
##################################################
function tabPrefCB73!()
    for l=0:lmax # Loop over the harmonic indices. ATTENTION, harmonic index starts at l=0
        for np=1:nradial # Loop over the radial basis numbers
            alpha = tabalphaCB73[l+1,np] # Reading the value of alpha. ATTENTION, l starts at l=0
            beta  = tabbetaCB73[l+1,np]  # Reading the value of beta.  ATTENTION, l starts at l=0
            #####
            A = sqrt(G/Rbasis)*alpha # Value of the prefactor A_n^\ell
            B = 1.0/(sqrt(G)*Rbasis^(5/2))*beta # Value of the prefactor B_n^\ell
            #####
            tabPrefCB73_Ulnp[l+1,np] = A # Filling in the array. ATTENTION, l starts at l=0
            tabPrefCB73_Dlnp[l+1,np] = B # Filling in the array. ATTENTION, l starts at l=0
        end
    end
end
##################################################
tabPrefCB73!() # Filling in the arrays of prefactors
##################################################
# Function that returns the parameter -1 <= rho <= 1
# for a given dimensionless radius x=r/Rbasis
##################################################
function rhoCB73(x::Float64)
    return (x^(2) - 1.0)/(x^(2) + 1.0) # Value of rho
end
##################################################
# Definition of the Gegenbauer polynomials
# These coefficients are computed through
# an upward recurrence
# @ IMPROVE -- Of course, it is much better
# to compute all the basis elements (n)_{1<=n<=nradial} at once
##################################################
function ClnCB73(alpha::Float64,n::Int64,rho::Float64)
    #####
    v0 = 1.0 # Initial value for n=0
    if (n == 0)
        return v0 # No need for a recurrence for n=0
    end
    #####
    v1 = 2.0*alpha*rho # Initial value for n=1
    if (n == 1)
        return v1 # No need for a recurrence for n=1
    end
    #####
    ic = 2 # Iteration counter that gives the index of the value that is about to be computed
    v = 0.0 # Initialisation of the temporary variable
    #####
    while (ic <= n) # Applying the recurrence as many times as needed
        v = (2.0*(ic+alpha-1.0)*rho*v1 - (ic+2.0*alpha-2.0)*v0)/(ic) # Applying the recurrence
        v0, v1 = v1, v # Updating the temporary variables
        ic += 1 # Updating the counter of iteration
    end
    #####
    return v # Output of the value
end
##################################################
# Definition of the radial basis elements
# from Clutton-Brock (1973)
# Arguments are:
# + l:harmonic index
# + np: radial index. ATTENTION, this index starts at n=1
# + r: radius
##################################################
function UlnpCB73(l::Int64,np::Int64,r::Float64)
    pref = tabPrefCB73_Ulnp[l+1,np] # Value of the prefactor. ATTENTION, l starts at l=0
    x = r/Rbasis # Dimensionless radius
    rho = rhoCB73(x) # Value of the rescaled parameter rho
    valR = ((x/(1.0+x^(2)))^(l))/(sqrt(1.0+x^(2))) # Value of the multipole factor
    valC = ClnCB73(l+1.0,np-1,rho) # Value of the Gegenbauer polynomials
    res = pref*valR*valC # Value of the radial function
    return res # Output
end
#####
function DlnpCB73(l::Int64,np::Int64,r::Float64)
    pref = tabPrefCB73_Dlnp[l+1,np] # Value of the prefactor. ATTENTION, l starts at l = 0
    x = r/Rbasis # Dimensionless radius
    rho = rhoCB73(x) # Value of the rescaled parameter rho
    valR = ((x/(1.0+x^(2)))^(l))/((1.0+x^(2))^(5/2)) # Value of the multipole factor
    valC = ClnCB73(l+1.0,np-1,rho) # Value of the Gegenbauer polynomials
    res = pref*valR*valC
    return res # Output
end
##################################################
# Function that computes the values of Ulnp(r)
# for a given l and r,
# and for 1 <= np <= nradial
# This is particularly useful for the Gegenbauer polynomials
# as those are computed through a recurrence
##################################################
function tabUlnpCB73!(l::Int64,r::Float64,
                       tabUlnp::Array{Float64,1})
    #####
    x = r/Rbasis # Dimensionless radius
    rho = rhoCB73(x) # Value of the rescaled parameter rho
    valR = ((x/(1.0+x^(2)))^(l))/(sqrt(1.0+x^(2))) # Value of the multipole factor
    #####
    alpha = l+1.0 # Value of alpha, the index of the Gegenbauer polynomials
    #####
    v0 = 1.0 # Initial value of the Gegenbauer polynomials for n=0
    v1 = 2.0*alpha*rho # Initial value of the Gegenbauer polynomials for n=1
    #####
    tabUlnp[1] = tabPrefCB73_Ulnp[l+1,1]*valR*v0 # Filling in the value for np=1. ATTENTION, l starts at l=0
    tabUlnp[2] = tabPrefCB73_Ulnp[l+1,2]*valR*v1 # Filling in the value for np=2. ATTENTION, l starts at l=0
    #####
    for np=3:nradial # Loop over the radial indices
        v = (2.0*(np+alpha-2.0)*rho*v1 - (np+2.0*alpha-3.0)*v0)/(np-1) # Applying the recurrence
        v0, v1 = v1, v # Updating the temporary variables
        tabUlnp[np] = tabPrefCB73_Ulnp[l+1,np]*valR*v # Filling in the value for np. ATTENTION, l starts at l=0
    end
end
