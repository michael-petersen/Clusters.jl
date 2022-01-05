#=
cluttonbrock.jl

All the machinery needed to implement the spherical Clutton-Brock (1973) basis.

Currently, only

-potential_function
-density_function
-accumulate_ln

are exported.

=#

using Statistics                  # access to mean
#using SphericalHarmonics # access to Spherical harmonics (but watch out for conventions!)

function mapping_r_xi(r::Float64,bc::Float64)
    # function to compute the reduced mapping
    red_r = r/bc
    return (red_r^2 - 1)/(red_r^2+1)
end

function mapping_xi_r(r::Float64,bc::Float64)
    # function to undo the reduced mapping
    return bc*sqrt((1+xi)/(1-xi))
end


function potential_function(r::Float64,l::Int64,n::Int64,bc::Float64=1.)
    # function to return an arbitrary Clutton-Brock potential function
    rs = r/bc
    xi = mapping_r_xi(r,bc)
    C = upward_gegenbauer(xi,l+1,n+1)
    return -(rs^l)*C[n+1]/((1+rs^2)^(l+0.5))
end

function density_function(r::Float64,l::Int64,n::Int64,bc::Float64=1.)
    # function to return an arbitrary Clutton-Brock density function
    rs = r/bc
    xi = mapping_r_xi(r,bc)
    C = upward_gegenbauer(xi,l+1,n+1)
    return (rs^l)*C[n+1]*Kln(l,n)/((1+rs^2)^(l+2.5))/(4pi)
end


function accumulate_ln(xvals,
                       yvals,
                       zvals,
                       l::Int64,n::Int64,m::Int64,
                       bc::Float64)
    # function to accumulate the monopole of the Clutton-Brock model
    r = ((xvals).^2 + (yvals).^2 + (zvals).^2).^0.5
    wln = Wln(l,n,m)

    C = 0.

    for a=1:length(r)
        C += potential_function(r[a],l,n,bc)
    end

    return wln*C/length(r)
end


function Kln(l::Int64,n::Int64)
    # the normalising constant for the Clutton-Brock basis
    return 4n*(n+2l+2) + (2l+1)*(2l+3)
end

function Wln(l::Int64,n::Int64,m::Int64=0)
    # the normalisation for the potential-weighted accumulation
    # in this nomenclature, n goes from 0->nmax-1
    m==0        ? delta=1 : delta=2
    knl         = Kln(l,n)
    numerator   = (l+n+1)*factorial(n)*factorial(l)^2
    denominator = knl*factorial(2l+n+1)
    prefac      = 2^(4l+6)*delta
    return prefac*numerator/denominator
end

function next_gegenbauer(rho::Float64,alpha::Int64,n::Int64,C1::Float64,C2::Float64)
    # recursive formula for the next Gegenbauer function
    return (1/n)*(2(n-1+alpha)*rho*C1 - (n-1+2alpha-1)*C2)
end

function upward_gegenbauer(rho::Float64,alpha::Int64,nmax::Int64)
    # the function to accumulate Gegenbauer prefactors
    C = zeros(Float64,nmax)
    C[1] = 1
    if nmax>1
        C[2] = 2alpha*rho
        if nmax>2
            for n=3:nmax
                C[n] = next_gegenbauer(rho,alpha,n,C[n-1],C[n-2])
            end
        end
    end
    return C
end
