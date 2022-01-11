#=
coordintes.jl

a collection of transformations for basic orbit elements

defined terms
-------------------
a  : semi-major axis
e  : eccentricity
rp : pericentre radius
ra : apocentre radius
bc : scaling radius for the cluster

=#


function ae_from_rpra(rp::Float64,ra::Float64)
    #=
    function to translate pericentre and apocentre to semi-major axis and eccentricity
    =#
    return (rp+ra)/2,(ra-rp)/(rp+ra)
end

function rpra_from_ae(a::Float64,e::Float64)
    #=
    function to translate semi-major axis and eccentricity to pericentre and apocentre
    =#
    return a*(1-e),a*(1+e)
end

function spsa_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.)
    # reduced coordinate for pericentre and apocentre
    # includes (optional) scaling factor bc
    xp = rp/bc
    xa = ra/bc
   return sqrt(1+xp^2),sqrt(1+xa^2)
end
