#=
definitions for the isochrone potential for checking empirical frequencies and other determined properties


=#

function isochrone_psi(r::Float64,M::Float64=1.,bc::Float64=1.)
    # the isochrone potential definition
    rbc = r^2 + bc^2
    return -astronomicalG*M*(bc+sqrt(rbc))^(-1)
end

function isochrone_dpsi_dr(r::Float64,M::Float64=1.,bc::Float64=1.)
    # the analytic isochrone potential derivative
    rbc = r^2 + bc^2
    return astronomicalG*M*r*(sqrt(rbc)*(sqrt(rbc)+bc)^2)^(-1)
end

function isochrone_Omega0(M::Float64=1.,bc::Float64=1.)
    # isochrone frequency scale, from Fouvry 21 (appendix G)
    return sqrt(astronomicalG*M/bc^3)
end

function isochrone_E0(M::Float64=1.,bc::Float64=1.)
    # isochrone energy scale, from Fouvry 21 (appendix G)
    return sqrt(astronomicalG*M/bc)
end
