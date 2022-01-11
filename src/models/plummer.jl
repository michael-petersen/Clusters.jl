#=
definitions for the plummer potential

=#


function plummer_psi(r::Float64,M::Float64=1.,bc::Float64=1.)
    # the plummer potential definition
    rbc = r^2 + bc^2
    return -astronomicalG*M*(sqrt(rbc))^(-1)
end

function plummer_dpsi_dr(r::Float64,M::Float64=1.,bc::Float64=1.)
    # the analytic plummer potential derivative
    rbc = r^2 + bc^2
    return astronomicalG*M*r*((rbc)^(-3/2))
end
