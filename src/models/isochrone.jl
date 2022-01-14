#=
definitions for the isochrone potential for checking empirical frequencies and other determined properties

todo
-the Jacobian can definitely be improved

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


function isochrone_omega_ae(rp::Float64,ra::Float64,bc::Float64=1.)
    # isochrone \omega (Fouvry 21 eq. G5)
    sp,sa = spsa_from_rpra(rp,ra)
    return (2/(sp+sa))^(3/2)
end

function isochrone_eta_ae(rp::Float64,ra::Float64,bc::Float64=1.)
    # isochrone \eta (Fouvry 21 eq. G7)
    xp = rp/bc
    xa = ra/bc
    sp,sa = spsa_from_rpra(rp,ra)
    return (1/2)*(1+(xp*xa)/((1+sp)*(1+sa)))
end

function isochrone_Omega0(M::Float64=1.,bc::Float64=1.)
    # isochrone frequency scale, from Fouvry 21 (appendix G)
    return sqrt(astronomicalG*M/bc^3)
end

function isochrone_E0(M::Float64=1.,bc::Float64=1.)
    # isochrone energy scale, from Fouvry 21 (appendix G)
    return sqrt(astronomicalG*M/bc)
end

function isochrone_L0(M::Float64=1.,bc::Float64=1.)
    # isochrone action scale, from Fouvry 21 (appendix G)
    return sqrt(astronomicalG*M*bc)
end

function isochrone_Omega_1_2(rp::Float64,ra::Float64,bc::Float64=1.)
    # wrapper for returning both \Omega_1 and \Omega_2

    Omega0   = isochrone_Omega0()
    omega_ae = isochrone_omega_ae(rp,ra)
    eta_ae   = isochrone_eta_ae(rp,ra)
    return omega_ae*Omega0,omega_ae*eta_ae*Omega0

end

function isochrone_Omega_1_2_ae(a::Float64,ecc::Float64,bc::Float64=1.)
    # wrapper for returning both \Omega_1 and \Omega_2

    rp,ra    = rpra_from_ae(a,ecc)
    Omega0   = isochrone_Omega0()
    omega_ae = isochrone_omega_ae(rp,ra)
    eta_ae   = isochrone_eta_ae(rp,ra)
    return omega_ae*Omega0,omega_ae*eta_ae*Omega0

end



function isochrone_E_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.)
    # isochrone analytic energy, Fouvry 21 G9
    E0 = isochrone_E0()
    sp,sa = spsa_from_rpra(rp,ra)
    return E0/(sp+sa)
end

function isochrone_L_from_rpra(rp::Float64,ra::Float64,bc::Float64=1.)
    # isochrone analytic energy, Fouvry 21 G9
    xp = rp/bc
    xa = ra/bc
    L0 = isochrone_L0()
    sp,sa = spsa_from_rpra(rp,ra)
    return sqrt(2)*L0*xp*xa/sqrt((1+sp)*(1+sa)*(sp+sa))
end

function isochrone_dthetadu_from_rpra(r::Float64,u::Float64,rp::Float64,ra::Float64,bc::Float64=1.)
    # the isochrone analytic Jacobian, Fouvry 21 G10
    xp = rp/bc
    xa = ra/bc
    xr = r/bc
    sr = sqrt(1+xr^(2))
    Omega0 = isochrone_Omega0()
    Omega1,Omega2 = isochrone_Omega_1_2(rp,ra)
    sp,sa = spsa_from_rpra(rp,ra)
    return (3/sqrt(2))*(Omega1/Omega0)*(xr/sqrt(4-u^(2)))*(sqrt((sr+sp)*(sr+sa)*(sp+sa))/sqrt((xr+xp)*(xr+xa)))
end
