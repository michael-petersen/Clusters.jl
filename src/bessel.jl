# -- originally written by Jean-Baptiste Fouvry --
# In this section, we re-implement our own Bessel functions
# In essence, we write explicit expressions for the first few harmonics
# This code follows from bessel_j.c from GSL
##################################################
# Spherical Bessel function for \ell=0
# Copied from GSL/bessel_j.c
##################################################
function spherical_bessel_j0(x::Float64)
    if (abs(x) < 0.5) # The argument is close to 0.0
        xsq = x*x
        c1 = -1.0/6.0
        c2 =  1.0/120.0
        c3 = -1.0/5040.0
        c4 =  1.0/362880.0
        c5 = -1.0/39916800.0
        c6 =  1.0/6227020800.0
        return 1.0 + xsq*(c1 + xsq*(c2 + xsq*(c3 + xsq*(c4 + xsq*(c5 + xsq*c6)))))
    else # We are fine with the full expression
        return sin(x)/x
    end
end
##################################################
# Spherical Bessel function for \ell=1
# Copied from GSL/bessel_j.c
##################################################
function spherical_bessel_j1(x::Float64)
    if (abs(x) < 0.25) # The argument is close to 0.0
        xsq = x*x
        c1 = -1.0/10.0
        c2 =  1.0/280.0
        c3 = -1.0/15120.0
        c4 =  1.0/1330560.0
        c5 = -1.0/172972800.0
        ss = 1.0 + xsq*(c1 + xsq*(c2 + xsq*(c3 + xsq*(c4 + xsq*c5))))
        return x/3.0 * ss
    else # We are fine with the full expression
        sin_x, cos_x = sincos(x)
        return (sin_x/x - cos_x)/x
    end
end
##################################################
# Spherical Bessel function for \ell=2
# Copied from GSL/bessel_j.c
##################################################
function spherical_bessel_j2(x::Float64)
    if (abs(x) < 1.3) # The argument is close to 0.0
        xsq  = x*x
        c1 = -1.0/14.0
        c2 =  1.0/504.0
        c3 = -1.0/33264.0
        c4 =  1.0/3459456.0
        c5 = -1.0/518918400.0
        c6 =  1.0/105859353600.0
        c7 = -1.0/28158588057600.0
        c8 =  1.0/9461285587353600.0
        c9 = -1.0/3916972233164390400.0
        ss = 1.0+xsq*(c1+xsq*(c2+xsq*(c3+xsq*(c4+xsq*(c5+xsq*(c6+xsq*(c7+xsq*(c8+xsq*c9))))))))
        return xsq/15.0 * ss
    else # We are fine with the full expression
        sin_x, cos_x = sincos(x)
        ff = (3.0/(x*x) - 1.0)
        return (ff*sin_x - 3.0*cos_x/x)/x
    end
end
##################################################
# Spherical Bessel function for \ell=3
# The values of the coefficients can be found in the notebook Explicit_Expressions_Bessel
# @ATTENTION, THE FUNCTION MISBEHAVES FOR VERY LARGE x ARGUMENTS
##################################################
function spherical_bessel_j3(x::Float64)
    if (abs(x) < 1.3) # The argument is close to 0.0
        xsq = x*x
        c1 = -1.0/18.0
        c2 =  1.0/792.0
        c3 = -1.0/61776.0
        c4 =  1.0/7413120.0
        c5 = -1.0/1260230400.0
        c6 =  1.0/287332531200.0
        c7 = -1.0/84475764172800.0
        c8 =  1.0/31087081215590400.0
        c9 = -1.0/13989186547015680000.0
        ss = 1.0+xsq*(c1+xsq*(c2+xsq*(c3+xsq*(c4+xsq*(c5+xsq*(c6+xsq*(c7+xsq*(c8+xsq*c9))))))))
        return (x*xsq)/105.0 * ss
    else # We are fine with the full expression
        sin_x, cos_x = sincos(x)
        ff = -15.0/(x*x)
        return ((ff+1.0)cos_x - (ff+6.0)*sin_x/x)/x
    end
end
##################################################
# Spherical Bessel function for \ell=4
# The values of the coefficients can be found in the notebook Explicit_Expressions_Bessel
# @ATTENTION, THE FUNCTION MISBEHAVES FOR VERY LARGE x ARGUMENTS
##################################################
function spherical_bessel_j4(x::Float64)
    if (abs(x) < 1.3) # The argument is close to 0.0
        xsq = x*x
        c1 = -1.0/22.0
        c2 =  1.0/1144.0
        c3 = -1.0/102960.0
        c4 =  1.0/14002560.0
        c5 = -1.0/2660486400.0
        c6 =  1.0/670442572800.0
        c7 = -1.0/215882508441600.0
        c8 =  1.0/86353003376640000.0
        c9 = -1.0/41967559641047040000.0
        ss = 1.0+xsq*(c1+xsq*(c2+xsq*(c3+xsq*(c4+xsq*(c5+xsq*(c6+xsq*(c7+xsq*(c8+xsq*c9))))))))
        return (xsq*xsq)/945.0 * ss
    else # We are fine with the full expression
        sin_x, cos_x = sincos(x)
        xsq = x*x
        ff = 105.0/xsq
        return ((-ff+10.0)*cos_x/x + ((ff-45.0)/xsq + 1.0)*sin_x)/x
    end
end
##################################################
# Spherical Bessel function for \ell=5
# The values of the coefficients can be found in the notebook Explicit_Expressions_Bessel
# @ATTENTION, THE FUNCTION MISBEHAVES FOR VERY LARGE x ARGUMENTS
##################################################
function spherical_bessel_j5(x::Float64)
    if (abs(x) < 1.3) # The argument is close to 0.0
        xsq = x*x
        c1 = -1.0/26.0
        c2 =  1.0/1560.0
        c3 = -1.0/159120.0
        c4 =  1.0/24186240.0
        c5 = -1.0/5079110400.0
        c6 =  1.0/1401834470400.0
        c7 = -1.0/490642064640000.0
        c8 =  1.0/211957371924480000.0
        c9 = -1.0/110641748144578560000.0
        ss = 1.0+xsq*(c1+xsq*(c2+xsq*(c3+xsq*(c4+xsq*(c5+xsq*(c6+xsq*(c7+xsq*(c8+xsq*c9))))))))
        return (xsq*xsq*x)/10395.0 * ss
    else # We are fine with the full expression
        sin_x, cos_x = sincos(x)
        xsq = x*x
        ff = 945.0/xsq
        return (((-ff+105.0)/xsq - 1.0)*cos_x + ((ff-420.0)/xsq+15.0)*sin_x/x)/x
    end
end
##################################################
# Wrapped function for the spherical Bessel functions
##################################################
function spherical_bessel_jl(l::Int64,x::Float64)
    if (l==0) # For \ell=0, we have re-implemented the Bessel function
        return spherical_bessel_j0(x)
    elseif (l==1) # For \ell=1, we have re-implemented the Bessel function
        return spherical_bessel_j1(x)
    elseif (l==2) # For \ell=2, we have re-implemented the Bessel function
        return spherical_bessel_j2(x)
    elseif (l==3) # For \ell=3, we have re-implemented the Bessel function
        return spherical_bessel_j3(x)
    elseif (l==4) # For \ell=4, we have re-implemented the Bessel function
        return spherical_bessel_j4(x)
    elseif (l==5) # For \ell=4, we have re-implemented the Bessel function
        return spherical_bessel_j5(x)
    else # For larger values of \ell, we use GSL generic implementation that uses a (slow) recurrence relation
        return sf_bessel_jl(l,x)
    end
end

