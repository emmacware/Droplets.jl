#---------------------------------------------------------------------------
#Condensation functions in this file (more description throughout the file)
#---------------------------------------------------------------------------
# using DifferentialEquations
export esat,sat,drdtcondensation1,drdtcondensation2,drdtcondensation3#,condense_and_calc_Sv!
export FK,FD,drkohler,θcondenseupdate!,qvcondenseupdate!,dXkohler_function_of_radius
export dXkohler_function_of_radius_activated,drkohler_activated
export set_X_crit!,find_equilibrium_radius
export dM_dt
export drkappakohler
export condensation_rhs!, condensation_rhs_single_cell, condensation_rhs_single_droplet
export calc_cond_rad_term, dXkappakohler_bisection
export limitsat, condensation_time_step_spatial!, apply_thermo_feedback!
######################################################################
# 

"""
    FK(T, constants)
    FD(T, constants)
    Functions in the Denominator of the Köhler Equation

Uses the constants structs defined in Droplets
# Arguments
- `T`: Temperature in K
"""
function FK(T,constants) 
    Fk = (constants.L /(constants.Rv *T) -1)*(constants.L*constants.ρl)/(constants.k *T) 
    return Fk
end

function FD(T,constants)
    Fd = constants.ρl*constants.Rv *(T)/(constants.Dv * esat(T))
    return Fd
end

@inline akk(T) = 3.3 * 10^(-7) / T #(m K/T)
@inline bkk(m) = 4.3 * 2 / m / 1e6 #m^3 for NaCL

######################################################################
# Saturation Functions

#saturation vapor pressure
@inline esat(T) = 610.94 * exp(17.625 * (T - 273.15) / ((T - 273.15) + 243.04)) ## August–Roche–Magnus approximation, hPa to Pa, T in K, converted to C

#environmntal vapor pressure based on the mixing ratio of water vapor to air
@inline sat(qv::FT,P::FT) where FT<:AbstractFloat = FT(qv*P/(FT(0.378)*qv + FT(0.622))) ##qvarray is the specific humidity: mixing ratio of water vapor to moist air, Pa


######################################################################
# The following ode functions are for the condensation radial growth of a droplet
#examples for using DifferentialEquations
#drdtcondensation1 p = (a,b,S,M,denom)
#drdtcondensation2 p = (M,m,T,qv,P)
#drdtcondensation3 p = (M,m,T,Senv)
#where S is environmental saturation, M is the mass of the solute,
#m is the molecular weight of the solute, T is the temperature, 
#qv is the mixing ration, and P is the pressure

function drdtcondensation1(u, p, t)
    a, b, S, M,denom = p
    return (S-1-(a ./u) .+b.*M./(u.^3))./(denom.*u)
end

function drdtcondensation2(u,p,t)
    M,m,T,qv,P = p
    du = drkohler(u,M,m,T,qv,P,t)
    return du
end

function drdtcondensation3(u,p,t)
    M,m,T,Senv = p
    du = drkohler(u,M,m,T,Senv,t)
    return du
end

######################################################################
#This function calculate the RHS of the Köhler equation
# drdt = ___
# Can take either the (mixing ratio and pressure) or (environmental saturation) as an argument
"""
Methods:
    -drkohler(R, M, m, T, qv, P, timestep)
    -drkohler(R, M, m, T, Senv, timestep)

RHS of the Köhler equation.

# Arguments
- `R`: Droplet radius (m)
- `M`: Molecular weight of the solute (kg/mol)
- `m`: Molar mass of the solute (g/mol)
- `T`: Temperature (K)
- `qv`: Water vapor mixing ratio (kg/kg)
- `P`: Atmospheric pressure (Pa)
- `Senv`: Environmental saturation (dimensionless)
- `timestep`: Time step (s)

# Returns
- `dr`: Change in droplet radius over the timestep (m)
"""
function drkohler(R, M, m, T, qv, P,constants, timestep)
    S = sat.(qv, P) ./ esat(T)
    denom = (FK(T, constants) + FD(T, constants))
    dr = (S - 1 .- (akk(T) ./ R) .+ bkk.(m) .* M ./(R .^ 3)) ./(denom .* R)
    return dr #R + dr * timestep > 0 ? dr : -R / timestep
end

function drkohler(R, M, m, T, Senv,constants, timestep)
    denom = (FK(T, constants) + FD(T, constants))
    dr = (Senv - 1 .- (akk(T) ./ R) .+ bkk.(m) .* M ./(R .^ 3)) ./(denom .* R)
    return dr #R + dr * timestep > 0 ? dr : -R / timestep
end

function drkappakohler(R,dry_r3,kappa,T,Senv,constants,timestep;rad_term=0.0)
    b = kappa * dry_r3
    # M = 4/3 * π * dry_r3 * ρ_solute
    fk = FK(T, constants)
    denom = (fk + FD(T, constants))
    # dr = (Senv - 1 .- (akk(T) ./ R) .+ b .* M ./(R .^ 3)) ./(denom .* R)
    dr = (Senv - 1.0 - (akk(T) / R) + (kappa * dry_r3) / R^3) / (denom * R)
    
    #if rad_term * fk / (constants.L*constants.ρl * denom) > 1e-6
    #    println(dr,' ',rad_term * fk / (constants.L*constants.ρl * denom))
    #end

    # dr += rad_term * fk / (constants.L*constants.ρl * denom) 
        

    return dr#R + dr * timestep > 0 ? dr : -R / timestep
end

function dXkappakohler_bisection(REM::DynOFF, droplets::droplet_attributes{FT}, i::Int, kappa::FT, T::FT, Senv::FT, c::FT, radcoeff::FT, constants::Constants{FT},
    raddata::Rad,
    absliq_r_interp::Abs,
    timestep::FT, iters::Int, depth::Int=0, max_depth::Int=6) where {FT, Abs,Rad}
    R      = volume_to_radius(droplets.X[i])
    dry_r3 = droplets.dry_r3[i]
    dry_r  = cbrt(dry_r3)
    A      = akk(T)
    b      = kappa * dry_r3
    S1     = Senv - one(FT)
    R2     = R * R

    # Scrit = critical (haze-max) supersaturation for this droplet's dry size/kappa.
    # Past it there is no equilibrium root at all (the droplet is unconditionally
    # activated), so don't waste recursion trying to bracket one — drop the
    # Kelvin/solute terms (negligible once past the hump) and integrate the
    # remaining smooth diffusional growth law analytically.
    r_crit = sqrt(3 * b / A)
    Scrit  = one(FT) + sqrt(4 * A^3 / (27 * b))
    if Senv > Scrit
        r_new = sqrt(max(R2 + FT(2) * c * S1 * timestep, dry_r * dry_r))
        droplets.X[i] = FT(4π/3) * r_new^3
        return
    end

    F0     = c * (S1 -A/R + b/R^3)

    # bisect in r-space: g(r) = r² - R²_old - Δt·c·(S1 + (b - A·r²)/r³)
    # one sqrt for the non-trivial bracket bound; none inside the loop
    if F0 >= 0
        lo_r = max(dry_r, R)
        hi_r = sqrt(R2 + FT(2) * F0 * timestep)
        # only cap at r_crit while still searching below the hump (unactivated
        # haze droplet): once already past it (R > r_crit), this is just an
        # ordinary large droplet doing normal diffusional growth and r_crit no
        # longer bounds anything — capping there would invert the bracket
        lo_r < r_crit && (hi_r = min(hi_r, r_crit))
    else
        g_dry = dry_r * dry_r - R2 - timestep * c * (S1 - A / dry_r + kappa)
        g_R   = -timestep * c * F0
        if g_dry * g_R >= 0
            droplets.X[i] = FT(4π/3) * dry_r^3
            return
        end
        lo_r = dry_r
        hi_r = R
    end

    lo2  = lo_r * lo_r;#  lo3 = lo2 * lo_r
    g_lo = lo2 - R2 - timestep * c * (S1 - A / lo_r + b / (lo2 * lo_r))
    hi2  = hi_r * hi_r;  hi3 = hi2 * hi_r
    g_hi = hi2 - R2 - timestep * c * (S1 - A / hi_r + b / (hi2 * hi_r))
    r_new = (lo_r + hi_r) * FT(0.5)


    if g_lo * g_hi >= 0

        if depth < max_depth
            half = timestep * FT(0.5)
            dXkappakohler_bisection(REM, droplets, i, kappa, T, Senv, c, radcoeff, constants, raddata, absliq_r_interp, half, iters, depth + 1, max_depth)
            dXkappakohler_bisection(REM, droplets, i, kappa, T, Senv, c, radcoeff, constants, raddata, absliq_r_interp, half, iters, depth + 1, max_depth)
            return
        end
        r_new = max(r_new, dry_r)
        droplets.X[i] = FT(4π/3) * r_new^3
        return
    end

    for _ in 1:iters
        r2 = (lo2 + hi2) * FT(0.5)
        # r    = sqrt(r2)
        r3    = r2 * sqrt(r2)
        g_mid = r2 - R2 - timestep * c * (S1 + (b - A * r2) / r3)
        if g_lo * g_mid <= 0
            hi2 = r2;  g_hi = g_mid
        else
            lo2 = r2;  g_lo = g_mid
        end
        (hi2 - lo2) < 1e-4 * r2 && break
    end
    r2 = (lo2 + hi2) * FT(0.5)

    r_new = max(sqrt(r2), dry_r)
    droplets.X[i] = FT(4π/3) * r_new^3
    return
end


function calc_cond_rad_term(R::FT,z::Int,constants::Constants{FT},raddata::Rad,absliq_r_interp::Abs) where {FT,Rad,Abs}
    if R < 2.5e-6 || R > 4e-5
        return 0.0
    end
    abs_fn = 0.0
    n_bnd = raddata.CArad.nband_lw
    # tot_wn = 0.0
    tot_flux = 0.0

    # if R > 21.5e-6
    #     R = 21.5e-6
    # end

    for ibnd in 1:n_bnd
        abs_drop = absliq_r_interp[ibnd](R*1.0e6)
        # Qa = abs_drop * 4/3 * R * constants.ρl * kg_to_g
        # Δwn = raddata.CArad.lookup_lw_cld.bnd_lims_wn[2, ibnd] - raddata.CArad.lookup_lw_cld.bnd_lims_wn[1, ibnd]

        abs_fn += abs_drop * raddata.flux_net_droplet[z,ibnd]
        # tot_flux += abs(raddata.flux_net_droplet[z,ibnd])
    end
    # abs_fn /= tot_flux

    rad_term = abs_fn * 4/3 * R * constants.ρl * kg_to_g #implicit ρlconversion for the multiplication and division of ρl

    return rad_term * 2 * R #/ constants.L
end

function dXkappakohler_bisection(REM::DynON, droplets::droplet_attributes{FT}, i::Int, kappa::FT, T::FT, Senv::FT, c::FT, radcoeff::FT, constants::Constants{FT},
    raddata::Rad,
    absliq_r_interp::Abs,
    timestep::FT, iters::Int, depth::Int=0, max_depth::Int=6) where {FT, Abs,Rad}

    depth == 0 && (raddata.cond_rad_term[i] = zero(FT))

    X_old  = droplets.X[i]
    R      = volume_to_radius(X_old)
    dry_r3 = droplets.dry_r3[i]
    dry_r  = cbrt(dry_r3)
    z      = droplets.cell_id[i]
    A      = akk(T)
    b      = kappa * dry_r3
    S1     = Senv - one(FT)
    R2     = R * R

    r_crit = sqrt(3 * b / A)
    Scrit  = one(FT) + sqrt(4 * A^3 / (27 * b))
    if Senv > Scrit
        rad0 = radcoeff*calc_cond_rad_term(R, z,constants, raddata, absliq_r_interp)
        Fcond0 = c * S1
        F0 = Fcond0 + rad0
        r_new = sqrt(max(R2 + FT(2) * F0 * timestep, dry_r * dry_r))
        X_new = FT(4π/3) * r_new^3

        dX_total = X_new - X_old
        dX_cond  = iszero(F0) ? dX_total : dX_total * Fcond0 / F0

        droplets.X[i] = X_new
        raddata.cond_rad_term[i] += dX_total - dX_cond
        return
    end

    rad0 = radcoeff*calc_cond_rad_term(R, z,constants, raddata, absliq_r_interp)
    Fcond0   = c * (S1 + (b - A * R2) / (R2 * R))
    F0 = Fcond0 + rad0

    if F0 >= 0
        lo_r = max(dry_r, R)
        hi_r = sqrt(R2 + FT(2) * F0 * timestep)
        lo_r < r_crit && (hi_r = min(hi_r, r_crit))
    else
        g_dry = dry_r * dry_r - R2 - timestep * (c * (S1 - A / dry_r + kappa) +
                radcoeff*calc_cond_rad_term(dry_r, z, constants, raddata, absliq_r_interp))
        g_R   = -timestep * F0
        if g_dry * g_R >= 0
            X_new = FT(4π/3) * dry_r^3
            dX_total = X_new - X_old
            dX_cond  = iszero(F0) ? dX_total : dX_total * Fcond0 / F0
            droplets.X[i] = X_new
            raddata.cond_rad_term[i] += dX_total - dX_cond
            return
        end
        lo_r = dry_r
        hi_r = R
    end

    lo2  = lo_r * lo_r;  lo3 = lo2 * lo_r
    g_lo = lo2 - R2 - timestep * (c * (S1 + (b - A * lo2) / lo3) + radcoeff*calc_cond_rad_term(lo_r, z,constants, raddata, absliq_r_interp))
    hi2  = hi_r * hi_r;  hi3 = hi2 * hi_r
    g_hi = hi2 - R2 - timestep * (c * (S1 + (b - A * hi2) / hi3) + radcoeff*calc_cond_rad_term(hi_r, z, constants, raddata, absliq_r_interp))
    r2_new = (lo2 + hi2) * FT(0.5)

    if g_lo * g_hi >= 0
        if depth < max_depth
            half = timestep * FT(0.5)
            dXkappakohler_bisection(REM, droplets, i, kappa, T, Senv, c, radcoeff, constants, raddata, absliq_r_interp, half, iters, depth + 1, max_depth)
            dXkappakohler_bisection(REM, droplets, i, kappa, T, Senv, c, radcoeff, constants, raddata, absliq_r_interp, half, iters, depth + 1, max_depth)
            return
        end
        r_new  = max(sqrt(r2_new), dry_r)
        X_new  = FT(4π/3) * r_new^3

        dX_total = X_new - X_old
        dX_cond  = iszero(F0) ? dX_total : dX_total * Fcond0 / F0

        droplets.X[i] = X_new
        raddata.cond_rad_term[i] += dX_total - dX_cond
        return
    end

    for _ in 1:iters
        r_1 = sqrt(r2_new)
        r3    = r2_new * r_1
        rad_new = radcoeff * calc_cond_rad_term(r_1, z,constants, raddata, absliq_r_interp)
        g_mid = r2_new - R2 - timestep * (c * (S1 + (b - A * r2_new) / r3) + rad_new)
        if g_lo * g_mid <= 0
            hi2 = r2_new;  g_hi = g_mid
        else
            lo2 = r2_new;  g_lo = g_mid
        end
        (hi2 - lo2) < 1e-4 * r2_new && break
        r2_new = (lo2 + hi2) * FT(0.5)
    end
    r2_new = (lo2 + hi2) * FT(0.5)

    r_new  = sqrt(r2_new)
    r_new  = max(r_new, dry_r)
    r2_new = r_new * r_new
    X_new  = FT(4π/3) * r2_new*r_new

    dX_total = X_new - X_old
    dX_cond  = iszero(F0) ? dX_total : dX_total * Fcond0 / F0

    droplets.X[i] = X_new
    raddata.cond_rad_term[i] += dX_total - dX_cond
    return
end

"""
    drkohler_activated(R, T, Senv, timestep)

Compute the rate of change of droplet radius for activated droplets, neglecting
    differences in solutes.

# Arguments
- `R`: Radius of the droplets.
- `T`: Temperature.
- `Senv`: Environmental supersaturation.
- `timestep`: simulation Time step.

# Returns
dr: rate of change of droplet radius.

"""
function drkohler_activated(R,T,Senv,constants,timestep)
    # S = Senv/esat(T)
    denom = (FK(T, constants)+FD(T, constants))
    dr = (Senv-1) ./(denom.*R)
    return dr#R + dr*timestep > 0 ? dr : -R/timestep
end

"""
    dXkohler_function_of_radius(R, M, m, T, qv, P, timestep)

Calculate the change in droplet volume due to condensation using the Kohler equation.

# Arguments
- `R`: Droplet radius
- `M`: Molecular weight of the droplet substance
- `m`: Molecular weight of the dry air
- `T`: Temperature
- `qv`: Water vapor mixing ratio
- `P`: Pressure
- `timestep`: Time step

# Returns
- `dX`: Change in droplet volume

"""
function dXkohler_function_of_radius(R, M, m, T, qv, P,constants, timestep)
    dX = 4 * π * R^2 * drkohler(R, M, m, T, qv, P)
    return dX
end

"""
    dXkohler_function_of_radius(R, M, m, T, Senv, timestep)

Calculate the change in droplet volume due to condensation using the Kohler equation.

# Arguments
- `R`: Droplet radius
- `M`: Molecular weight of the droplet substance
- `m`: Molecular weight of the dry air
- `T`: Temperature
- `Senv`: Saturation of the environment
- `timestep`: Time step

# Returns
- `dX`: Change in droplet mass

"""
function dXkohler_function_of_radius(R, M, m, T, Senv,constants, timestep)
    dX = 4 * π * R^2 * drkohler(R, M, m, T, Senv,constants, timestep)
    return dX
end

"""
    dXkohler_function_of_radius_activated(R, T, Senv, timestep)

Calculate the change in droplet volume due to condensation using the Kohler equation for activated droplets,
    neglecting solute.

# Arguments
- `R`: Droplet radius
- `T`: Temperature
- `Senv`: Saturation of the environment
- `timestep`: Time step

# Returns
- `dX`: Change in droplet mass

"""
function dXkohler_function_of_radius_activated(R, T, Senv,constants, timestep)
    dX = 4 * π * R^2 * drkohler_activated(R, T, Senv,constants, timestep)
    return dX
end

# function FK_vent(T) #does k contain these already?
#     Fk = (constants.L ./(constants.Rv.*T) .-1).*(constants.L*constants.ρl)./(constants.k .*T)
#     return Fk/(Fa*Fq)
# end

K(T) = (2.38+0.00703(T - constants.T0))*1e-2
const β_diff = 0.04
β(T) = β_diff * exp(-(T-constants.T0)/85)
const α_diff = 0.7
D(T,p) = 2.11e−5 * (T/constants.T0)*1.94*(constants.P0/p)

function F_α(R,T,ρ_air,Re,K_th,constants)
    l_q = (2π/(constants.Rd*T))^0.5*K(T)*f_q(Re,ρ_air,T)/(ρ_air*constants.Cp_air*2*α_diff*(2-α_diff)^(-1))
    # l_q = K_th * f_q(Re,ρ_air,T) / (ρ_air*constants.Cp_air*1/4 * α_diff* (8*(constants.Rd*T/pi))^0.5)
    return R/(R+ l_q)
end

function F_β(R,T,ρ_air,P,Dv,Re,constants)
    l_m =  (2π/(constants.Rv*T))^0.5 * Dv * f_m(Re,ρ_air,T,P)/(2*β(T)*(2-β(T))^(-1))
    return R/(R + l_m)
end

function fvmol(Re,ScPr,constants)
    re_half_scpr_third = Re^0.5 * ScPr^(1/3)
    if re_half_scpr_third < 1.4
        # fvmol = 1+0.108*(re_half_scpr_third^2)
        fvmol = 1+0.14*(re_half_scpr_third^2)
    else
        # fvmol = 0.78 + 0.308*re_half_scpr_third
        fvmol = 0.86 + 0.28*re_half_scpr_third
    end
    return fvmol
end


function Sc(ρ_air,T,P,constants)
    return η_air(T,ρ_air,constants)/(ρ_air*D(T,P))
end

function Pr(ρ_air,T,constants)
    return η_air(T,ρ_air,constants)*constants.Cp_air/K(T)
end
f_m(Re_,ρ_air,T,P,constants) = fvmol(Re_,Sc(ρ_air,T,P,constants))
f_q(Re_,ρ_air,T,constants) = fvmol(Re_,Pr(ρ_air,T,constants))

function Reynoldsnumber(R,ρ_air,T,constants)
    return 2*R*ρ_air*terminal_v(R)/η_air(T,ρ_air,constants)
end

function η_air(T,ρ_air,constants) #kinematic viscosity of air
    μ = constants.μ*(T/296.16)^1.5 * (T + 120)/(T+296.16) #Hall and Pruppracher 1976 term index
    return μ
end

function C_drag(Re_)
    return 24*(1 + 0.15*Re_^0.687)/Re_ + 0.42/(1 + 42500*Re_^(-1.16))
end


# function v_term(radius_m)
#     radius_cm = radius_m * 1e2
#     if radius_m < 30e-6
#         vt = 1.19e6 * radius_cm^2
#     elseif radius_m < 60e-6
#         vt = 8e3 * radius_cm
#     elseif radius_m < 2e-3
#         vt = 2.01e3 * radius_cm^0.5
#     else
#         vt = 2.01e3 * 2e-3^0.5
#     end 
#     return vt * 1e-2 #convert back to m/s
# end



function dM_dt(R,T,Senv,p_air,ρ_air,η,ε_r,Fs,ε_0,constants) #Zeng et al 2018
    Dv = D(T,p_air)
    Re_ = Reynoldsnumber(R,ρ_air,T)
    K_th = K(T)
    Fafq = F_α(R,T,ρ_air,Re_,K_th)*f_q(Re_,ρ_air,T)
    Fbfm = F_β(R,T,ρ_air,p_air,Dv,Re_)*f_m(Re_,ρ_air,T,p_air)

    lvrtm1 = (constants.L/(constants.Rv*T)) - 1
    term2 = 1 - η*ε_r - Fs/(ε_0*constants.σSB*T^4)
    term3 = R*(ε_0*constants.σSB*T^3)/(K_th*Fafq)

    Hwc = 1 - lvrtm1*term2*term3
    Aw = lvrtm1 * constants.L/(Fafq*K_th*T)
    Bw = constants.Rv*T/(Dv*esat(T)*Fbfm)

    HAB = Hwc*Aw +Bw
    dm_term = 4pi*R*(Senv - 1)/HAB
    dm_rad_term = term2*lvrtm1 * 4pi*R^2*ε_0*constants.σSB*T^3/(K_th*Fafq*HAB)

    return dm_term + dm_rad_term
end


"""
    dq_liq_cond(R, M, m, T, Senv, timestep, ρ_air)

Calculate the change in q, liquid mixing ratio due to condensation of droplets,
using droplet solute information.

# Arguments
- `R`: Droplet radius, meters
- `M`: Molecular weight of the droplet
- `m`: Molecular weight of the solute
- `T`: Temperature
- `Senv`: Environmental saturation
- `timestep`: Time step
- `ρ_air`: Density of air

# Returns
- `dql`: Change in liquid water mass

"""
function dq_liq_cond(R, M, m, T, Senv,constants, timestep, ρ_air)
    dX = dXkohler_function_of_radius(R, M, m, T, Senv, timestep)
    dql = sum(dX .* ξ .* constants.ρl / ρ_air)
    return dql
end

"""
    dq_liq_cond_activated(R, T, Senv, timestep, ρ_air)
    dq_liq_cond_activated(R, M, m, T, Senv, timestep, ρ_air)

Calculate the change in liquid water mass due to condensation of activated droplets.

# Arguments
- `R`: Droplet radius, meters
- `T`: Temperature
- `Senv`: Environmental saturation
- `timestep`: Time step
- `ρ_air`: Density of air

# Returns
- `dql`: Change in liquid water mass

"""
function dq_liq_cond_activated(R, T, Senv,constants, timestep, ρ_air)
    dX = dXkohler_function_of_radius_activated(R, T, Senv, timestep)
    dql = sum(dX .* ξ .* constants.ρl / ρ_air)
    return dql
end

function dX_droplets!(X,dry_r3, kappa, qv, T, P, constants,dt)
    # Calculate the change in droplet volume due to condensation
    R = volume_to_radius.(X)  # Convert volume to radius
    dX = 4 * π * R^2 * drkappakohler(R,dry_r3,kappa,Senv,constants)
    return dX
end


function set_X_crit!(droplets,i,kappa,T)
    # Set the critical volume for condensation
    b = kappa * droplets.dry_r3[i]
    a = akk(T)
    R_crit = sqrt(3 * b / a) / 2
    droplets.X[i] = 0.98 * radius_to_volume(R_crit)  # Convert radius to volume
end



# function find_equilibrium_radius(droplets,drop_idx, kappa, T, S_env,constants; max_iter=100, tol=1e-12)
#     #println("Finding equilibrium radius for droplet $drop_idx with S_env=$S_env, T=$T")
#     dry_r3 = droplets.dry_r3[drop_idx]
#     dry_r = (dry_r3)^(1/3)
#     R_guess = max(dry_r * 2, 1e-8)
    
#     for i in 1:4
        
#         kelvin_term = akk(T) / R_guess
#         kappa_term = kappa * dry_r3 / (R_guess^3)
        

#         f = S_env - 1.0 - kelvin_term + kappa_term
        
#         if abs(f) < tol
#             droplets.X[drop_idx] = radius_to_volume(R_guess)
#             #println("Equilibrium radius converged for droplet $drop_idx: dry_r= $dry_r R = $R_guess m")
#             return
#         end

#         dS_dR = (-kelvin_term/R_guess + 3*kappa_term/R_guess)
#         R_new = R_guess - f / dS_dR
          
#         R_new = max(R_new, dry_r * 1.01)  # Must be larger than dry radius
#         R_new = min(R_new, 1e-3)          # Cap at 1mm
        
#         if abs(R_new - R_guess) < tol
#             droplets.X[drop_idx] = radius_to_volume(R_new)
#             #println("Equilibrium radius converged for droplet $drop_idx: dry_r= $dry_r R = $R_new m")
#             return
#         end
        
#         R_guess = R_new
#     end
    
#     return
# end

function find_equilibrium_radius(droplets,drop_idx, kappa, T, S_env,constants; max_iter=100, tol=1e-12)
    #println("Finding equilibrium radius for droplet $drop_idx with S_env=$S_env, T=$T")
    if S_env >=1.0
        error("S_env must be less than 1.0 for equilibrium radius calculation.")
    end
    dry_r3 = droplets.dry_r3[drop_idx]
    dry_r = (dry_r3)^(1/3)
    Sm1 = S_env - 1.0
    

    #use rcrit and do bisection
    a_term = akk(T)
    b_term = kappa * dry_r3
    lo_r = dry_r
    hi_r = sqrt(3 * b_term / a_term)
    r_guess = (lo_r + hi_r) / 2
    g_hi = Sm1 - a_term/hi_r + b_term/hi_r^3
    g_lo = Sm1 - a_term/lo_r + b_term/lo_r^3
    if g_hi * g_lo >= 0
        r_guess = max(r_guess, dry_r)
        droplets.X[drop_idx] = radius_to_volume(dry_r)
        return
    end

    for _ in 1:max_iter
       
        g_mid = Sm1 - a_term/r_guess + b_term/r_guess^3
        if g_lo * g_mid <= 0
            hi_r = r_guess;  g_hi = g_mid
        else
            lo_r = r_guess;  g_lo = g_mid
        end
        (hi_r - lo_r) < 1e-15 * r_guess && break
        r_guess = (lo_r + hi_r) * (0.5)
    end
    r_guess = (lo_r + hi_r) * (0.5)

    r_new = max(dry_r,r_guess)
    droplets.X[drop_idx] = radius_to_volume(r_new)
    return
end




# function condensation_rhs!(du, u, p, t)
#     du .= 0#zero.(u)  # Initialize du to zero, preserving the structure and types of u
#     FT = eltype(u.lnR)
#     lnR = u.lnR
#     qvap_col = u.qvap
#     T_col = u.T
#     nz,drops,grid,constants,condsettings,spatialsettings = p
#     R = exp.(lnR)
#     T_v = T_virtual.(T_col, qvap_col)
#     Rd = constants.Rd
#     gridV = spatialsettings.z_grid_height * spatialsettings.area_per_grid

#     for k in 1:nz
#         P::FT = grid.P[k]
#         T::FT  = T_col[k]
#         qv::FT = qvap_col[k]
#         Tv::FT = T_v[k]
#         ρ_air::FT = P / Rd / Tv
#         S_env::FT = sat(qv, P) / esat(T)
#         #find all indexes where drops.cell_id == k
#         R_idx = findall(drops.cell_id .== k)
#         if !isempty(R_idx)
#             dR = drkappakohler.(R[R_idx],drops.dry_r3[R_idx],condsettings.kappa,T,S_env,constants, t)
#             dX = 4 * π .* R[R_idx].^2 .* dR

#             dqvap = - sum(dX .* drops.ξ[R_idx] .* constants.ρl / ρ_air) / gridV
#             dT = - dqvap * constants.L / constants.Cp_air

#             du.lnR[R_idx] .= FT.(dR ./ R[R_idx])
#             du.qvap[k] = FT(dqvap)
#             du.T[k] = FT(dT)
#         end
#     end
# end


limitsat(::DynON, S) = min(S, 1.01)
limitsat(::DynOFF, S) = S

function apply_thermo_feedback!(::DynON, state, z, dqv, constants)
    state.T_tmp[z] -= dqv * constants.L / constants.Cp_air
    state.θ[z] = theta_from_T(state.T_tmp[z], state.P[z],state.qv[z], constants)
    state.ρ[z] = ρ_ideal_gas(state.P[z], state.T_tmp[z], state.qv[z], constants)
end
apply_thermo_feedback!(::DynOFF, state, z, dqv, constants) = nothing
function condensation_time_step_spatial!(::DynOFF,droplets::droplet_attributes{FT}, state::states, nz::Int, Δtg::FT, conddata::cdat, constants::Constants{FT},
    condsettings::condensation_settings{FT}, spatialsettings::spatial_settings{FT}, raddata::Rad,scmsettings::scm_settings{FT},absliq_r_interp, m_t::FT) where {FT<:AbstractFloat,Rad,states,cdat}
end
function condensation_time_step_spatial!(::DynON,droplets::droplet_attributes{FT}, state::states, nz::Int, Δtg::FT, conddata::cdat, constants::Constants{FT},
    condsettings::condensation_settings{FT}, spatialsettings::spatial_settings{FT}, raddata::Rad,scmsettings::scm_settings{FT},absliq_r_interp, m_t::FT) where {FT<:AbstractFloat,Rad,states,cdat}
    # Calculate condensation time step for each droplet
    # Ns = spatialsettings.Nz
    # FT = eltype(droplets.X)

    
    for z in 1:nz
        isempty(droplets.grid_range[z]) && continue
        T_k = state.T_tmp[z]
        conddata.S_env_cell[z] = limitsat(scmsettings.spinupsaturation, sat(state.qv[z], state.P[z]) / esat(T_k))
        fk = FK(T_k, constants)
        fd = FD(T_k, constants)
        conddata.c_cell[z] = FT(2) / (fk + fd)
        conddata.radcoeff_cell[z] = (fk / (constants.L * constants.ρl)) / (fk + fd)
    end

    lo_cell = findfirst(!isempty, droplets.grid_range)
    if lo_cell !== nothing
        hi_cell = findlast(!isempty, droplets.grid_range)
        i_lo = first(droplets.grid_range[lo_cell])
        i_hi = last(droplets.grid_range[hi_cell])
        Threads.@threads for i in i_lo:i_hi
            idrop = droplets.I[i]
            if droplets.ξ[idrop] == 0
                conddata.droplet_vol_change[idrop] = zero(FT)
                continue
            end
            z = droplets.cell_id[idrop]
            T_k      = state.T_tmp[z]
            S_env    = conddata.S_env_cell[z]
            c        = conddata.c_cell[z]
            radcoeff = conddata.radcoeff_cell[z]
            X_before = droplets.X[idrop]
            dXkappakohler_bisection(scmsettings.REM, droplets, idrop, condsettings.kappa, T_k, S_env, c, radcoeff, constants, raddata, absliq_r_interp, Δtg, 20)
            conddata.droplet_vol_change[idrop] = (droplets.X[idrop] - X_before) * droplets.ξ[idrop]
        end
    end

    for z in 1:nz
        isempty(droplets.grid_range[z]) && continue
        k = @view droplets.I[droplets.grid_range[z]]

        ΔVdrop = zero(FT)
        for idrop in k
            ΔVdrop += conddata.droplet_vol_change[idrop]
        end

        inv_vol = constants.ρl / (state.ρ[z] * spatialsettings.area_per_grid * spatialsettings.z_grid_height)
        dqv = -ΔVdrop * inv_vol
        conddata.condensation_src[z] -= dqv
        state.qv[z] += dqv
        apply_thermo_feedback!(scmsettings.thermo_feedback, state, z, dqv, constants)

        if scmsettings.REM == DynON()
            conddata.condensation_rad_net[z] += inv_vol * sum(i -> raddata.cond_rad_term[i] * droplets.ξ[i], k)
            conddata.condensation_rad_abs[z] += inv_vol * sum(i -> abs(raddata.cond_rad_term[i]) * droplets.ξ[i], k)
        end
    end
    
    # dqvap = - vol_change .* constants.ρl ./ (state.ρ .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
    # state.qv .+= dqvap
    # state.T_tmp .-= dqvap * constants.L ./ (constants.Cp_air)
    # state.θ .= theta_from_T(state.T_tmp, state.P, constants)
    # state.ρ .= ρ_ideal_gas(state.P, state.T_tmp, state.qv, constants)

    return
end
