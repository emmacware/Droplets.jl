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
export drrad_term, dXkappakohler_newtonraphson
export calc_cond_rad_term
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
    Fk = (constants.L ./(constants.Rv.*T) .-1).*(constants.L*constants.ρl)./(constants.k .*T) 
    return Fk
end

function FD(T,constants)
    Fd = constants.ρl*constants.Rv .*(T)./(constants.Dv .* esat(T))
    return Fd
end

akk(T) = 3.3 * 10^(-7) / T #(m K/T)
bkk(m) = 4.3 * 2 / m / 1e6 #m^3 for NaCL

######################################################################
# Saturation Functions

#saturation vapor pressure
esat(T) = 100*6.1094.*exp.(17.625.*(T.-273)./((T.-273) .+243.04)) ## August–Roche–Magnus approximation, hPa to Pa, T in K, converted to C

#environmntal saturation based on the mixing ratio of water vapor to air
sat(qvarray,P) = qvarray.*P./(0.378.*qvarray .+ 0.622) ##qvarray is the specific humidity: mixing ratio of water vapor to moist air, Pa


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
    return R + dr * timestep > 0 ? dr : -R / timestep
end

function drkohler(R, M, m, T, Senv,constants, timestep)
    denom = (FK(T, constants) + FD(T, constants))
    dr = (Senv - 1 .- (akk(T) ./ R) .+ bkk.(m) .* M ./(R .^ 3)) ./(denom .* R)
    return R + dr * timestep > 0 ? dr : -R / timestep
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
        

    return R + dr * timestep > 0 ? dr : -R / timestep
end

function dXkappakohler_newtonraphson(R, dry_r3, kappa, T, Senv, constants,timestep, iters)
    # Fully implicit solve for R²_new via Newton-Raphson.
    # Find u = R²_new such that g(u) = u - R²_old - dt·F(u) = 0,
    # where F(u) = dR²/dt = 2(S-1 - A/√u + B/u^{3/2}) / denom.
    # g'(u) = 1 - dt·F'(u),  F'(u) = (A/r³ - 3B/r⁵) / denom  (r = √u)
    A = akk(T)
    b = kappa * dry_r3
    denom = FK(T, constants) + FD(T, constants)
    R2_old = R^2

    # Initial guess: explicit Euler step in R²
    F0 = 2 * (Senv - 1.0 - A/R + b/R^3) / denom
    R2 = max(R2_old + F0 * timestep, R2_old * 1e-4)

    for _ in 1:iters
        r  = sqrt(R2)
        F  = 2 * (Senv - 1.0 - A/r + b/r^3) / denom
        g  = R2 - R2_old - timestep * F
        R2r = R2 * r
        gp = 1.0 - timestep * (A/R2r - 3*b/(R2*R2r)) / denom
        δ  = g / gp
        R2 -= δ
        R2  = max(R2, R2_old * 1e-6)   # keep positive
        abs(δ) < 1e-18 * R2 && break
    end

    return max(radius_to_volume(sqrt(R2)), 4/3*pi *dry_r3)# - radius_to_volume(R)
end

function drrad_term(R, T, rad_term,constants, timestep)
    fk = FK(T, constants)
    denom = (fk + FD(T, constants))
    dr = rad_term * fk / (constants.L*constants.ρl * denom)
    return R + dr * timestep > 0 ? dr : -R / timestep
end

function calc_cond_rad_term(R,z,constants,raddata,absliq_r_interp)
    if R < 2.5e-6 || R > 4e-5
        return 0.0
    end
    condradterm = 0.0
    n_bnd = raddata.CArad.nband_lw
    for ibnd in 1:n_bnd
        abs_drop = absliq_r_interp[ibnd](R*1.0e6)
        Qa = abs_drop * radius_to_volume(R) * constants.ρl * kg_to_g / (pi *R^2) #grams
        condradterm += Qa * raddata.flux_net_droplet[z,ibnd]
    end
    return condradterm
end

function dXkappakohler_newtonraphson(droplets,i, kappa, T, Senv, constants,
    raddata,   # r -> dR/dt contribution from radiation
    absliq_r_interp,
    timestep, iters)
R = volume_to_radius(droplets.X[i])
dry_r3 = droplets.dry_r3[i]
z = droplets.cell_id[i]
A = akk(T)
b = kappa * dry_r3
fk = FK(T, constants)
denom = fk + FD(T, constants)
R2_old = R^2

r0 = R
# F0 = 2*(Senv - 1.0 - A/r0 + b/r0^3)/denom + 2*r0*drrad_term(r0,T,radterm,constants,timestep)/denom
F0 = 2*(Senv - 1.0 - A/r0 + b/r0^3)/denom + 2*r0*calc_cond_rad_term(r0,z,constants,raddata,absliq_r_interp) * fk / (constants.L*constants.ρl * denom)
R2 = max(R2_old + F0*timestep, R2_old*1e-4)

for _ in 1:iters
r  = sqrt(R2)
F_cond = 2*(Senv - 1.0 - A/r + b/r^3) / denom
F_rad  = 2*r * calc_cond_rad_term(r,z,constants,raddata,absliq_r_interp) * fk / (constants.L*constants.ρl * denom)      # re-evaluated at current r
g  = R2 - R2_old - timestep*(F_cond + F_rad)
gp = 1.0 - timestep*(A/(R2*r) - 3b/(R2^2*r)) / denom  # Jacobian: Köhler only, rad treated explicit
δ  = g / gp
R2 -= δ
R2  = max(R2, R2_old*1e-6)
abs(δ) < 1e-18*R2 && break
end

r_new   = sqrt(R2)
F_cond  = 2*(Senv - 1.0 - A/r_new + b/r_new^3) / denom
F_rad   = 2*r_new * calc_cond_rad_term(r_new,z,constants,raddata,absliq_r_interp) * fk / (constants.L*constants.ρl * denom)
f_total = F_cond + F_rad

X_new    = max(radius_to_volume(r_new), 4/3*π*dry_r3)
dX_total = X_new - radius_to_volume(R)
dX_cond  = iszero(f_total) ? dX_total : dX_total * F_cond / f_total
dX_rad   = dX_total - dX_cond

droplets.X[i] = X_new
raddata.cond_rad_term[i] = dX_rad
return #dX_cond, dX_rad
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
    return R + dr*timestep > 0 ? dr : -R/timestep
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

K(T) = (2.38+0.00703(T −constants.T0))*1e-2
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
    dX = 4 * π * R^2 * drkappakohler(R,dry_r3,kappa,Senv,constants,timestep)
    return dX
end


function set_X_crit!(droplets,i,kappa,T)
    # Set the critical volume for condensation
    b = kappa * droplets.dry_r3[i]
    a = akk(T)
    R_crit = sqrt(3 * b / a) / 2
    droplets.X[i] = 0.98 * radius_to_volume(R_crit)  # Convert radius to volume
end



function find_equilibrium_radius(droplets,drop_idx, kappa, T, S_env,constants; max_iter=100, tol=1e-12)
    #println("Finding equilibrium radius for droplet $drop_idx with S_env=$S_env, T=$T")
    dry_r3 = droplets.dry_r3[drop_idx]
    dry_r = (dry_r3)^(1/3)
    R_guess = max(dry_r * 2, 1e-8)
    
    for i in 1:4
        
        kelvin_term = akk(T) / R_guess
        kappa_term = kappa * dry_r3 / (R_guess^3)
        

        f = S_env - 1.0 - kelvin_term + kappa_term
        
        if abs(f) < tol
            droplets.X[drop_idx] = radius_to_volume(R_guess)
            #println("Equilibrium radius converged for droplet $drop_idx: dry_r= $dry_r R = $R_guess m")
            return
        end

        dS_dR = (-kelvin_term/R_guess + 3*kappa_term/R_guess)
        R_new = R_guess - f / dS_dR
          
        R_new = max(R_new, dry_r * 1.01)  # Must be larger than dry radius
        R_new = min(R_new, 1e-3)          # Cap at 1mm
        
        if abs(R_new - R_guess) < tol
            droplets.X[drop_idx] = radius_to_volume(R_new)
            #println("Equilibrium radius converged for droplet $drop_idx: dry_r= $dry_r R = $R_new m")
            return
        end
        
        R_guess = R_new
    end
    
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
