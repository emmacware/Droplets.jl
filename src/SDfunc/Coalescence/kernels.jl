

#---------------------------------------------------------
# KERNELS
#---------------------------------------------------------
export terminal_v,hydrodynamic,golovin,long1974, terminal_v_X
include("hall_davis.jl")

# terminal velocity of droplets (this should move..)
"""
    terminal_v(r::FT) where FT<:AbstractFloat

Compute the terminal velocity of a droplet of radius r, using tables from
Gunn and Kinzer, (1949), https://doi.org/10.1175/1520-0469(1949)006<0243:TTVOFF>2.0.CO;2

# Arguments
- `r::FT`: The radius of the droplet in meters.

# Returns
The terminal velocity of the droplet in meters/second.

"""
const _tv_d_table = [0.078,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,
    2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8] ./ 10  # diameter in cm

const _tv_v_table = [18,27,72,117,162,206,247,287,327,367,403,464,517,565,609,649,690,
    727,757,782,806,826,844,860,872,883,892,898,903,907,909,912,914,916,917]  # velocity in cm/s

const _tv_interp = linear_interpolation(_tv_d_table, _tv_v_table, extrapolation_bc=Line())

function terminal_v(r::FT)::FT where FT<:AbstractFloat  # terminal velocity
    d_cm = 2 * r * 100  # radius in meters → diameter in cm
    if d_cm < 0.0078
        # return FT(1.2e7 * (r * 100)^2 / 100)  # Stokes regime (note: 10e6 == 1e7)
        return 1.2e6 * (d_cm/2)^2 * 1e-2 #from Adele slides
    elseif d_cm > 0.58
        return FT(_tv_interp(0.58) / 100)  # Gunn and Kinzer, need?
    else
        return FT(_tv_interp(d_cm) / 100)  # cm/s → m/s
    end
end

function terminal_v_X(X::FT)::FT where FT<:AbstractFloat  # terminal velocity
    d_cm = 2 * volume_to_radius(X) * 100  # radius in meters → diameter in cm
    if d_cm < 0.0078
        # return FT(1.2e7 * (r * 100)^2 / 100)  # Stokes regime (note: 10e6 == 1e7)
        return 1.2e6 * (d_cm/2)^2 * 1e-2 #from Adele slides
    elseif d_cm > 0.58
        return FT(_tv_interp(0.58) / 100)  # Gunn and Kinzer, need?
    else
        return FT(_tv_interp(d_cm) / 100)  # cm/s → m/s
    end
end


# #not working correctly:
# #collision efficiency function
# function collision_efficiency(R1::FT,R2::FT)::FT where FT<:AbstractFloat
#     #Parameterization from Berry 1967
#     #https://doi.org/10.1175/1520-0469(1967)024<0688:CDGBC>2.0.CO;2
#     r = max(R1,R2)*1e6
#     rs = min(R1,R2)*1e6

#     p = rs/r
#     D = (-27)/(r^1.65)
#     E = (-58)/(r^1.9)
#     F = (15/r)^4 +1.13 
#     G = (16.7/r)^8 +1 +0.004*r

#     Y = 1+p+D/(p^F)+E/((1-p)^G)
#     Y < 0 && (Y = 0)
#     Y < 1 && (Y = 1)
#     # Y > 1 && error("Collision efficiency cannot be greater than 1. Check the input radii.")

#     return Y
# end

function collision_efficiency(R1::FT, R2::FT)::FT where FT <: AbstractFloat
    r_max = FT(1100)

    # meters → micrometers
    r1 = R1 * FT(1e6)
    r2 = R2 * FT(1e6)

    r1 = min(r1, r_max - FT(1e-6))
    r2 = min(r2, r_max - FT(1e-6))

    # lower grid point and spacing for each radius
    if r1 >= FT(100)
        x0 = floor(r1 / FT(10)) * FT(10);  dx = FT(10)
    else
        x0 = floor(r1);                     dx = FT(1)
    end

    if r2 >= FT(100)
        y0 = floor(r2 / FT(10)) * FT(10);  dy = FT(10)
    else
        y0 = floor(r2);                     dy = FT(1)
    end

    x1 = x0 + dx
    y1 = y0 + dy

    # radius (μm) → 0-based table index
    idx(R) = R <= FT(100) ? Int(R) : 100 + Int((R - FT(100)) / FT(10))

    # (i, j) → 1-based index into triangular-packed vector
    tidx(i, j) = let ii = max(i, j), jj = min(i, j)
        ii * (ii + 1) ÷ 2 + jj + 1
    end

    i00 = tidx(idx(x0), idx(y0))
    i10 = tidx(idx(x1), idx(y0))
    i01 = tidx(idx(x0), idx(y1))
    i11 = tidx(idx(x1), idx(y1))

    # bilinear interpolation
    wx = r1 - x0
    wy = r2 - y0

    return (
        hall_davis_efficiencies[i00] * (dx - wx) * (dy - wy) +
        hall_davis_efficiencies[i10] *       wx  * (dy - wy) +
        hall_davis_efficiencies[i01] * (dx - wx) *       wy  +
        hall_davis_efficiencies[i11] *       wx  *       wy
    ) / (dx * dy)
end

# function collision_efficiency(r1::FT, r2::FT) where {FT<:AbstractFloat}

#     # ------------------------------------------------------------
#     # Ensure ordering:
#     # R = collector radius (larger)
#     # r = collected radius (smaller)
#     # ------------------------------------------------------------
#     R = max(r1, r2)
#     r = min(r1, r2)

#     # convert to microns
#     Rμ = R * FT(1e6)
#     rμ = r * FT(1e6)

#     # radius ratio
#     p = rμ / Rμ

#     # ------------------------------------------------------------
#     # Tiny droplets: Davis/Klett-style suppression
#     # ------------------------------------------------------------
#     if Rμ < 20

#         # strong hydrodynamic deflection
#         Ec = FT(0.001) * (Rμ / 10)^2 * p^1.5

#         return clamp(Ec, zero(FT), one(FT))
#     end

#     # ------------------------------------------------------------
#     # Intermediate cloud droplets: Hall-like growth
#     # ------------------------------------------------------------
#     if Rμ < 100

#         # empirical smooth growth
#         Ec = FT(0.01) *
#              (Rμ / 20)^1.8 *
#              p^0.7

#         return clamp(Ec, zero(FT), one(FT))
#     end

#     # ------------------------------------------------------------
#     # Rain-drop regime
#     # ------------------------------------------------------------
#     Ec = FT(0.8) + FT(0.2) * tanh((Rμ - 100) / 50)

#     return clamp(Ec, zero(FT), one(FT))
# end


# function collision_efficiency(R1::FT,R2::FT)::FT where FT<:AbstractFloat
#     #Parameterization from Simmel et al., 2002
#     Rj_cm = max(R1,R2)*1e4
#     Rk_cm = min(R1,R2)*1e4
#     Ecoal = 1
#     if Rk_cm < 50e-4
#         Ecol = max(4.5*1e4*Rk_cm^2(1-3e-4/Rj_cm),1e-3)
#     else
#         Ecol = 1
#     end
#     return Ecol
# end
#---------------------------------------------------------
# Coalescence Kernels
#---------------------------------------------------------
"""
    hydrodynamic(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat

Hydronynamic kernel function for droplet coalescence, used to calculate the probabilities of two 
    droplets colliding.
    K(r,r') = π * (r+r')^2 * |v(r)-v(r')| * E(r,r')
    where E(r,r') is the collision efficiency function, (currently not implemented), and v is the terminal velocity of 
    a droplet of radius size r.

# Arguments
- `droplets::droplet_attributes`: The attributes of the droplets.
- `(j,k)::Tuple{Int,Int}`: The indices of the droplets.
- `settings::coag_settings{FT}`: The coagulation settings.

"""
@inline function hydrodynamic(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    Rj, Rk = volume_to_radius(droplets.X[j]), volume_to_radius(droplets.X[k])
    if settings.hydrodynamic_collision_eff_func == true
        E = collision_efficiency(Rj,Rk)
    else    
        E = 1
    end
    Rsum = (Rj+Rk) # sum of radius is in meters
    vdif = abs(terminal_v(Rj)-terminal_v(Rk))
    return E *π * Rsum^2 *vdif
end

"""
    golovin(droplets, (j, k), settings)

Golovin (or additive) kernel function for droplet coalescence, used to calculate the probabilities of two droplets colliding.
    K(x,x') = b(x+x'),where b is the Golovin kernel coefficient, and x is the volume of the droplet.
Taken from Golovin (1963), this kernel has an analytic solution for the Smulochowski Coagulation Equation.

# Arguments
- `droplets`: droplet attributes
- `(j, k)`: indices of the droplets
- `settings`: coagulation settings

"""
@inline function golovin(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    return settings.golovin_kernel_coeff *(droplets.X[j] + droplets.X[k])# Xsum
end

"""
    long1974(droplets, (j, k), settings)

Piecewise collision kernel from Eq. (11) in Long (1974),
https://doi.org/10.1175/1520-0469(1974)031<1040:STTDCE>2.0.CO;2

    K(x, x') = sq_coeff * (x_lg² + x_sm²)   if r_lg < r_thresh
    K(x, x') = lin_coeff * (x_lg + x_sm)     if r_lg ≥ r_thresh

where x is droplet volume and r_lg is the radius of the larger droplet.
Default parameters: lin_coeff = 5.78e3 s⁻¹, sq_coeff = 9.44e15 m⁻³ s⁻¹, r_thresh = 5e-5 m.
"""
@inline function long1974(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    v_lg = max(droplets.X[j], droplets.X[k])
    v_sm = min(droplets.X[j], droplets.X[k])
    r_lg = volume_to_radius(v_lg)
    if r_lg < settings.long_r_thresh
        return settings.long_sq_coeff * (v_lg^2 + v_sm^2)
    else
        return settings.long_lin_coeff * (v_lg + v_sm)
    end
end
