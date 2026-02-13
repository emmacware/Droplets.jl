##################################################
# Constants
##################################################

export Constants,constants
export kg_to_g, m_to_μm, volume_to_radius, radius_to_volume, ccm_to_cm

"""
    Constants{FT<:AbstractFloat}

A struct representing constants used in the Droplets package.

# Fields
- `FT`: The type of the constants, which must be a subtype of `AbstractFloat`.

- `k`: thermal conductivity (default: 0.024).
- `μ`: dynamic viscosity (default: 1.6e-5).
- `Dv`: diffusion coefficient of water vapor in air (default: 2.26e-5).
- `ρl`: density of liquid water (default: 1000.0).
- `Rd`: gas constant for dry air (default: 287.0).
- `Rv`: gas constant for water vapor (default: 461.0).
- `gconst`: gravitational constant (default: 9.8).
- `L`: latent heat of vaporization (default: 22.6e5).
- `Cp`: specific heat of dry air at constant pressure (default: 4181).
- `σSB`: Stefan-Boltzmann constant (default: 5.670374419e-8).

"""
Base.@kwdef struct Constants{FT<:AbstractFloat}
    λ = FT(1.5625)        # 
    κ = FT(1.5625)         # 
    k = FT(0.024)         # from clima K_therm
    μ = FT(1.6e-5)         # from clima ν_air
    Dv = FT(2.26e-5)        # from clima D_vapor water vapor in air m2/s???
    ρl = FT(1000.0)        # 
    Rd = FT(287.0)        # 
    Rv = FT(461.0)        # 
    gconst = FT(9.8)    # gravitational constant, m/s2
    L = FT(22.6e5)         # Latent Heat of Vaporization J/kg
    Cp_water = FT(4181)        # Specific Heat of Dry air at constant pressure J/kgK
    Cp_vapor = FT(1859)     # ClimaParams "isobaric_specific_heat_vapor"
    Cp_air = FT(1005)
    σSB = FT(5.670374419e-8) # Stefan-Boltzmann constant W⋅m−2⋅K−4
    T0 = FT(273.15)       # reference temperature K
    P0 = FT(100000.0)     # reference pressure Pa for potential temperature
    P_SLP = FT(101325.0)     # standard sea level pressure Pa
    # μ = FT(1.81e-5)         # Hall and Pruppracher 1976
    ϵ = FT(1.6080793637401138)           # molar mass ratio air and vapor
end

const constants = Constants{Float32}()

#Conversions
const kg_to_g = 1e3
const m_to_μm = 1e6
const ccm_to_cm = 1e6 # conversion from cm^-3 to m^-3

"""
    volume_to_radius(V::AbstractFloat)

radius based on spherical volume

# Arguments
- `V::AbstractFloat`: Spherical volume in meters^3.

# Returns
- The radius of the droplet in meters.

"""
@inline function volume_to_radius(V::AbstractFloat)
    return (3*V/(4*pi))^(1/3)
end

"""
    radius_to_volume(R::AbstractFloat)
Volume of sphere with radius R
# Arguments
- `R::AbstractFloat`: Radius in meters.
# Returns
- The volume of the droplet in meters^3.
"""
@inline function radius_to_volume(R::AbstractFloat)
    return 4/3*pi*R^3
end





