using Revise
using Droplets
using Distributions
# using DifferentialEquations

using OrdinaryDiffEq

using Plots
# using StatsPlots
# using ComponentArrays
using JLD2
using Interpolations
using Random
# using OrdinaryDiffEqBDF
include("../radiation_call.jl")
include("../initial_state_functions.jl")
const FT = Float64



# Numerical settings
const Z_max = 1500.0 #m
const dz = 10.0 #m
const nz = Int(Z_max/dz)
const dt = 1.0 #s
const dt_coag = 0.1
const dt_cond = 0.1
const dt_rad = 1.0 #s,
const Ns_per_grid = 64
seed = 42
Random.seed!(seed)
const t_max = 3600*6#s
const t_output = 60*6 #s



#DyCOMS-II RF02 case specifications
12
const P_surface = 101780.0 # Pa, Ackerman et al., 2009
# From Wyant et al., 2007
const inv_height = 795.0  # m
const ρ_inv = 1.12 # kg/m^3, air density @ inversion height
const D = 3.75e-6 # s^-1, constant Divergence
const u_star = 0.25 # m/s, friction velocity
const geostrophic_u(z)::FT = 3 + 4.3*(z/1000) # m/s
const geostrophic_v(z)::FT = -9 + 5.6*z/1000        # m/s
const prescribed_w(z)::FT = min(-D*z,0) # m/s
θl_initial(z) = z <= inv_height ? 288.3 : 295.0 + (z - inv_height)^(1/3) #boundary layer liquid water potential temperature
qt_initial(z) = z <= inv_height ? 9.45e-3 : 5e-3 - 3e-3(1-exp(-(z - inv_height)/500)) # q total
θl_init_inv(inv_height,z) = z <= inv_height ? 288.3 : 295.0 + (z - inv_height)^(1/3) #boundary layer liquid water potential temperature
qt_init_inv(inv_height,z) = z <= inv_height ? 9.45e-3 : 5e-3 - 3e-3(1-exp(-(z - inv_height)/500)) # q total
const surface_latent_heat_flux = 93.0 # W/m^2
const surface_sensible_heat_flux = 16.0 # W/m^2
#Bimodal Aerosol Distribution, Ackerman et al., 2009
#ammonium sulfate, Dziekan et al 2019 (from Petters and Kreidenweis (2007))
kappa_ammonium_sulfate = 0.61
ammonium_sulfate_density = 1.78e3 # kg/m^3 
n1 = 125 * ccm_to_cm # m^-3, number concentration of mode 1
n2 = 65 * ccm_to_cm # m^-3, number concentration of mode 2
n0 = n1 + n2 # m^-3, total number concentration
m1,σ1 = (0.011e-6), log(1.2) # mode 1
m2,σ2 = (0.06e-6), log(1.7) # mode 2
# dist = MixtureModel(LogNormal, [(log(m1) + σ1^2,σ1),(log(m2) + σ2^2,σ2)], [n1/n0, n2/n0])
dist = MixtureModel(LogNormal, [(log(m1), σ1), (log(m2), σ2)], [n1/n0, n2/n0])

# inv_idx = findfirst(k -> grid.centers_z[k] >= inv_height, 1:nz)




#Settings Structs
spatialsettings = spatial_settings_1d{FT}(Nz=nz, Z_max=Z_max,dt=dt, t_max=t_max, dt_output=t_output,area_per_grid=1.0,
    weighted_droplet_allocation = true, # uniform Ns/nz per cell instead of qv-weighted seeding
)

coagsettings = coag_settings{FT}(Ns=Ns_per_grid*nz,ΔV=dz*spatialsettings.area_per_grid, n0=n0,Δt=dt_coag,kernel=hydrodynamic,hydrodynamic_collision_eff_func=true)

condensationsettings = condensation_settings{FT}(kappa=kappa_ammonium_sulfate,ρ_solute = ammonium_sulfate_density,Δt=dt_cond)

mpdatasettings = mpdata_settings_1d(nz,nonoscillatory=true, vertical_boundary_condition=NoFlux(),infinite_gauge=true,
    # thermo_variable = ThetalQtVar(), # advect (θ_l,qt) instead of (θ,qv)
)

tkesettings = tke_settings{FT}(u_star=u_star, geostrophic_u=geostrophic_u, geostrophic_v=geostrophic_v,
    LHF=surface_latent_heat_flux, SHF=surface_sensible_heat_flux,
    mixing_length_scheme = BottMixing(), # or DeardorffMixing(); default EDMFXMixing()
    average_e_l_3pt = true, # use e[z] (and l[z]) alone in droplet diffusion instead of averaging [z-1,z,z+1]
    # droplet_diffusion_length_dz = true, # use the diagnosed mixing length l[z] instead of dz in droplet diffusion
    # density_weighted_diffusion = false, # mass-weight implicit_diffuse! (ρK face diffusivity) instead of plain ∂/∂z(K∂ϕ/∂z)
)

diagnosticsettings = diagnostic_settings{FT}()

#Default dynamics are all on
base_dynamics = (
    turbulence                  = DynON(),
    motion                      = DynON(),
    advection                   = DynON(),
    radiation                   = DynON(),
    condensation                = DynON(),
    turbulent_droplet_diffusion_on = DynON(),
    # droplet_diffusion_scheme    = WellMixedDropletDiffusion(), # has the well-mixed drift correction + no qv-threshold inversion wall (see nonlocal.jl); or OUDropletDiffusion() (no drift correction, has the wall), WeilDropletDiffusion(), VisserDropletDiffusion(), NoDropletDiffusion()
    keep_grid_filled            = DynON(), # disable donor-cell SD reseeding of thin layers
)

dyn_spinup = dynamic_settings(; base_dynamics...,
    settling         = DynOFF(),
    spinupsaturation = DynON(),
    coalescence       = DynOFF(),
)

dyn_main = dynamic_settings(; base_dynamics...,
    settling         = DynON(),
    spinupsaturation = DynOFF(),
    coalescence       = DynON(),
)

base_settings = (
    spatial            = spatialsettings,
    coagsettings       = coagsettings,
    condsettings       = condensationsettings,
    tkesettings        = tkesettings,
    mpdatasettings     = mpdatasettings,
    diagnosticsettings = diagnosticsettings,
    Δt                 = dt,
    n_cond             = round(Int, dt / dt_cond),
    n_coag             = round(Int, dt / dt_coag),
    spinup_time        = 3600.0,
    coag_threading     = Parallel(), # or Serial() for single-threaded
)

scmspinupsettings = scm_settings{FT}(; base_settings..., dynamics = dyn_spinup, n_rad = round(Int, dt / dt))
scmsettings        = scm_settings{FT}(; base_settings..., dynamics = dyn_main, n_rad = round(Int, dt_rad / dt))



# Build environment and droplets
scmgrid, droplets, scmdata = initialize_scm_environment(nz, dz, P_surface, θl_initial, qt_initial, prescribed_w, dist,
    scmsettings, constants;
    build_raddata = radiation_data)

#initialize with some turbulence
scmgrid.states.e .= 0.1

run_scm!(scmgrid, droplets, scmdata, constants, scmsettings, prescribed_w;
    scmspinupsettings = scmspinupsettings,
    )


ptime = plot!(plot_output_timeseries(scmgrid),plot_title="Droplets.jl DyCOMS")