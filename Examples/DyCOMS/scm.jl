
using Droplets
using DifferentialEquations
using Plots
using ComponentArrays
using LSODA
include("initial_state_functions.jl")
include("forward_solve.jl")
include("radiation_call.jl")
FT = Float64



# Numerical settings
Z_max = 1500.0 #m
dz = 10.0 #m
nz = Int(Z_max/dz)
dt = 1.0 #s
Ns_per_grid = 16
seed = 30



#DyCOMS-II RF02 case specifications

P_surface = 101780.0 # Pa, Ackerman et al., 2009
# From Wyant et al., 2007
inv_height = 795.0  # m
ρ_inv = 1.12 # kg/m^3, air density @ inversion height
D = 3.75e-6 # s^-1, constant Divergence
u_star = 0.25 # m/s, friction velocity
prescribed_u(z) = 3 + 4*(z/1000) # m/s
prescribed_v(z) = -9 + 5.6*z/1000        # m/s
prescribed_w(z) = min(-D*z,0) # m/s
θl(z) = z <= inv_height ? 288.3 : 295.0 + (z - inv_height)^(1/3) #boundary layer liquid water potential temperature
qt(z) = z <= inv_height ? 9.45e-3 : 5e-3 - 3e-3(1-exp(-(z - inv_height)/500)) # q total
surface_latent_heat_flux = 93 # W/m^2
surface_sensible_heat_flux = 16 # W/m^2
surface_zonal_momentum_flux = prescribed_u(0)*u_star^2 / sqrt(prescribed_u(0)^2 + prescribed_v(0)^2)
surface_meridional_momentum_flux = prescribed_v(0)*u_star^2 / sqrt(prescribed_u(0)^2 + prescribed_v(0)^2)

#Bimodal Aerosol Distribution, Ackerman et al., 2009
#ammonium sulfate, Dziekan et al 2019 (from Petters and Kreidenweis (2007))
kappa_ammonium_sulfate = 0.61
ammonium_sulfate_density = 1.78e3 # kg/m^3 
n1 = 125 * ccm_to_cm # m^-3, number concentration of mode 1
n2 = 65 * ccm_to_cm # m^-3, number concentration of mode 2
n0 = n1 + n2 # m^-3, total number concentration
m1,σ1 = 0.011e-6, log(1.2) # mode 1
m2,σ2 = 0.06e-6, log(1.7) # mode 2
initial_aerosol_dist = MixtureModel(LogNormal, [(log(m1) + σ1^2,σ1),(log(m2) + σ2^2,σ2)], [n1/n0, n2/n0])


#Settings Structs
spatialsettings = spatial_settings_1d{FT}(Nz=nz, Z_max=Z_max)
coagsettings = coag_settings{FT}(Ns=Ns_per_grid*nz,ΔV=dz*spatialsettings.area_per_grid, n0=n0,Δt=dt)
condensationsettings = condensation_settings{FT}(kappa=kappa_ammonium_sulfate,ρ_solute = ammonium_sulfate_density,Δt=dt)
mpdatasettings = mpdata_settings_1d(nz, vertical_boundary_condition=NoFlux())
scmsettings = scm_settings{FT}(Δt=dt, surface_latent_heat_flux=surface_latent_heat_flux, surface_sensible_heat_flux=surface_sensible_heat_flux)
diagnosticsettings = diagnostic_settings()

#Read Radiation Files 
lookup_lw, lookup_lw_cld, lookup_sw, lookup_sw_cld, idx_gases = read_radiation_tables()

# Create environmnent
grid = initialize_scm_environment(nz, dz, P_surface, θl, qt, prescribed_u, prescribed_v, prescribed_w)
#to plot initial profiles: plot_env_profiles(grid)

# Create superdroplets
droplets = init_droplets_dycoms_scm(initial_aerosol_dist,coagsettings, spatialsettings)


# Create other needed data structures
coagdata = coagulation_run_spatial{FT}(nz, coagsettings.Ns,droplets)
condensation_integrator = create_condensation_integrator(grid, droplets, condensationsettings, coagsettings, spatialsettings,constants)
mpdata_tmp = mpdata_tmp_1d(grid.states.qv, grid.faces_z)


#set to critical radius to begin, need to change to eq radius
for k in 1:nz
    drop_idx = coagdata.I[droplets.grid_range[k]]
    set_X_crit!.(droplets,drop_idx,kappa_ammonium_sulfate,grid.states.T[k])
end


for i in 1:2
    if i % 1 == 0
        println("Timestep: ", i)
    end
    single_column_timestep(grid,dt,droplets,coagsettings,spatialsettings,condensationsettings,condensation_integrator,
    coagdata,diagnosticsettings,prescribed_w, mpdata_tmp, mpdatasettings,constants,scmsettings)
end

plot_env_profiles(grid)
