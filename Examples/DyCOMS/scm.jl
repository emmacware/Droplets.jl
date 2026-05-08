using Revise
using Droplets
# using DifferentialEquations

using OrdinaryDiffEq

using Plots
using ComponentArrays
using JLD2
# using OrdinaryDiffEqBDF
include("radiation_call.jl")
include("initial_state_functions.jl")
include("forward_solve.jl")
include("prescribed_tke.jl")
FT = Float64





# Numerical settings
Z_max = 1500.0 #m
dz = 10.0 #m
nz = Int(Z_max/dz)
dt = 1.0 #s
Ns_per_grid =32
seed = 42
t_max = 3600*2#s
t_output = 60*5 #s



#DyCOMS-II RF02 case specifications

P_surface = 101780.0 # Pa, Ackerman et al., 2009
# From Wyant et al., 2007
inv_height = 795.0  # m
ρ_inv = 1.12 # kg/m^3, air density @ inversion height
D = 3.75e-6 # s^-1, constant Divergence
u_star = 0.25 # m/s, friction velocity
geostrophic_u(z) = 3 + 4*(z/1000) # m/s
geostrophic_v(z) = -9 + 5.6*z/1000        # m/s
prescribed_w(z) = min(-D*z,0) # m/s
θl_initial(z) = z <= inv_height ? 288.3 : 295.0 + (z - inv_height)^(1/3) #boundary layer liquid water potential temperature
qt_initial(z) = z <= inv_height ? 9.45e-3 : 5e-3 - 3e-3(1-exp(-(z - inv_height)/500)) # q total
surface_latent_heat_flux = 93 # W/m^2
surface_sensible_heat_flux = 16 # W/m^2
# surface_zonal_momentum_flux = geostrophic_u(0)*u_star^2 / sqrt(geostrophic_u(0)^2 + geostrophic_v(0)^2)
# surface_meridional_momentum_flux = geostrophic_v(0)*u_star^2 / sqrt(geostrophic_u(0)^2 + geostrophic_v(0)^2)

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
spatialsettings = spatial_settings_1d{FT}(Nz=nz, Z_max=Z_max,dt=dt, t_max=t_max, dt_output=t_output)
coagsettings = coag_settings{FT}(Ns=Ns_per_grid*nz,ΔV=dz*spatialsettings.area_per_grid, n0=n0,Δt=dt/10,kernel=hydrodynamic,hydrodynamic_collision_eff_func=true)
condensationsettings = condensation_settings{FT}(kappa=kappa_ammonium_sulfate,ρ_solute = ammonium_sulfate_density,Δt=dt)
mpdatasettings = mpdata_settings_1d(nz,nonoscillatory=true, vertical_boundary_condition=NoFlux(),infinite_gauge=true)
scmsettings = scm_settings{FT}(Δt=dt, surface_latent_heat_flux=surface_latent_heat_flux, surface_sensible_heat_flux=surface_sensible_heat_flux,
    turbulence_on=true, condensation_on=true, radiation_on=true, coalescence_on=true, settling=true,
    coag_threading=Serial(), scheme=none(), dt_cond=dt/10, dt_coag=dt/10, spinup_time=3600.0,
    rho_weighted_diffusion=false, turbulent_droplet_diffusion_on=false)
tkesettings = tke_settings{FT}(u_star=u_star, geostrophic_u=geostrophic_u, geostrophic_v=geostrophic_v,turbulence_scheme=MellorYamada(),mixing_length_function=Deardorff(), A_e=FT(0.2))
diagnosticsettings = diagnostic_settings()


#Read Radiation Files 
# lookup_lw, lookup_lw_cld, lookup_sw, lookup_sw_cld, idx_gases = read_radiation_tables()

# Create superdroplets
droplets = init_droplets_dycoms_scm(initial_aerosol_dist,coagsettings, spatialsettings)

# Create environmnent
grid = initialize_scm_environment(nz, dz, P_surface, θl_initial, qt_initial, geostrophic_u, geostrophic_v, prescribed_w,droplets,spatialsettings,condensationsettings,constants)
grid.states.e .= 0.06#1.0 #from Ackerman et al 2009
#add random perturbation -0.1 to 0.1 K to temperature field
grid.states.θ .+= rand(Uniform(-0.1,0.1), length(grid.states.θ))

# Create other needed data structures
coagdata = coagulation_run_spatial{FT}(nz, coagsettings.Ns,droplets)
# condensation_integrator = create_condensation_integrator(grid, droplets, condensationsettings, coagsettings, spatialsettings,constants)
conddata = condensation_data(FT,nz,coagsettings.Ns)
raddata = radiation_data(FT,length(droplets.X),grid,constants)
mpdatatmp = mpdata_tmp_1d(grid.states.qv, grid.faces_z)
Xinit= droplets.X .+ 0

# plot_env_profiles(grid)
CArad = raddata.CArad
R_array = range(CArad.lookup_lw_cld.bounds[1], CArad.lookup_lw_cld.bounds[2],
                length=CArad.lookup_lw_cld.dims[3])  # μm
const absliq_r_interp = ntuple(CArad.nband_lw) do ibnd
    ext, ssa, _ = LookUpTables.getview_liqdata(CArad.lookup_lw_cld, ibnd)
    absliq = collect(ext .* (1 .- ssa))
    linear_interpolation(R_array, absliq, extrapolation_bc=Flat())
end



for i in 1:Int(t_max/dt)
    if i*dt % 1 == 0
        println("Timestep: ", i)
    end
 
    single_column_timestep(grid,dt,droplets,coagsettings,spatialsettings,condensationsettings,
    coagdata,conddata,raddata,diagnosticsettings,prescribed_w, mpdatatmp, mpdatasettings,constants,scmsettings,tkesettings,i)

    if i % 1200 == 0
        plot_env_profiles(grid)
        title!("Timestep: $(i*dt) seconds")
    end
end

penv = plot_env_profiles(grid)
#put tkesettings.turbulence_scheme in the title
plot!(plot_output_timeseries(grid),plot_title="$(tkesettings.turbulence_scheme),mixing_length=$(tkesettings.mixing_length_function)")
