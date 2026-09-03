using Revise
using Droplets
# using DifferentialEquations

using OrdinaryDiffEq

using Plots
# using ComponentArrays
using JLD2
using Interpolations
using Random
# using OrdinaryDiffEqBDF
include("radiation_call.jl")
include("initial_state_functions.jl")
include("forward_solve.jl")
const FT = Float64



# Numerical settings
const Z_max = 1500.0 #m
const dz = 5.0 #m
const nz = Int(Z_max/dz)
const dt = 1.0 #s
const dt_coag = 0.1
const dt_cond = 0.1
const Ns_per_grid = 64
seed = 42
Random.seed!(seed)
const t_max = 3600*6#s
const t_output = 60*6 #s



#DyCOMS-II RF02 case specifications

const P_surface = 101780.0 # Pa Ackerman et al., 2009
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
spatialsettings = spatial_settings_1d{FT}(Nz=nz, Z_max=Z_max,dt=dt, t_max=t_max, dt_output=t_output,area_per_grid=5.0)

coagsettings = coag_settings{FT}(Ns=Ns_per_grid*nz,ΔV=dz*spatialsettings.area_per_grid, n0=n0,Δt=dt_coag,kernel=hydrodynamic,hydrodynamic_collision_eff_func=true)

condensationsettings = condensation_settings{FT}(kappa=kappa_ammonium_sulfate,ρ_solute = ammonium_sulfate_density,Δt=dt_cond)

mpdatasettings = mpdata_settings_1d(nz,nonoscillatory=true, vertical_boundary_condition=NoFlux(),infinite_gauge=true)

tkesettings = tke_settings{FT}(u_star=u_star, geostrophic_u=geostrophic_u, geostrophic_v=geostrophic_v,
    LHF=surface_latent_heat_flux, SHF=surface_sensible_heat_flux)
    
#Default dynamics are all on
base_scm = (
    Δt                          = dt,
    turbulence                  = DynON(),
    motion                      = DynON(),
    advection                   = DynON(),
    radiation                  = DynON(),
    condensation                = DynON(),
    n_cond                      = round(Int, dt / dt_cond),
    n_coag                      = round(Int, dt / dt_coag),
    spinup_time                 = 3600.0,
    turbulent_droplet_diffusion_on = DynON(),
)

scmspinupsettings = scm_settings{FT}(; base_scm...,
REM              = DynOFF(),
settling         = DynON(),
spinupsaturation = DynON(),
coalescence       = DynOFF(),
)
scmsettings = scm_settings{FT}(; base_scm...,
REM              = DynOFF(),
settling         = DynON(),
spinupsaturation = DynOFF(),
coalescence       = DynON(),
)
scmradsettings = scm_settings{FT}(; base_scm...,
REM              = DynON(),
settling         = DynON(),
spinupsaturation = DynOFF(),
coalescence       = DynON(),
)
diagnosticsettings = diagnostic_settings()



# Create environmnent
grid, droplets, coagdata,conddata,raddata,mpdatatmp,turbdata = initialize_scm_environment(
    nz, dz, P_surface, θl_initial, qt_initial, prescribed_w, dist,
    coagsettings,spatialsettings,condensationsettings,tkesettings,constants
    )

inv_idx = findfirst(k -> grid.centers_z[k] >= inv_height, 1:nz)
grid.states.e[1:inv_idx] .= 0.1#1e-6

CArad = raddata.CArad
R_array = range(CArad.lookup_lw_cld.bounds[1], CArad.lookup_lw_cld.bounds[2],
                length=CArad.lookup_lw_cld.dims[3])  # μm

const absliq_r_interp = ntuple(CArad.nband_lw) do ibnd
    ext, ssa, _ = LookUpTables.getview_liqdata(CArad.lookup_lw_cld, ibnd)
    absliq = collect(ext .* (1 .- ssa))
    linear_interpolation(R_array, absliq, extrapolation_bc=Flat())
end


spinup_step = round(Int, scmsettings.spinup_time / dt)

droplets_snapshots = Dict{Int, droplet_attributes}()

for i in 1:spinup_step
    if i*dt % 100 == 0
        println("Timestep: ", i*dt)
        droplets_snapshots[Int(div(i*dt,600))] = deepcopy(droplets)
    end

    single_column_timestep(grid,dt,droplets,coagsettings,spatialsettings,condensationsettings,
    coagdata,conddata,raddata,turbdata,
    diagnosticsettings,prescribed_w, mpdatatmp, mpdatasettings,constants,scmspinupsettings,tkesettings,absliq_r_interp,i)

end

for i in (spinup_step+1):Int(spatialsettings.t_max / dt)
    if i*dt % 100 == 0
        println("Timestep: ", i*dt)
        droplets_snapshots[Int(div(i*dt,600))] = deepcopy(droplets)
    end

    single_column_timestep(grid,dt,droplets,coagsettings,spatialsettings,condensationsettings,
    coagdata,conddata,raddata,turbdata,
    diagnosticsettings,prescribed_w, mpdatatmp, mpdatasettings,constants,scmsettings,tkesettings,absliq_r_interp,i)
end


using StatsBase,StatsPlots
num_bins = 100
radius_bins_edges = 10 .^range(log10(1e-7), log10(1e-3), length=num_bins+1)
runsettings = run_settings{FT}(num_bins=num_bins,radius_bins_edges=radius_bins_edges,binning_method=number_density,normalize_bins_dlnr=false)
# # seeds = 1    
# mids = (radius_bins_edges[1:end-1] .* radius_bins_edges[2:end]) .^ 0.5
time_idx = 20
x1 = droplets_snapshots[time_idx].grid_range[50][1]
x2 = droplets_snapshots[time_idx].grid_range[85][2]
range_dr = droplets_snapshots[time_idx].I[x1:x2]
r_um = (volume_to_radius.(droplets_snapshots[time_idx].X[range_dr]))*1e6
density(r_um, weights=Weights(droplets_snapshots[time_idx].ξ[range_dr]), bandwidth=0.1,normalize_weights=true,
        xlims=(1,80), 
#         # yscale=:log10,
#         # xscale=:log10,
        )
using KernelDensity
k = kde(r_um, weights=Weights(droplets_snapshots[time_idx].ξ[range_dr]), bandwidth=0.5)
plot(k.x, max.(k.density, 1e-10), yscale=:log10)


penv = plot_env_profiles(grid)
# #put tkesettings.turbulence_scheme in the title
ptime = plot!(plot_output_timeseries(grid),plot_title="REM: $(scmsettings.REM)")#, Coalescence: $(scmsettings.coalescence), Turbulent Droplet Diffusion: $(scmsettings.turbulent_droplet_diffusion_on)")

bins = binning_func(droplets_snapshots[time_idx], 1.0, runsettings, coagsettings, indices=x1:x2)
mids = 0.5*(radius_bins_edges[1:end-1] + radius_bins_edges[2:end])*1e6
plot(mids, max.(bins,1e-15), label="Binned DSD", yscale=:log10, xlims=(1,80))

heatmap(grid.output.ql)