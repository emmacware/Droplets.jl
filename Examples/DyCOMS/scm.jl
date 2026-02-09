using Droplets
using DifferentialEquations
using Plots
using ComponentArrays
include("hydrostatic_ode.jl")
include("init_dycoms_drops.jl")
include("condensation.jl")
include("scm_motion.jl")
include("forward_solve.jl")
FT = Float64


# From Wyant et al., 2007
inv_height = 795.0  # m, inversion height, or where qt = 8 g/kg (should this be changing?)

# The inversion height is found
# by interpolating between levels to the height where qt drops
# to 8 g kg 1. The cloud base is determined by interpolating
# to find the lowest level where cloud fraction reaches 0.5.

ρ_inv = 1.12#e-3 # kg/m^3, air density @ inversion height
N0 = 55e6      #cm-3 -> m-3 initial number concentration
D = 3.75e-6 # s^-1, constant Divergence
u_star = 0.25 # m/s, friction velocity
prescribed_u(z) = 3 + 4*(z/1000) # m/s
prescribed_v(z) = -9 + 5.6*z/1000        # m/s
prescribed_w(z) = -D*z # m/s
θl(z) = z <= inv_height ? 288.3 : 295.0 + (z - inv_height)^(1/3) #boundary layer liquid water potential temperature
qt(z) = z <= inv_height ? 9.45e-3 : 5e-3 - 3e-3(1-exp(-(z - inv_height)/500)) # q total
surface_latent_heat_flux = 93 # W/m^2
surface_sensible_heat_flux = 16 # W/m^2
surface_zonal_momentum_flux = prescribed_u(0)*u_star^2 / sqrt(prescribed_u(0)^2 + prescribed_v(0)^2)
surface_meridional_momentum_flux = prescribed_v(0)*u_star^2 / sqrt(prescribed_u(0)^2 + prescribed_v(0)^2)


halo=1
max_height = 1500.0
dz = 2.5
nz = Int(max_height/dz)
dt = 1.0
grid = create_scm_grids(nz, dz, halo)
ρ_init = zeros(length(grid.centers_z))
#pref =1000 mb, P_surface = 101780 Pa from Ackerman et al 2009
P_surface = 101780.0

constants_dycoms_rf02 = Constants{Float32}(P0=100000.0)

press = ODEProblem(dP_dz,P_surface,(0,maximum(grid.centers_z)+1),constants_dycoms_rf02)
P_init = solve(press,Tsit5(), reltol = 1e-10, abstol = 1e-10,saveat = grid.centers_z)
grid.states.P[2:nz+1] .= P_init.u


#######################################


#fill in initial profiles
for z_idx in 1:length(grid.centers_z)

    k = z_idx + halo
    z = grid.centers_z[z_idx]
    grid.states.θ[k] = θl(z)
    grid.states.qv[k] = qt(z)
    
    T = in_situ_temperature(θl(z),qt(z),grid.states.P[k],constants_dycoms_rf02)
    T_virtual = T * (1 + 0.61 * grid.states.qv[k])
    ρ = grid.states.P[k] / (constants_dycoms_rf02.Rd * T_virtual)
    grid.states.ρ[k] = ρ
    grid.states.T[k] = T
end

for z_idx in 1:length(grid.faces_z)
    k = z_idx + halo
    z = grid.faces_z[z_idx]
    grid.wind.u[k] = prescribed_u(z)
    grid.wind.v[k] = prescribed_v(z)
    grid.wind.w[k] = prescribed_w(z)
end

# interior_idx = CartesianIndices(halo+1:nz+halo)
interior_idx = halo+1:nz+halo


p1 = plot(grid.states.θ[interior_idx],grid.centers_z, ylabel="Height (m)", xlabel="Potential Temperature liquid", title="θl initial", legend=false)
p2 = plot(grid.states.qv[interior_idx]*1000,grid.centers_z, ylabel="Height (m)", xlabel="q_tot (g/kg)", title="q_tot initial", legend=false)
p3 = plot(grid.states.P[interior_idx],grid.centers_z, ylabel="Height (m)", xlabel="Pressure (Pa)", title="P", legend=false)
p4 = plot(grid.states.ρ[interior_idx],grid.centers_z, ylabel="Height (m)", xlabel="density (kg/m3)", title="ρ initial", label = "ρ")
p4 = vline!([ρ_inv], line=:dash, color=:red, label="ρ_inv")
p5 = plot(grid.states.T[interior_idx],grid.centers_z, ylabel="Height (m)", xlabel="Temperature(K)", title="T initial", legend=false)

plot(p1, p2, p3,p4,p5, layout=(3,2), size=(800,900))



# κ = 0.61 #ammonium sulfate, Petters and Kreidenweis (2007)

#θl(p,T) = (pref/p)^(Rd/Cp_air)*(T-L*ql/cp)
# surface air density of 1.21 is implicit

#from Ackerman et al., 2009: 
# with mo-lecular weight 115 g mol21, dry density 1.78 g cm23, and
# two ions dissolved per molecule). The total number,
# mode radius, and geometric standard deviation for the
# two modes are 125 and 65 cm23, 0.011 and 0.06 mm, and
# 1.2 and 1.7, respectively

Ns = 20*nz
spatial = spatial_settings_1d{Float64}(Nz=nz, Z_max=max_height)
n0 = (65+ 125) * 1e6 # total number concentration in m^-3, converted from cm^-3
settings = coag_settings{Float64}(Ns=Ns, ΔV=spatial.z_grid_height*spatial.area_per_grid, n0=n0)
drops = init_droplets_dycoms_scm(settings, spatial)
diagnosticsettings = diagnostic_settings{Float64}(aerosol_cloud_cuttoff=1e-6, cloud_rain_cuttoff=40e-6)
coag_data = coagulation_run_spatial{FT}(nz, Ns)

#create other needed data structures
lnR = log.(volume_to_radius.(drops.X))
Y = ComponentVector{FT}(lnR = lnR, qvap=grid.states.qv, T = grid.states.T)
p = (nz,drops, grid.states, constants)
condensation_prob = ODEProblem(condensation_rhs, Y, (0.0, 1.0), p)
condensation_integrator = init(condensation_prob, Tsit5(),dt=0.5, reltol = 1e-10, abstol = 1e-10)

condensationsettings = condensation_settings{FT}()
#advection 
tmp = mpdata_tmp_1d(grid.states.qv, grid.wind.w)
mpdata_settings = mpdata_settings_1d(nz, n_corr=2, vertical_boundary_condition=NoFlux())
scmsettings = scm_settings{Float64}(init_random_seed=30, coag_threading=Serial(), scheme=none(), Δt=1.0)
diagnosticsettings = diagnostic_settings{Float64}(aerosol_cloud_cuttoff=1e-6, cloud_rain_cuttoff=40e-6)


single_column_timestep(grid,dt,drops,settings,spatial,condensationsettings,condensation_integrator,
    coag_data,diagnosticsettings,prescribed_w, tmp, mpdata_settings,constants,scmsettings)