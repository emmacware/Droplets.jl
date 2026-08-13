using Revise
using Droplets
# using DifferentialEquations

using OrdinaryDiffEq

using Plots
using StatsPlots
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
const dz = 10.0 #m
const nz = Int(Z_max/dz)
const dt = 1.0 #s
const dt_coag = 0.1
const dt_cond = 0.1
const Ns_per_grid = 200
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
settling         = DynOFF(),
spinupsaturation = DynON(),
coalescence       = DynOFF(),
# turbulent_droplet_diffusion_on = DynOFF(),


)
scmsettings = scm_settings{FT}(; base_scm...,
REM              = DynOFF(),
settling         = DynON(),
spinupsaturation = DynOFF(),
coalescence       = DynON(),
turbulent_droplet_diffusion_on = DynON(),


)
scmradsettings = scm_settings{FT}(; base_scm...,
REM              = DynON(),
settling         = DynON(),
spinupsaturation = DynOFF(),
coalescence       = DynON(),
turbulent_droplet_diffusion_on = DynON(),

)
diagnosticsettings = diagnostic_settings()



# Create environmnent
grid, droplets, coagdata,conddata,raddata,mpdatatmp,turbdata = initialize_scm_environment(
    nz, dz, P_surface, θl_initial, qt_initial, prescribed_w, dist,
    coagsettings,spatialsettings,condensationsettings,tkesettings,constants
    )

radgrid, droplets, coagdata,conddata,raddata,mpdatatmp,turbdata = initialize_scm_environment(
    nz, dz, P_surface, θl_initial, qt_initial, prescribed_w, dist,
    coagsettings,spatialsettings,condensationsettings,tkesettings,constants
    )

# inv_idx = findfirst(k -> grid.centers_z[k] >= inv_height, 1:nz)
# grid.states.e[1:inv_idx] .= 0.01#1e-6

CArad = raddata.CArad
R_array = range(CArad.lookup_lw_cld.bounds[1], CArad.lookup_lw_cld.bounds[2],
                length=CArad.lookup_lw_cld.dims[3])  # μm

const absliq_r_interp = ntuple(CArad.nband_lw) do ibnd
    ext, ssa, _ = LookUpTables.getview_liqdata(CArad.lookup_lw_cld, ibnd)
    absliq = collect(ext .* (1 .- ssa))
    linear_interpolation(R_array, absliq, extrapolation_bc=Flat())
end


spinup_step = round(Int, scmsettings.spinup_time / dt)
# bins_output = Dict()
# ensemble_output = Dict()
# To skip re-running, load a previous save instead:
# ensemble_output, bins_output = load("ensemble_output.jld2", "ensemble_output", "bins_output")
num_bins = 50
radius_bins_edges = 10 .^ range(log10(1*1e-8), log10(1e2*1e-6), length=num_bins+1) 
# radius_bins_edges = range(0.5*1e-6,100e-6, length=num_bins+1) 

runsettings = run_settings{FT}(num_bins=num_bins,radius_bins_edges=radius_bins_edges,normalize_bins_dlnr=false,binning_method=mass_density_lnr)
seeds = 10


for num_seeds in 1:seeds
    # for n0_change in [n0,1.4e8,n0/2,2.7e8, n0*2]
    for n0_change in [n0,n0/2, n0*2]
        for rad in [false, true]
            coagsettings = coag_settings{FT}(Ns=Ns_per_grid*nz,ΔV=dz*spatialsettings.area_per_grid, n0=n0_change,Δt=dt_coag,kernel=hydrodynamic,hydrodynamic_collision_eff_func=true)
            scmrunsettings = rad ? scmradsettings : scmsettings
            
            

            if (rad, n0_change, num_seeds) in keys(ensemble_output)
                continue
            end
            println("Running simulation with rad=", rad, " n0=", n0_change, " seed=", num_seeds)

            droplets_snapshots = Vector{Any}(undef, Int(spatialsettings.t_max / 150))

            Random.seed!(seed + num_seeds)

            grid, droplets, coagdata,conddata,raddata,mpdatatmp,turbdata = initialize_scm_environment(
            nz, dz, P_surface, θl_initial, qt_initial, prescribed_w, dist,
            coagsettings,spatialsettings,condensationsettings,tkesettings,constants
            )
            inv_idx = findfirst(k -> grid.centers_z[k] >= inv_height, 1:nz)
            grid.states.e[1:inv_idx] .= 0.01#1e-6
            ensemble_output[rad,n0_change,num_seeds] = grid.output
            bins_output[rad,n0_change,num_seeds] = droplets_snapshots
            try
                for i in 1:spinup_step
                    if i*dt % 3600 == 0
                        println("Timestep: ", i*dt)
                    end
                    if i*dt % 150 == 0
                        droplets_snapshots[Int(div(i*dt,150))] = deepcopy(droplets)
                    end

                    single_column_timestep(grid,dt,droplets,coagsettings,spatialsettings,condensationsettings,
                    coagdata,conddata,raddata,turbdata,
                    diagnosticsettings,prescribed_w, mpdatatmp, mpdatasettings,constants,scmspinupsettings,tkesettings,absliq_r_interp,i)

                end

                let calls = sum(conddata.bisection_calls), hits = sum(conddata.bisection_maxdepth_hits),
                    growth = sum(conddata.bisection_maxdepth_hits_growth), shrink = sum(conddata.bisection_maxdepth_hits_shrink)
                    pct = calls > 0 ? round(100 * hits / calls, digits=4) : 0.0
                    println("Spin-up bisection max-depth fallback: ", hits, " / ", calls, " (", pct, "%)  [growth: ", growth, ", shrink: ", shrink, "]")
                end
                conddata.bisection_calls .= 0
                conddata.bisection_maxdepth_hits .= 0
                conddata.bisection_maxdepth_hits_growth .= 0
                conddata.bisection_maxdepth_hits_shrink .= 0

                for i in (spinup_step+1):Int(spatialsettings.t_max / dt)
                    if i*dt % 3600 == 0
                        println("Timestep: ", i*dt)
                    end
                    if i*dt % 150 == 0
                        droplets_snapshots[Int(div(i*dt,150))] = deepcopy(droplets)
                    end

                    single_column_timestep(grid,dt,droplets,coagsettings,spatialsettings,condensationsettings,
                    coagdata,conddata,raddata,turbdata,
                    diagnosticsettings,prescribed_w, mpdatatmp, mpdatasettings,constants,scmrunsettings,tkesettings,absliq_r_interp,i)
                end

                let calls = sum(conddata.bisection_calls), hits = sum(conddata.bisection_maxdepth_hits),
                    growth = sum(conddata.bisection_maxdepth_hits_growth), shrink = sum(conddata.bisection_maxdepth_hits_shrink)
                    pct = calls > 0 ? round(100 * hits / calls, digits=4) : 0.0
                    println("Main-run bisection max-depth fallback: ", hits, " / ", calls, " (", pct, "%)  [growth: ", growth, ", shrink: ", shrink, "]")
                end
            catch e
                @warn "Simulation failed for rad=$(rad), n0=$(n0_change), seed=$(num_seeds) with error: $(e)"
                # ensemble_output[rad,n0_change,num_seeds] = nothing
                # bins_output[rad,n0_change,num_seeds] = nothing
                continue
            end


            ensemble_output[rad,n0_change,num_seeds] = grid.output
            bins_output[rad,n0_change,num_seeds] = droplets_snapshots
            jldsave("ensemble_output.jld2"; ensemble_output, bins_output)


        end
    end
end

for field in fieldnames(typeof(ensemble_output[false,n0,1]))
    getfield(grid.output, field) .= getfield(ensemble_output[false,n0,1], field)
end

#change 0 to NaN for plotting
#find the column where all values are 0, and replace with NaN
for field in fieldnames(typeof(ensemble_output[false,n0,1]))
    data = getfield(grid.output, field)
    for t in 1:size(data,2)
        if all(data[:,t] .== 0.0)
            data[:,t] .= NaN
        end
    end
end

penv = plot_env_profiles(grid)
# #put tkesettings.turbulence_scheme in the title
ptime = plot!(plot_output_timeseries(grid),plot_title="REM: $(scmsettings.REM)")#, Coalescence: $(scmsettings.coalescence), Turbulent Droplet Diffusion: $(scmsettings.turbulent_droplet_diffusion_on)")

num = ensemble_output[(false, n0,1)].number * 1e-6/50
nummasked = ifelse.(num .< 20, NaN, num)
heatmap(nummasked,clims=(0,100))

plot()
plot_ensemble_field(:number, ensemble_output, grid; key_filter=(false, n0),   scale=1e-6/50, axis=:height, show_ribbon=false, t_window=(0,1), label="", color=:darkturquoise,  linestyle=:dot,   linewidth=2)


tCpR = mapslices(x -> begin v = filter(!isnan, x); isempty(v) ? NaN : mean(v) end, nummasked, dims=1)
#time series mean of only values above 20, N_cloud (threshold applied after run-averaging)
plot(tCpR', label="No REM", color=:aquamarine4, linewidth=2)

#time series of integrated tke (ensemble_output[(false, n0,1)].e) over Height
tke_integrated = mapslices(x -> begin v = filter(!isnan, x); isempty(v) ? NaN : sum(v) end, ensemble_output[(false, n0,1)].e*dz, dims=1)
plot(tke_integrated', label="No REM", color=:aquamarine4, linewidth=2, xlabel="Time (hr)", ylabel="Integrated TKE (m²/s²)", title="Integrated TKE over Height")