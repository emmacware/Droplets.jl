using Revise
using Droplets
using OrdinaryDiffEq
using Plots
using Distributions
using Random

include("../DyCOMS/radiation_call.jl")
include("../DyCOMS/initial_state_functions.jl")
include("../DyCOMS/forward_solve.jl")
# using OrdinaryDiffEqBDF

const FT = Float64

# ── Numerical settings ────────────────────────────────────────────────────────
const Z_max = 3000.0   # m
const dz    = 25.0     # m
const nz    = Int(Z_max / dz)   # 120
const dt    = 1.0      # s
const dt_coag = 0.1
const dt_cond = 0.1
const Ns_per_grid = 32
seed = 42
Random.seed!(seed)
const t_max    = 3600  # s  (60 min)
const t_output = 60.0  # s

# ── Shipway & Hill 2012 initial profiles ─────────────────────────────────────
θ_initial(z) = z <= 740.0 ? 297.9 :
               297.9 + (312.66 - 297.9) * (z - 740.0) / (3260.0 - 740.0)

qv_initial(z) = z <= 740.0 ? 0.015 + (0.0138 - 0.015) * z / 740.0 :
                              0.0138 + (0.0024 - 0.0138) * (z - 740.0) / (3260.0 - 740.0)
# θ_initial(z) = calc_θ_dry(constants, θs_initial(z), qv_initial(z))
                    
# ── KiD kinematic updraft: ρw = ρw₁ sin(πt/t₁) for t < t₁, else 0 ──────────
const rho_times_w_1 = 2.0    # kg m⁻² s⁻¹
# const rho_STP_kid   = 1.225  # kg m⁻³  (used for reservoir depth estimate)
# const apprx_w1      = rho_times_w_1 / rho_STP_kid   # ≈ 1.63 m/s
const t_w           = 600.0  # s
prescribed_rho_w(t) = t < t_w ? rho_times_w_1 * sin(π * t / t_w) : 0.0

# ── Lognormal aerosol: mode = 0.04 μm, σ_geom = 1.4, N = 50 cm⁻³ ───────────
const n0         = 50.0 * ccm_to_cm   # m⁻³
const r_mode_kid = 0.04e-6            # m
const sigma_geom = 1.4
const sigma_ln   = log(sigma_geom)
const mu_ln      = log(r_mode_kid) + sigma_ln^2   # ensures mode = r_mode_kid
initial_aerosol_dist = LogNormal(mu_ln, sigma_ln)

const kappa_kid               = 0.9
const ammonium_sulfate_density = 1.78e3   # kg m⁻³
const P_surface_kid            = 100000.0 # Pa  (Shipway & Hill 2012 reference)

# ── Settings structs ──────────────────────────────────────────────────────────



#Settings Structs
spatialsettings = spatial_settings_1d{FT}(Nz=nz, Z_max=Z_max,dt=dt, t_max=t_max, dt_output=t_output,area_per_grid=10.0)#5.0)

coagsettings = coag_settings{FT}(Ns=Ns_per_grid*nz,ΔV=dz*spatialsettings.area_per_grid, n0=n0,Δt=dt_coag,kernel=hydrodynamic,hydrodynamic_collision_eff_func=false)

condensationsettings = condensation_settings{FT}(kappa=kappa_kid,ρ_solute = ammonium_sulfate_density,Δt=dt_cond)

mpdatasettings = mpdata_settings_1d(nz,nonoscillatory=true, vertical_boundary_condition=Extrapolated(),infinite_gauge=true)

tkesettings = tke_settings{FT}()#u_star=u_star, geostrophic_u=geostrophic_u, geostrophic_v=geostrophic_v,
    #LHF=surface_latent_heat_flux, SHF=surface_sensible_heat_flux)
    
diagnosticsettings = diagnostic_settings()
   
#Default dynamics are all on
base_scm = (
    Δt                          = dt,
    turbulence                  = DynOFF(),
    motion                      = DynON(),
    advection                   = DynON(),
    radiation                  = DynOFF(),
    condensation                = DynON(),
    n_cond                      = round(Int, dt / dt_cond),
    n_coag                      = round(Int, dt / dt_coag),
    spinup_time                 = 0.0,
    turbulent_droplet_diffusion_on = DynOFF(),
)

# scmspinupsettings = scm_settings{FT}(; base_scm...,
# REM              = DynOFF(),
# settling         = DynOFF(),
# spinupsaturation = DynON(),
# coalescence       = DynOFF(),
# )
scmsettings = scm_settings{FT}(; base_scm...,
REM              = DynOFF(),
settling         = DynON(),
spinupsaturation = DynOFF(),
coalescence       = DynON(),
top_escape = DynON(),
thermo_feedback = DynOFF(),
)
diagnosticsettings = diagnostic_settings()


# Create environmnent
grid, droplets, coagdata,conddata,raddata,mpdatatmp,turbdata = initialize_scm_environment(
    nz, dz, P_surface_kid, θ_initial, qv_initial, z -> zero(FT), initial_aerosol_dist,
    coagsettings,spatialsettings,condensationsettings,tkesettings,constants
    )

absliq_r_interp = 0.0


for i in 1:Int(spatialsettings.t_max / dt)
    if i*dt % 100 == 0
        println("Timestep: ", i*dt)
    end

    rho_w_t = prescribed_rho_w(i * dt)
    w_const = rho_w_t #/ rho_STP_kid   # constant w at all faces for conservative anelastic MPDATA
    ρ_faces = vcat(grid.states.ρ[1], 0.5 * (grid.states.ρ[1:end-1] .+ grid.states.ρ[2:end]), grid.states.ρ[end]) #
    ρ_faces[1] = 2 * grid.states.ρ[1] - grid.states.ρ[2]  # Extrapolate ρ to bottom face
    ρ_faces[end] = 2 * grid.states.ρ[end] - grid.states.ρ[end-1]  # Extrapolate ρ to top face
    for k in 1:nz+1
        grid.wind.w[k] = w_const / ρ_faces[k]
    end
    prescribed_w(z) = rho_w_t / grid.states.ρ[clamp(floor(Int, z / dz) + 1, 1, nz)]

    single_column_timestep(grid,dt,droplets,coagsettings,spatialsettings,condensationsettings,
    coagdata,conddata,raddata,turbdata,
    diagnosticsettings,prescribed_w, mpdatatmp, mpdatasettings,constants,scmsettings,tkesettings,absliq_r_interp,i)

end

# for i in (spinup_step+1):Int(spatialsettings.t_max / dt)
#     if i*dt % 100 == 0
#         println("Timestep: ", i)
#     end

#     single_column_timestep(grid,dt,droplets,coagsettings,spatialsettings,condensationsettings,
#     coagdata,conddata,raddata,turbdata,
#     diagnosticsettings,prescribed_w, mpdatatmp, mpdatasettings,constants,scmsettings,tkesettings,absliq_r_interp,i)
#     grid.wind.z .= prescribed_w.(grid.centers_z, i*dt)
# end

# ensemble_output[num_seeds] = grid.output

penv = plot_env_profiles(grid)
#put tkesettings.turbulence_scheme in the title
ptime = plot!(plot_output_timeseries(grid),tskips=20)#, Coalescence: $(scmsettings.coalescence), Turbulent Droplet Diffusion: $(scmsettings.turbulent_droplet_diffusion_on)")

# savefig("DyCOMS_SCMno_wprime.pdf")

p1 = heatmap(range(1, 3600,length(grid.output.ql[1,:])),grid.centers_z,     
    grid.output.cloud_LWC,color=:viridis,clims=(0,1.5e-3),colorbar=false,title="Cloud LWC")
xlabel!("t (s)")
yaxis!("z (m)")

p2 = heatmap(range(1, 3600,length(grid.output.ql[1,:])),grid.centers_z,     
    grid.output.rain_LWC,color=:viridis,clims=(0,1.5e-3),colorbar=false,title="Rain LWC")
xlabel!("t (s)")
yaxis!("z (m)")
    
cbar = heatmap([0],range(0, 1.5e-3, length=100),  reshape(range(0, 1.5e-3, length=100), 1, :), color=:viridis, ylabel="LWC (kg/m³)",colorbar=false ,xticks=false)

l = @layout [a b c{0.05w}]
plot(p1, p2, cbar, layout=l, size=(900,400),plot_title="K1D Droplets.jl             ",left_margin=3Plots.mm, bottom_margin=3Plots.mm,top_margin=3Plots.mm)

# savefig("DropletsKid1d_LWC_heatmaps.pdf")

# heatmap(range(1, 3600,length(grid.output.ql[1,:])),grid.centers_z,     
# grid.output.rain_LWC,color=:viridis,clims=(0,1.5e-3))
# # t_skips=20
# p6 = plot(sat.(grid.output.qv[:,1:t_skips:end],grid.output.P[:,1:t_skips:end]) ./ esat.(T_from_theta.(grid.output.θ[:,1:t_skips:end],grid.output.P[:,1:t_skips:end],constants)), grid.centers_z, ylabel="z (m)", xlabel="ql (g/kg)", title="ql", legend=false)


# plot(sat.(grid.states.qv,grid.states.P) ./ esat.(T_from_theta.(grid.states.θ,grid.states.P,constants)), grid.centers_z, ylabel="z (m)", xlabel="ql (g/kg)", title="ql", legend=false)
