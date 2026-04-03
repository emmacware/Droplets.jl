
using Droplets
using OrdinaryDiffEq
using Plots
using ComponentArrays
using OrdinaryDiffEqBDF
include("initial_state_functions.jl")
include("forward_solve.jl")
FT = Float64


# -----------------------------------------------------------------------
# Shipway and Hill (2012) kinematic 1D (KiD) test case
#
# Reference: Shipway & Hill, 2012, GMD 5, 1369–1384
# -----------------------------------------------------------------------

# --- Numerical settings -------------------------------------------------
Z_max        = 1500.0   # m, domain top
dz           = 25.0     # m, grid spacing (60 levels)
nz           = Int(Z_max / dz)
dt           = 1.0      # s
Ns_per_grid  = 16       # superdroplets per grid cell


# --- Kinematic velocity field -------------------------------------------
# Prescribed sinusoidal updraft; evaluated on face levels
w_max        = 2.0      # m/s, S&H 2012 Table 1
t1           = 600.0    # s, time of max updraft
prescribed_w(t) = w_max * sin(π * t / t1)


# --- Initial thermodynamic profiles (S&H 2012) -------------------------
# Dry adiabat: θ uniform throughout column
θ_0          = 283.15   # K, surface potential temperature
H_qv         = 2200.0   # m, q_v scale height (sets sub-cloud RH ≈ 80 % at surface)
q_v0         = 7.5e-3   # kg/kg, surface specific humidity
P_surface    = 101325.0 # Pa

θl(z)   = θ_0
qt(z)   = q_v0 * exp(-z / H_qv)


# --- Aerosol distribution (S&H 2012 Table 1) ---------------------------
# Ammonium sulfate, single-mode lognormal
kappa_ammonium_sulfate  = 0.61
ammonium_sulfate_density = 1.78e3  # kg/m^3
N_a  = 100 * ccm_to_cm            # m^-3, 100 cm^-3
r_dry_mode = 0.04e-6               # m, median dry radius
σ_g  = 1.4                         # geometric standard deviation
n0   = N_a

aerosol_dist = LogNormal(log(r_dry_mode), log(σ_g))


# --- Settings structs ---------------------------------------------------
spatialsettings      = spatial_settings_1d{FT}(Nz=nz, Z_max=Z_max)
coagsettings         = coag_settings{FT}(Ns=Ns_per_grid*nz, ΔV=dz*spatialsettings.area_per_grid, n0=n0, Δt=dt)
condensationsettings = condensation_settings{FT}(kappa=kappa_ammonium_sulfate, ρ_solute=ammonium_sulfate_density, Δt=dt)
mpdatasettings       = mpdata_settings_1d(nz, vertical_boundary_condition=NoFlux())
# Kinematic model: no surface fluxes, no TKE
scmsettings          = scm_settings{FT}(Δt=dt, surface_latent_heat_flux=0.0, surface_sensible_heat_flux=0.0)
tkesettings          = tke_settings{FT}(u_star=0.0)
diagnosticsettings   = diagnostic_settings()


# --- Initialise Eulerian grid and superdroplets -------------------------
grid     = initialize_kid_environment(nz, dz, P_surface, θl, qt, prescribed_w)
droplets = init_droplets_kid1d(aerosol_dist, coagsettings, spatialsettings)


# --- Initialise droplets to Köhler equilibrium -------------------------
for k in 1:nz
    drop_idx = findall(i -> droplets.cell_id[i] == k, 1:length(droplets.X))
    T_k  = grid.states.T[k]
    qv_k = grid.states.qv[k]
    P_k  = grid.states.P[k]
    S_env = min(sat(qv_k, P_k) / esat(T_k), 0.95)
    find_equilibrium_radius.(droplets, drop_idx, kappa_ammonium_sulfate, T_k, S_env)
end


# --- Pre-run diagnostics ------------------------------------------------
coagdata             = coagulation_run_spatial{FT}(nz, coagsettings.Ns, droplets)
condensation_integrator = create_condensation_integrator(grid, droplets, condensationsettings, coagsettings, spatialsettings, constants)
mpdata_tmp           = mpdata_tmp_1d(grid.states.qv, grid.faces_z)

sd_fill_diagnostics(droplets, grid, spatialsettings, diagnosticsettings)
plot_env_profiles(grid)


# --- Time loop (3600 s = 1 hour) ----------------------------------------
n_steps = 3600
for i in 1:n_steps
    if i % 600 == 0
        println("Timestep: $i  /  $n_steps  ($(i*dt/60) min)")
    end
    single_column_timestep(grid, dt, droplets, coagsettings, spatialsettings,
        condensationsettings, condensation_integrator, coagdata, diagnosticsettings,
        prescribed_w, mpdata_tmp, mpdatasettings, constants, scmsettings, tkesettings)
end

plot_env_profiles(grid)
