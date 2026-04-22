using Droplets
using Printf

FT = Float64
nz  = 30
dz  = FT(50.0)
dt  = FT(5.0)

drops = droplet_attributes_1d{FT}(
    Int[], FT[], FT[], FT[], Int[],
    Vector{Union{UnitRange{Int},Nothing}}([nothing for _ in 1:nz]),
    Int[]
)
spatial = spatial_settings_1d{FT}(Nz=nz, Z_max=FT(nz*dz))
grid    = create_scm_grids(nz, dz, drops; spatial=spatial)

inv_height = 795.0
P_ref = FT(101780.0)
for k in 1:nz
    z = (k - 0.5) * dz
    grid.states.P[k]  = P_ref
    grid.states.θ[k]  = z <= inv_height ? FT(288.3) : FT(295.0) + (z - inv_height)^(1/3)
    grid.states.qv[k] = z <= inv_height ? FT(9.45e-3) : FT(5e-3)
    grid.states.e[k]  = FT(1.0)
    grid.wind.u[k]    = 3.0 + 4.0*(z/1000)
    grid.wind.v[k]    = -9.0 + 5.6*z/1000
end
grid.states.ql_tmp .= zeros(FT, nz)
grid.states.ρ .= ρ_calc_θ.(P_ref, grid.states.θ, grid.states.qv, constants)

tke = tke_settings{FT}(u_star=FT(0.25),
    geostrophic_u = z -> 3.0 + 4.0*(z/1000),
    geostrophic_v = z -> -9.0 + 5.6*z/1000)

delta = dz
l  = zeros(FT, nz)
mixing_length_deardorff!(l, grid, tke.l_inf, constants)

println("\n=== TKE BUDGET PER CELL (before timestep) ===")
println(@sprintf("%-4s %-7s %-7s %-8s %-8s %-10s %-10s %-10s %-8s %-8s",
    "k","z(m)","l(m)","K_m","K_h","P_shear","P_buoy","epsilon","net","τ_decay(s)"))

for k in 1:nz
    z      = (k - 0.5) * dz
    e_k    = grid.states.e[k]
    sqrte  = sqrt(max(e_k, tke.e_min))
    K_m_k  = tke.c_m * l[k] * sqrte
    K_h_k  = K_m_k * (1 + 2*l[k]/delta)
    N2     = calculate_buoyancy_frequency(grid, k, constants)
    S2     = S2_flow_deformation(grid, k, constants)
    c_eps  = 0.19 + 0.51 * (l[k] / delta)
    P_shear = K_m_k * S2
    P_buoy  = -K_h_k * N2
    epsilon = c_eps * e_k^(3/2) / l[k]
    net     = P_shear + P_buoy - epsilon
    tau     = e_k / max(epsilon, 1e-30)   # decay timescale if no production
    println(@sprintf("%-4d %-7.1f %-7.4f %-8.5f %-8.5f %-10.3e %-10.3e %-10.3e %-8.3e %-8.1f",
        k, z, l[k], K_m_k, K_h_k, P_shear, P_buoy, epsilon, net, tau))
end

# Run 10 timesteps and track e
println("\n=== e after each timestep (min / mean / max) ===")
for i in 1:10
    no_prog_tke_timestep!(grid, tke, constants, dt)
    println(@sprintf("step %2d:  min=%.3e  mean=%.3e  max=%.3e  e[1]=%.3e",
        i, minimum(grid.states.e), sum(grid.states.e)/nz,
        maximum(grid.states.e), grid.states.e[1]))
end
