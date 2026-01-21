using Droplets
FT = Float64


# From Wyant et al., 2007
inv_height = 875.0  # m, inversion height, or where qt = 8 g/kg (should this be changing?)
ρ_inv = 1.12e-3 # kg/m^3, air density @ inversion height
N0 = 55e6      #cm-3 -> m-3 initial number concentration
D = 3.75e-6 # s^-1, constant Divergence
u_star = 0.25 # m/s, friction velocity
prescribed_u(z) = 3 + 4*(z/1000) # m/s
prescribed_v(z) = -9 + 5.6*z/1000        # m/s
prescribed_w(z) = -D*z # m/s
θl(z) = z <= inv_height ? 288.0 : 295.0 + (z - inv_height)^(1/3)
qt(z) = z <= inv_height ? 9.45e-3 : 5e-3 - 3e-3(1-exp(-(z - inv_height)/500)) # q total
surface_latent_heat_flux = 93 # W/m^2
surface_sensible_heat_flux = 16 # W/m^2
surface_zonal_momentum_flux = prescribed_u(0)*u_star^2 / sqrt(prescribed_u(0)^2 + prescribed_v(0)^2)
surface_meridional_momentum_flux = prescribed_v(0)*u_star^2 / sqrt(prescribed_u(0)^2 + prescribed_v(0)^2)

#######################################

grid = create_scm_grids(100, 20.0, 1)
#fill in initial profiles
for k in 1:length(grid.centers_z)
    z = grid.centers_z[k]
    grid.states.θ[k] = θl(z)
    grid.states.qv[k] = qt(z)
    grid.states.P[k] = 101325.0 * (1 - 2.25577e-5 * z)^5.25588
    # grid.wind.u[k1] = prescribed_u(z)
    # grid.wind.v[k+1] = prescribed_v(z)
    # grid.wind.w[k] = prescribed_w(z)
end


# p1 = plot(grid.states.θ, xlabel="Height (m)", ylabel="Potential Temperature (K)", title="Initial Potential Temperature Profile", legend=false)
# p2 = plot(grid.states.qv*1000, xlabel="Height (m)", ylabel="Water Vapor Mixing Ratio (g/kg)", title="Initial Water Vapor Mixing Ratio Profile", legend=false)
# p3 = plot(grid.states.P/100, xlabel="Height (m)", ylabel="Pressure (hPa)", title="Initial Pressure Profile", legend=false)
# p4 = plot(grid.centers_z, grid.wind.u[2:end-1], xlabel="Height (m)", ylabel="Zonal Wind (m/s)", title="Initial Zonal Wind Profile", legend=false)
# p5 = plot(grid.centers_z, grid.wind.v[2:end-1], xlabel="Height (m)", ylabel="Meridional Wind (m/s)", title="Initial Meridional Wind Profile", legend=false)
# p6 = plot(grid.centers_z[1:end-1], grid.wind.w[2:end-1], xlabel="Height (m)", ylabel="Vertical Wind (m/s)", title="Initial Vertical Wind Profile", legend=false)
plot(p1, p2, p3, layout=(3,1), size=(800,900))
