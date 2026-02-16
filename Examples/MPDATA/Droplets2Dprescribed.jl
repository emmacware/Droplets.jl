using Droplets


grid = (20, 20)
dt = 1
dx = 1500/ grid[1]
dy = 1500 / grid[2]
omega = 0.6
xc = 0.5 * grid[1] * dx
yc = 0.5 * grid[2] * dy
# nsteps = 300
# n_corr = 2  
# output = 10
# n_output = nsteps/output
sf = make_stream_function(grid, dx, dy, omega, xc, yc)
u,v,GCx,GCy = nondivergent_vector_field_2d(dx,dy,grid, dt, sf)

# create fields to be advected:
qvarray = fill!(7.5, zeros(grid))
θ = fill!(289,  zeros(grid))
P0 = 1015 # hPa

function rhod_of_zZ(self, zZ)
    p =p_of_zconst_th_and_initial_wvmr.(P0, θ, self.initial_water_vapour_mixing_ratio, zZ)
    rhod = rho_d(p, self.initial_water_vapour_mixing_ratio, θ)
end

function p_of_zconst_th_and_initial_wvmr(p0, θ, water_vapour_mixing_ratio, z)
    z0 = 0
    Rq = constants.Rv / (1 / water_vapour_mixing_ratio + 1) + constants.Rd / (
        1 + water_vapour_mixing_ratio
    )
    Rd_over_c_pd = constants.Rd / constants.cp_d
    arg = ((p0 / 1000).^ (Rd_over_c_pd) - (z - z0) * Rd_over_c_pd * constants.g / θ / Rq)
    return 1000 * arg.^ (1 / Rd_over_c_pd)
end
# The initial air-density profile corresponds to the hydrostatic equilibrium
# with a pressure of 1015 hPa at the bottom of the domain.

# find initial saturation
# sat = saturation.(qvarray, temp, pressure)

# create a droplet field
droplets = init_spatial_uniformvol_dry_sd(dist,settings::coag_settings{FT},
    spatial::spatial_settings{FT})

for cell in numgrid
    i = findall(droplets.cell_id .== cell)
    T = Temp[cell]
    set_X_crit!.(droplets,i,kappa,T)
end

# create velocity field
u,v,GCx,GCy = nondivergent_vector_field_2d()

tmp = mpdata_tmp(ϕ,GCx,GCy)
run_settings = mpdata_settings(grid, n_corr=n_corr,
    horizontal_boundary_condition=Periodic(),
    vertical_boundary_condition=Periodic())
