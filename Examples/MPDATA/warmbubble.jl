using Droplets, Plots, LinearAlgebra
include("warmbubble_funcs.jl")


########################
# DOMAIN + CONSTANTS
########################
nx, ny = 110,100
x_domain, y_domain = 2000,2000
dx, dy = x_domain/nx, y_domain/ny
dt = 1.5
nsteps = 600#Int(round(600/dt))

xc, yc = 1000,260
rc = 250.0
Δθamp = 0.5

# physical constants
R_d = 287.0
cp  = 1004.0
cv  = 717.0
g   = 9.81
p0  = 100000.0
θ_ref = 300.0

########################
# Scalars at cell centers
########################
θ = fill(θ_ref, nx, ny)
ρ = zeros(nx, ny)
T = zeros(nx, ny)
exner = ones(nx, ny)

# initialize warm bubble
for i in 1:nx, j in 1:ny
    x = (i-0.5)*dx
    y = (j-0.5)*dy
    r = sqrt((x-xc)^2 + (y-yc)^2)
    # Δθ = r ≤ rc ? Δθamp * (1 - r/rc) : 0.0
    Δθ = r ≤ rc ? Δθamp : 0.0

    θ[i,j] = θ_ref + Δθ
    exner[i,j] = 1 - g/(cp*θ[i,j])*y
    T[i,j] = θ[i,j]*exner[i,j]
    ρ[i,j] = p0/(R_d*θ[i,j]) * exner[i,j]^(cv/R_d)
end

########################
# Velocities on staggered grid
########################
u = zeros(nx+1, ny)    # x-faces
v = zeros(nx, ny+1)    # y-faces

# MPDATA tmp arrays (full size including halos)
tmpθ = mpdata_tmp(θ, u, v)
u_centers = 0.5 .* (u[1:end-1, :] + u[2:end, :])
v_centers = 0.5 .* (v[:, 1:end-1] + v[:, 2:end])
tmp_u = mpdata_tmp(u_centers, u, v)
tmp_v = mpdata_tmp(v_centers, u, v)

# MPDATA run settings
runsettings = mpdata_settings((nx, ny), n_corr=2,
    horizontal_boundary_condition=Periodic(),
    vertical_boundary_condition=Periodic(),
    nonoscillatory=true,
    infinite_gauge=true)
    

ARRAYS = pressure_loop_tmp_struct()
    # interpolate staggered velocities to cell centers

avg_u = 0.5 .* (u[1:end-1, :] .+ u[2:end, :])
avg_v = 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])
GCx = zeros(nx+1, ny)
GCy = zeros(nx, ny+1)

GCx_centers = 0.5 .* (GCx[1:end-1, :] .+ GCx[2:end, :])
GCy_centers = 0.5 .* (GCy[:, 1:end-1] .+ GCy[:, 2:end])

# ARRAYS.Phi[1:end-1, 1:end-1] .= exner

########################
# TIME LOOP
########################

anim = @animate for n in 1:nsteps

    println(n)


    calc_gc_extrapolate_in_time(ARRAYS, avg_u, avg_v) 
    calc_gc_interpolate_in_space(ARRAYS, GCx, GCy, dx,dy,dt) 
    apply_half_rhs(ARRAYS, avg_u, avg_v,dt)
    fill_stash(ARRAYS, avg_u, avg_v)
    vip_rhs_apply(ARRAYS, avg_u, avg_v,dt,nx,ny)

    mpdata_step!(θ, GCx, GCy, tmpθ, runsettings,f=0.5)
    mpdata_step!(avg_u, GCx,GCy, tmp_u, runsettings,f=0.5)
    mpdata_step!(avg_v, GCx, GCy, tmp_v, runsettings,f=0.5)


    update_rhs(ARRAYS, θ,g,θ_ref)

    apply_half_rhs(ARRAYS, avg_u, avg_v,dt)

    finilize_rhs(ARRAYS, avg_u, avg_v,nx,ny)
    pressure_solver_update(ARRAYS,avg_u,avg_v)
    pressure_solver_apply(ARRAYS,avg_u,avg_v,dt,nx,ny)
    if n%5 == 0
        heatmap(θ'.-θ_ref, levels=5,c=:viridis, clims=(299.95-θ_ref, 300.55-θ_ref),
            xlabel="x", ylabel="y", title="Rising Warm Bubble",
            aspect_ratio=1, #zlims=(θ_ref, θ_ref+Δθamp)
            )
        # contour!(θ'.-θ_ref, levels=5,#c=:viridis, clims=(299.95-θ_ref, 300.55-θ_ref),
        #     xlabel="x", ylabel="y", title="Rising Warm Bubble",
        #     aspect_ratio=1, #zlims=(θ_ref, θ_ref+Δθamp)
        #     )
        plot!(grid=false)
    end
end

gif(anim, "cgrid_rising_bubble_corrected.gif", fps=20)








