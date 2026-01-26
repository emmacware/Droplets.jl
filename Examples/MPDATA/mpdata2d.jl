using Plots
using Droplets

dt = 0.1
dx = 1.0
dy = 1.0
omega = 0.1
h = 4.0
h0 = 1

grid = (100, 100)
r = 15.0 * dx
x0 = 50 * dx
y0 = 75 * dy
xc = 0.5 * grid[1] * dx
yc = 0.5 * grid[2] * dy
nsteps = 300
n_corr = 2  
output = 10
n_output = nsteps/output

ϕ = zeros(grid[1], grid[2])
for i in 1:grid[1], j in 1:grid[2]
    x,y = (i * dx, j * dy)
    tmp = (x - x0)^ 2 + (y - y0) ^ 2
    if tmp <= r^2
        pdf = h0 + h - sqrt(tmp / (r / h)^2)
    else
        pdf = 0
    end
    ϕ[i,j] = pdf 
end

sf = make_stream_function(grid, dx, dy, omega, xc, yc)

u,v,GCx,GCy = nondivergent_vector_field_2d(dx,dy,grid, dt, sf)
tmp = mpdata_tmp(ϕ,GCx,GCy)
run_settings = mpdata_settings(grid, n_corr=n_corr,
    horizontal_boundary_condition=Periodic(),
    vertical_boundary_condition=Periodic())

#################################################

anim = @animate for n in 1:n_output
    for _ in 1:output
        mpdata_step!(ϕ,GCx,GCy,tmp, run_settings)
        
        surface(ϕ', c=:viridis, clims=(0,1),
            xlabel="x", ylabel="y", title="Julia MPdata!!", aspect_ratio=1,zlims=(0, 4.2))
    end
end

gif(anim, "mpdata_periodic.gif", fps=10)
