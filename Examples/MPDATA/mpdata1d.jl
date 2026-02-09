using Plots
using Droplets

dt = 0.1
dz = 1.0
h = 4.0
h0 = 1

nz = 200
z_max = nz * dz
x0 = 40 * dz
a = 1.0
sigma = 10.0
nsteps = 100
n_corr = 2  
output = 5
n_output = nsteps/output
x = collect(1:nz) * dz
ϕ = a * exp.(-(x.-x0).^2 ./ 2 ./ sigma^2)

z = fill(0.5, nz+1)
GCz = z * dt / dz

tmp = mpdata_tmp_1d(ϕ,GCz)
run_settings = mpdata_settings_1d(nz, n_corr=n_corr,
    vertical_boundary_condition=NoFlux())

#################################################

anim = @animate for n in 1:n_output
    for _ in 1:output
        mpdata_step!(ϕ,GCz,tmp, run_settings)
        
        plot(x,ϕ)
    end
end

gif(anim, "mpdata_1d.gif", fps=10)
