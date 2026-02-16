using Droplets
using Distributions
using Plots


FT = Float64
#coagsettings
Ns = 2^17
Δt = FT(1.0)
ΔV = FT(1e6)
golovin_kernel_coeff = FT(1.5e3)
kernel = golovin 
n0 = FT(2^23) 
R0 = FT(30.531e-6) 

#runsettings
num_bins= 128
init_random_seed = 30 
output_steps = [0,1200,2400,3600]
init_method = init_logarithmic
binning_method = mass_density_lnr 



coagsettings = coag_settings{FT}(
    Ns=Ns,
    Δt=Δt,
    ΔV=ΔV,
    golovin_kernel_coeff=golovin_kernel_coeff,
    n0=n0,
    R0=R0,
    kernel=kernel,
)
runsettings=run_settings{FT}(
    num_bins=num_bins,
    init_random_seed=init_random_seed,
    output_steps=output_steps,
    init_method=init_method,
    binning_method=binning_method,
)
X0 = radius_to_volume(settings.R0)
dist = Exponential(X0)

droplets = alpha_sampling(dist,0,settings)


x = binning_func(droplets,runsettings.output_steps[1],runsettings,coagsettings)
plot(x, xlabel="Radius (m)", ylabel="Number of Superdroplets", title="Alpha Sampling Distribution", legend=false)

droppies = init_logarithmic(coagsettings)
x2 = binning_func(droppies,runsettings.output_steps[1],runsettings,coagsettings)
plot!(x2, xlabel="Radius (m)", ylabel="Number of Superdroplets", title="Logarithmic Distribution", legend=false)