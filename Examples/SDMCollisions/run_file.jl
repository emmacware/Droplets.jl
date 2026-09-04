
using Random
using Plots
using Droplets
using CPUTime
include("../testfunctions.jl")

FT = Float64

#coagsettings
Ns = 2^14
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


drops = runsettings.init_method(coagsettings)


bins,times = coag_runtime(1,drops,coagsettings,runsettings)
plot()
plot_dsd(bins*kg_to_g,runsettings)


#Save plot
savefig("Shima09.pdf")