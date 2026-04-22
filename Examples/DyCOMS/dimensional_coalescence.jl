
using Random
using Plots
using Droplets
using CPUTime
include("../testfunctions.jl")

FT = Float64

#coagsettings
Nz = 10
Ns_per_grid = 2^10 +1 
Ns = Ns_per_grid * Nz
Δt = FT(1.0)
ΔV = FT(1e6)
golovin_kernel_coeff = FT(1.5e3)
kernel = golovin 
n0 = FT(2^23) 
R0 = FT(30.531e-6) 
coagsettings = coag_settings{FT}(Ns=Ns_per_grid,ΔV=ΔV, n0=n0,Δt=Δt, golovin_kernel_coeff=golovin_kernel_coeff, kernel=kernel, R0=R0)

#runsettings
num_bins= 128
init_random_seed = 30 
output_steps = [0,1200,2400,3600]
init_method = init_ξ_const
binning_method = mass_density_lnr 
runsettings = run_settings{FT}(num_bins=num_bins, init_random_seed=init_random_seed, output_steps=output_steps, init_method=init_method, binning_method=binning_method)



grid_range = [nothing for i in 1:Nz]
droplet_1dtest = droplet_attributes_1d{FT}(zeros(Int,Ns), zeros(FT,Ns),zeros(FT,Ns),zeros(FT,Ns), zeros(Int,Ns),grid_range, collect(1:Ns))
# droplets = droplet_attributes_1d{FT}(ξstart, volumes,dry_r3,z_loc, cell_id,grid_range, dropidx)

for g in 1:Nz
    griddrops = init_method(coagsettings)
    start_idx = (g-1)*Ns_per_grid + 1
    end_idx = g*Ns_per_grid
    droplet_1dtest.ξ[start_idx:end_idx] .= griddrops.ξ
    droplet_1dtest.X[start_idx:end_idx] .= griddrops.X
    droplet_1dtest.dry_r3[start_idx:end_idx] .= 1
    droplet_1dtest.z_loc[start_idx:end_idx] .= g
    droplet_1dtest.cell_id[start_idx:end_idx] .= g
    droplet_1dtest.grid_range[g] = start_idx:end_idx
end


# coalescence_timestep!(scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings)
coagsettings = coag_settings{FT}(Ns=Ns,ΔV=ΔV, n0=n0,Δt=Δt, golovin_kernel_coeff=golovin_kernel_coeff, kernel=kernel, R0=R0)
coagdata = coagulation_run_spatial{FT}(Nz, coagsettings.Ns,droplet_1dtest)

function coag_runtime(randseed::Int,droplets::droplet_attributes,
    coag_settings::coag_settings{FT},run_settings::run_settings{FT},coag_data::coagulation_run_spatial{FT}) where FT<:AbstractFloat
    
    Random.seed!(randseed)
    println("Running simulation...")

    coal_func_time::FT = 0.0
    bins::Matrix{FT} = zeros(FT, run_settings.num_bins, length(run_settings.output_steps))
    threading,scheme = run_settings.coag_threading, run_settings.scheme
    simtime::FT = @CPUelapsed begin
        for i  in  1:length(run_settings.output_steps)
            # if i,seconds in enumerate(run_settings.output_steps)
    
            if i !=1
                timestepper = Int(round((run_settings.output_steps[i]-run_settings.output_steps[i-1])/coag_settings.Δt))
                ctime::FT = @CPUelapsed begin
                    for _ in 1:timestepper
                        coalescence_timestep!(threading,scheme,droplets,coag_data,coag_settings)
                    end
                end
                coal_func_time += ctime
            end
            bins[:,i] = binning_func(droplets,run_settings.output_steps[i],run_settings,coag_settings)
            # println("Time: ", run_settings.output_steps[i], " seconds")
        end
    end
    println("simtime =", simtime)
    println("coal_func_time =", coal_func_time)

    return bins, coal_func_time
end

bins, coal_func_time = coag_runtime(runsettings.init_random_seed,droplet_1dtest,coagsettings,runsettings,coagdata)

# plot1 = plot()
plot1 = plot_dsd(bins*kg_to_g,runsettings)
plot(plot1)

radius_bins_edges = runsettings.radius_bins_edges
mids = 0.5*(radius_bins_edges[1:end-1] + radius_bins_edges[2:end])*1e6
plot1 = plot(mids,bins*kg_to_g,xaxis=:log,label= runsettings.output_steps',legendtitle="Time Steps (s)")


