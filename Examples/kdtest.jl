using Droplets
using StaticArrays
using Distributions
using Random
using Plots

FT = Float64
settings = coag_settings{FT}(Ns=500)
Ns = settings.Ns
ΔV = settings.ΔV
n0 = settings.n0
R0 = settings.R0

ξs = SVector{Ns,FT}(div(n0*ΔV,Ns)*ones(Ns))
X0 = radius_to_volume(R0)# initial volume m3    
Xs = SVector{Ns,FT}(rand(Exponential(X0), Ns))

droplets = static_droplet_attributes(ξs, Xs)
coagdata = coagulation_run{FT}(settings.Ns)

# coalescence_timestep!(Serial(),KiD(), ξs, Xs,coagdata,settings)
plot()
function run(droplets)
    Random.seed!(42)
    runset = run_settings{FT}(scheme=KiD(), coag_threading=Serial())
    bins::Matrix{FT} = zeros(FT, runset.num_bins, length(runset.output_steps))
    threading,scheme = runset.coag_threading, runset.scheme
        coagdata = coagulation_run{FT}(settings.Ns)
        SD_Vol = droplets.X
        SD_Mult = droplets.ξ
            for i  in  1:length(runset.output_steps)
                
                if i !=1
                    timestepper = (runset.output_steps[i]-runset.output_steps[i-1])/settings.Δt
                        for _ in 1:timestepper
                            (;SD_Vol,SD_Mult) = coalescence_timestep!(threading,scheme,SD_Mult,SD_Vol,coagdata,settings)
                        end
                end
                # ξint = SVector{length(SD_Mult), Int}(map(Int, SD_Mult))
                droplets = droplet_attributes(Int.(Vector(SD_Mult)), Vector(SD_Vol))
                x = binning_func(droplets,runset.output_steps[i],runset,settings)
                plot!(x)
                println("Time: ", runset.output_steps[i], " seconds")
            end

    # plot1 = plot!(bins,lc="black",xaxis=:log,ylims=[1e-9,2e14],label=false)
end 

run(droplets)

# x =binning_func(droplets,2.0,runset,settings)