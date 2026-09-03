
using Random
using Distributions
using Plots
using PyCall
using Droplets
using CPUTime
include("../testfunctions.jl")
include("pysce.jl")


try
    pyimport("PySDM")
catch err
    println("PySDM is required to run this comparison script, pip install PySDM.")
    rethrow(err)
end

# A comparison with PySDM and Droplets.jl for box collision-coalescence

#-------------------------
# Variable Initialization
#-------------------------
n_sd = 2^15 # number of superdroplets
si = pyimport("PySDM.physics").si
rho_w = pyimport("PySDM.physics.constants_defaults").rho_w # kg/m3
ρ = rho_w # kg/m3.
n0 = Integer(2^23) # initial number density of droplets, 1/m^3
Δt = Float64(1.0) # seconds
ΔV = Float64(10^6) # m^3 box volume
R0 = Float64(30.531e-6) # m
X0 = radius_to_volume(R0)#Float64(4*π/3*R0^3) # m^3
Random.seed!(30)
kernel = golovin #hydrodynamic,golovin -- for Julia version, edit PySDM kernel on line 58

#-------------------------

#Set UP PYSDM
ConstantMultiplicity = pyimport("PySDM.initialisation.sampling.spectral_sampling").ConstantMultiplicity
Exponentialpy = pyimport("PySDM.initialisation.spectra").Exponential
initial_spectrum = Exponentialpy(norm_factor=n0*ΔV, scale=(X0*1e18) * si.um^3)
attributes = Dict()
attributes["volume"], attributes["multiplicity"] = ConstantMultiplicity(spectrum=initial_spectrum).sample(n_sd)
attributes["multiplicity"] = round.(attributes["multiplicity"])

# attributes["multiplicity"] .= 2.56e8
#-------------------------
# Copy for Julia arguments
X_start = attributes["volume"]
ξ_start = Int.(attributes["multiplicity"])


#-------------------------
# Build Test Environment for PySDM
Builder = pyimport("PySDM.builder").Builder
Box = pyimport("PySDM.environments").Box
Coalescence = pyimport("PySDM.dynamics").Coalescence
Golovin = pyimport("PySDM.dynamics.collisions.collision_kernels").Golovin
CPU = pyimport("PySDM.backends").CPU
ParticleVolumeVersusRadiusLogarithmSpectrum = pyimport("PySDM.products").ParticleVolumeVersusRadiusLogarithmSpectrum
radius_bins_edges = 10 .^ range(log10(10*si.um), log10(5e3*si.um), length=128) 

env = Box(dt=Δt * si.s, dv=ΔV * si.m^3)
builder = Builder(n_sd=n_sd, backend=CPU()) #, environment=env)
builder.set_environment(env)
builder.add_dynamic(Coalescence(collision_kernel=Golovin(b=1.5e3 / si.s)))
products = [ParticleVolumeVersusRadiusLogarithmSpectrum(radius_bins_edges=radius_bins_edges, name="dv/dlnr")] 
particulator = builder.build(attributes, products)


#-------------------------
# Run
#-------------------------
FT = Float64
coagsettings = coag_settings{FT}(Ns=n_sd,Δt=Δt,ΔV=ΔV,)
runsettings=run_settings{FT}(coag_threading =Serial(),num_bins=127, output_steps=collect(0:1200:3600))

plot()
#analytic solution

asol = analytic_soln([0.001,1200,2400,3600],radius_bins_edges,n0,X0,ΔV, 1500)
plot_dsd(asol*kg_to_g,runsettings,color="red",label=["Analytic Soln" false false false],legend=true)

#julia
drops = droplet_attributes{FT}(ξ_start,X_start)
bins,times = coag_runtime(1,drops,coagsettings,runsettings)
plot2 = plot_dsd(bins*kg_to_g,runsettings,color="black",label=["Droplets.jl" false false false],legend=true)

#pysdm
pytime = @elapsed begin
for step = 0:1200:3600
    println("pystep: ", step)
    particulator.run(step - particulator.n_steps)

    scope = 2
    water = copy(particulator.products["dv/dlnr"].get()[:] * rho_w)

    if step != 0
    new_water = copy(water)
    for j in 1:2
        # scope = 2
        for i in scope+1:length(water)-scope
            new_water[i] = mean(water[i-scope:i+scope])
        end
        scope = 1
        for i in scope+1:length(water)-scope
            water[i] = mean(new_water[i-scope:i+scope])
        end
    end
    end

    if step == 0
        label = "PySDM"
    else
        label = false
    end

    plot2 = plot!(
        (radius_bins_edges[1:end-1]+radius_bins_edges[2:end])/2 / si.um, water*kg_to_g,
        lc=:dodgerblue,xaxis=:log, label= label) #"t = $step s")   
end
end
println("PySDM time: ", pytime)
display(plot2)

#-------------------------
# Plot time stamps
annotate!(13, 1, text("t= 0s", 8, :left))
annotate!(70, 1, text("t= 1200s", 8, :left))
annotate!(270, 1, text("t= 2400s", 8, :left))
annotate!(1100, 1, text("t= 3600s", 8, :left))
title!("Comparison of Droplets.jl and PySDM")
xlabel!("Particle radius [µm]")
ylabel!("Mass Density [g/m^3/(unit dr/r)]")
savefig("box_collision_coalescence_droplets_pysdm.png")






