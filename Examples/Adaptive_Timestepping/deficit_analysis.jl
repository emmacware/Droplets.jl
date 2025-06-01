##################################################
# Coalescence example
##################################################
using Random
using Combinatorics
using Distributions
using Plots
using Interpolations
using CPUTime
using Droplets
using JSON
using PyCall
include("../PySDM/pysce.jl")
include("/Users/emmaware/.julia/dev/Droplets/Examples/Adaptive_Timestepping/deficit_funcs.jl")


# ##################################################
# # Coalescence example, with vectorized structur
# ##################################################
# #Settings 
# #-------------------------------------------------

FT = Float32

##Setup Analytic Solution for comparison
rset=run_settings{FT}()
cset = coag_settings{FT}()
X0 = radius_to_volume(cset.R0)
Gsolution = density_lnr_golovin_analytic(rset,cset)


#Set UP PYSDM for droplet initialization
si = pyimport("PySDM.physics").si
ConstantMultiplicity = pyimport("PySDM.initialisation.sampling.spectral_sampling").ConstantMultiplicity
Logarithmic = pyimport("PySDM.initialisation.sampling.spectral_sampling").Logarithmic
Linear = pyimport("PySDM.initialisation.sampling.spectral_sampling").Linear
Exponentialpy = pyimport("PySDM.initialisation.spectra").Exponential
initial_spectrum = Exponentialpy(norm_factor=cset.n0*cset.ΔV, scale=(X0))
init_py = (ConstantMultiplicity=ConstantMultiplicity,Logarithmic=Logarithmic,Linear=Linear)


#Run Droplets
log2_Ns = [13,14,15,16,17,18,19]
dts = [20.0,10.0,5.0,2.0,1]
seeds = 7
init_names = ["ConstantMultiplicity","Logarithmic","Linear"]


########################################


Droplets_Data = JSON.parsefile("DropletsDeficitData.json")
PySDM_Data = JSON.parsefile("/Users/emmaware/PySDM/test_runs_4_28.json")
#Global Norms
global_droplets_norm_error = 5e-5
global_droplets_norm_time = 2e-2

#No Adaptivity
#Droplets regular
droplets_regular_heatmaps = create_heatmap(Droplets_Data["regular"], global_droplets_norm_error, global_droplets_norm_time,log2_Ns,dts, init_names)
DropletsHeatmaps = plot(droplets_regular_heatmaps..., layout=(length(init_py),1),suptitle="Droplets No-Adaptive")
#PySDM regular

pysdm_regular_heatmaps = create_heatmap(PySDM_Data["regular"], global_droplets_norm_error, global_droplets_norm_time,log2_Ns,dts, init_names)
PySDMHeatmaps = plot(pysdm_regular_heatmaps..., layout=(length(init_py), 1),suptitle="PySDM No-Adaptive")
#Together
plot(DropletsHeatmaps,PySDMHeatmaps,layout=(2,1),size=(1200, 2400), margin=6Plots.mm)
# savefig("Droplets_no_adaptive.pdf")

##ADAPTIVE
#Droplets adaptive
droplets_adaptive_heatmaps = create_heatmap(Droplets_Data["adaptive"], global_droplets_norm_error, global_droplets_norm_time,log2_Ns,dts, init_names;Deficit = false)
DropletsHeatmapsAdaptive = plot(droplets_adaptive_heatmaps..., layout=(length(init_py), 2), size=(1200, 1200),suptitle="Droplets Adaptive", margin=5Plots.mm)
#PYSDM adaptive
pysdm_adaptive_heatmaps = create_heatmap(PySDM_Data["adaptive"], global_droplets_norm_error, global_droplets_norm_time,log2_Ns,dts, init_names;Deficit = false)
PySDMHeatmapsAdaptive = plot(pysdm_adaptive_heatmaps..., layout=(length(init_py), 2), size=(1200, 1200),suptitle="PySDM Adaptive", margin=5Plots.mm)
#Together
plot(DropletsHeatmapsAdaptive,PySDMHeatmapsAdaptive,layout=(2,1),size=(1200, 1900))
# savefig("Droplets_adaptive.pdf")




plot(title="No Adaptivivity")
p1 = error_versus_time_vanishingpoint(Droplets_Data["regular"]["ConstantMultiplicity"],log2_Ns,dts,:circle, line=:dash)
plot(title="Adaptivivity")
p2 = error_versus_time_vanishingpoint(Droplets_Data["adaptive"]["ConstantMultiplicity"],log2_Ns,dts,:square, line=:solid)
plot(p1,p2,layout=(2,1),size=(500,800),legend=:topright)

plot(title="Droplets No Adaptivivity")
d1 = error_versusdt(Droplets_Data["regular"]["ConstantMultiplicity"]["Error"],log2_Ns,dts,:circle)
plot(title="Droplets Adaptivivity")
d2 = error_versusdt(Droplets_Data["adaptive"]["ConstantMultiplicity"]["Error"],log2_Ns,dts,:square)
plot(d1,d2,layout=(1,2),size=(900,400), legend=:topleft)

plot(title = "PySDM No Adaptivivity")
p1 = error_versusdt(PySDM_Data["regular"]["ConstantMultiplicity"]["Error"],log2_Ns,dts,:circle)
plot(title = "PySDM Adaptivivity")
p2 = error_versusdt(PySDM_Data["adaptive"]["ConstantMultiplicity"]["Error"],log2_Ns,dts,:square)
plot(p1,p2,layout=(1,2),legend=:topright)

plot(title="PyPartMC Error")
pp1 = error_versusdt(pyparterror["PyPartMC_err"],log2_Ns,dts,:circle)


# fplot = Array{Plots.Plot{Plots.GRBackend}}(undef, length(init_names))
# for i in 1:3
#     plot()
#     f = init_names[i]
#     plot(title="Droplets No Adaptivivity")
#     d1 = error_versusdt(Droplets_Regular[f]["Error"],log2_Ns,dts,:circle)
#     plot(title="Droplets Adaptivivity")
#     d2 = error_versusdt(Droplets_Adaptive[f]["Error"],log2_Ns,dts,:square)
#     # plot(d1,d2,layout=(1,2),size=(900,400))

#     plot(title = "PySDM No Adaptivivity")
#     p1 = error_versusdt(PySDM_Data["regular"][f]["Error"],log2_Ns,dts,:circle)
#     plot(title = "PySDM Adaptivivity")
#     p2 = error_versusdt(PySDM_Data["adaptive"][f]["Error"],log2_Ns,dts,:square)
#     # plot(p1,p2,layout=(1,2),size=(900,400))

#     # plot(title="PyPartMC Error")
#     # pp1 = error_versusdt(pyparterror["PyPartMC_err"],log2_Ns,dts,:circle)

#     emptyplot = plot()

#     fplot[i] = plot(d1,p1,d2,p2,layout=(2,2),size=(1200,800),suptitle=f)
# end

# plot(fplot...,layout=(1,3),size=(1800, 600))
 