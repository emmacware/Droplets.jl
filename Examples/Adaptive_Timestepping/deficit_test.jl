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
seeds = 10
init_names = ["ConstantMultiplicity","Logarithmic","Linear"]



droplets_runs_regular = Dict()
droplets_runs_adaptive = Dict()
all_runs(none(),init_py,log2_Ns,dts,seeds,droplets_runs_regular)
all_runs(Adaptive(),init_py,log2_Ns,dts,seeds,droplets_runs_adaptive)


########################
#Create heatmap matrixes

Droplets_Regular = Dict()
Droplets_Adaptive = Dict()
heatmap_data!(none(),droplets_runs_regular,Droplets_Regular,init_names,log2_Ns,dts,seeds)
heatmap_data!(Adaptive(),droplets_runs_adaptive,Droplets_Adaptive,init_names,log2_Ns,dts,seeds)



#Global Norms
global_droplets_norm_error = mean(mean([droplets_runs_regular[(f,13, 20, i,)][3] for i in 2:seeds]) for f in init_names)
global_droplets_norm_time = mean(mean([droplets_runs_regular[(f,13, 20, i,)][2] for i in 2:seeds]) for f in init_names)

##############################
#HEATMAPS
############################

#No Adaptivity
#Droplets regular
droplets_regular_heatmaps = create_heatmap(Droplets_Regular, global_droplets_norm_error, global_droplets_norm_time,log2_Ns,dts, init_names)
DropletsHeatmaps = plot(droplets_regular_heatmaps..., layout=(length(init_py), 3),suptitle="Droplets No-Adaptive", margin=5Plots.mm)
#PySDM regular
PySDM_Data = JSON.parsefile("/Users/emmaware/PySDM/test_runs_4_27.json")
pysdm_regular_heatmaps = create_heatmap(PySDM_Data["regular"], global_droplets_norm_error, global_droplets_norm_time,log2_Ns,dts, init_names)
PySDMHeatmaps = plot(pysdm_regular_heatmaps..., layout=(length(init_py), 3),suptitle="PySDM No-Adaptive", margin=5Plots.mm)
#Together
plot(DropletsHeatmaps,PySDMHeatmaps,layout=(2,1),size=(1800, 2200))
# savefig("Droplets_no_adaptive.pdf")

##ADAPTIVE
#Droplets adaptive
droplets_adaptive_heatmaps = create_heatmap(Droplets_Adaptive, global_droplets_norm_error, global_droplets_norm_time,log2_Ns,dts, init_names;Deficit = false)
DropletsHeatmapsAdaptive = plot(droplets_adaptive_heatmaps..., layout=(length(init_py), 2), size=(1200, 1200),suptitle="Droplets Adaptive", margin=5Plots.mm)
#PYSDM adaptive
pysdm_adaptive_heatmaps = create_heatmap(PySDM_Data["adaptive"], global_droplets_norm_error, global_droplets_norm_time,log2_Ns,dts, init_names;Deficit = false)
PySDMHeatmapsAdaptive = plot(pysdm_adaptive_heatmaps..., layout=(length(init_py), 2), size=(1200, 1200),suptitle="PySDM Adaptive", margin=5Plots.mm)
#Together
plot(DropletsHeatmapsAdaptive,PySDMHeatmapsAdaptive,layout=(2,1),size=(1200, 1900))
# savefig("Droplets_adaptive.pdf")

#PYPARTMC
pyparterror = JSON.parsefile("/Users/emmaware/Downloads/pypartmc_error_matrix_.json")

PyPartMCErrorHeatmap = each_heatmap(log2_Ns,dts,pyparterror["PyPartMC_err"],"Error",0,10;norm = global_droplets_norm_error,cgrad=:cividis)
plot!(PyPartMCErrorHeatmap,title="PyPartMC Error")

# savefig("PyPartMCError.pdf")
Droplets_Data = JSON.parsefile("DropletsDeficitData.json")

# #write to json
# DropletsData = Dict()
# DropletsData["regular"] = Droplets_Regular
# DropletsData["adaptive"] = Droplets_Adaptive
# open("DropletsDeficitData.json", "w") do io
#     JSON.print(io, DropletsData)
# end