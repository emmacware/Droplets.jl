#---------------------------------------------------------
#Collision Coalescence Module
#---------------------------------------------------------
#currently, this collision-coalescence is only set up to handle, 
#multiplicity, volume, superdroplet attributes

#---------------------------------------------------------

using Distributions
using Random
using Interpolations

include("eq_types.jl")
include("sampling.jl")
include("kernels.jl")
include("SDM.jl")
include("golovin_analytic_soln.jl")
