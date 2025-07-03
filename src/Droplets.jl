module Droplets

include(joinpath("SDfunc", "constants.jl"))
include(joinpath("SDfunc", "droplet_types.jl"))
include(joinpath("SDfunc","Coalescence", "coalescence.jl"))
include(joinpath("SDfunc", "condensation.jl"))
include(joinpath("SDfunc", "setup.jl"))
include(joinpath("SDfunc", "updateposition.jl"))
include(joinpath("SDfunc", "density.jl"))
include(joinpath("SDfunc", "binning.jl"))

end
