# 
A Julia language implementation of the superdroplet method [(Shima et al., 2009)](https://doi.org/10.1002/qj.441).
[![Build Status](https://github.com/emmacware/Droplets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmacware/Droplets.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://codecov.io/gh/emmacware/Droplets.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/emmacware/Droplets.jl)    
![alt text](JuliaSDM.svg)

## Features

- Superdroplet coalescence scheme (Shima et al. 2009 SDM), with multiple kernel options
- Droplet activation/diffusional growth via the κ-Köhler equation
- Droplet motion and settling
- MPDATA scheme for fluid advection (Smolarkiewicz, 1998)
- Turbulence closures and functions (Abade and Grabowski., 2018; Mellor-Yamada 2.5)
- 0D and 1D spatial architecture for cloud simulations

## Examples

- Box model collision-coalescence and/or condensation
- Superdroplet microphysics coupling for CloudMicrophysics.jl parcel model
- 1D and 2D MPDATA advection examples
- MPDATA Warm Bubble Boussinesq solver
- 1D Kinematic Rainshaft with prescribed Thermodynamics (Shipway and Hill, 2012, PySDM)
- A single-column DyCOMSII-RF02 stratocumulus example with MPDATA advection and a TKE-based turbulence closure (Wyant et al., 2007, Ackerman et al., 2009)

![alt Text](/Examples/sediment.gif)

## Installation

```julia
using Pkg
Pkg.add("Droplets")
```
# Quick Example: SDM collisions

Example of setting up and running collisions, the [Shima et al., 2009](https://doi.org/10.1002/qj.441) box model collision-coalecence case 

```julia
using Droplets

FT = Float64

Shima2009box_settings = (
    Ns = 2^14,
    Δt = FT(1.0),
    ΔV = FT(1e6),
    kernel = golovin,
    n0 = FT(2^23),
    R0 = FT(30.531e-6),
)
coagsettings = coag_settings{FT}(Shima09_settings...)

#initialize droplets and helper struct
droplets = init_ξ_const(coagsettings) 
allocations = coagulation_run{FT}(coagsettings.Ns)

timesteps = 20
for _ in 1:timesteps
    # droplet volume (droplets.X) and multiplicity (droplets.ξ) update in place
    coalescence_timestep!(Serial(),none(),droplets,allocations,coagsettings)
end     
```
for more settings and plotting enabled, use the box model example in Droplets/Examples/SDMCollisions/run_file.jl 

or on Colab in a jupyter notebook (making sure to change runtime type to Julia):
[![launch on Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/emmacware/Droplets.jl/blob/main/Examples/SDMCollisions/box_collision_coalescence.ipynb)

## Documentation

- **[Droplets Documentation](https://emmacware.github.io/Droplets.jl/dev/)** 

Help for the Droplets functions and structs can searched with 

```julia
julia> ?
help>
```

