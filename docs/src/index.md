# Droplets

Superdroplet Method for Cloud Microphysics

Features
- Superdroplet coalescence scheme (Shima et al. 2009 SDM), with multiple kernel options
- Droplet activation/diffusional growth via the κ-Köhler equation
- Droplet motion and settling
- MPDATA scheme for fluid advection (Smolarkiewicz, 1998)
- Turbulence closures and functions (Abade and Grabowski., 2018; Mellor-Yamada 2.5)
- 0D and 1D spatial architecture for cloud simulations

Examples (Droplets/Examples)
- Box model collision-coalescence and/or condensation
- Superdroplet microphysics coupling for CloudMicrophysics.jl parcel model
- 1D and 2D MPDATA advection examples
- MPDATA Warm Bubble Boussinesq solver
- 1D Kinematic Rainshaft with prescribed Thermodynamics (Shipway and Hill, 2012, PySDM)
- Single column DyCOMSII-RF02 stratocumulus example with MPDATA advection and a TKE-based turbulence closure (Wyant et al., 2007, Ackerman et al., 2009)
