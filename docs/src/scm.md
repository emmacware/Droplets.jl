# Single-Column Model

Eulerian grid state and the routines that advance the single-column model
(condensation, coalescence, turbulence, advection, and diagnostics).

```@autodocs
Modules = [Droplets]
Public = true
Pages = [
    "src/SDfunc/SCM/scm_types.jl",
    "src/SDfunc/SCM/fill_sd_diagnostics.jl",
    "src/SDfunc/SCM/advance.jl",
]
```
